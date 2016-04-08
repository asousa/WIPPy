# calc_scattering.py
#
# Basically a port of readOneFileInC -- this is the bulk of the wipp code.
# 
#   1. Load ray files
#   2. Interpolate over space and frequency to get a realistic grid
#   3. Find crossings of the given field line
#   4. Determine scattering at resonance

import pandas as pd
import numpy as np 
from load_rayfile import load_rayfile
import sim_consts as sc
from ionoAbsorp import ionoAbsorp
import matplotlib.pyplot as plt
from scipy import interpolate
import itertools



def calc_scattering(directory, I0, center_lat, lower_freq, upper_freq, L_shells, dlat):
    ''' The main event! '''
    # Load rayfiles:
    ray_lf = load_rayfile(directory, frequency=lower_freq)
    ray_hf = load_rayfile(directory, frequency=upper_freq)

    dfreq  = abs(upper_freq - lower_freq)

    assert ray_lf.keys() == ray_hf.keys(), "Key mismatch between ray files"


    in_lats = sorted(ray_lf.keys())


    # Generate fine-scale interpolating factors
    lat_fine_grid = np.linspace(0, 1, sc.DIV_LAT_NUM)
    freq_fine_grid= np.linspace(0, 1, sc.DIV_FREQ_NUM)
    interp_grid = []
    for l in lat_fine_grid:
        for f in freq_fine_grid:
            interp_grid.append([l, f])
    
    interp_grid = np.array(interp_grid)
    BL_fact = 1 - interp_grid[:,0]-interp_grid[:,1] + interp_grid[:,0]*interp_grid[:,1]
    TL_fact = interp_grid[:,0] - interp_grid[:,0]*interp_grid[:,1]
    BH_fact = interp_grid[:,1] - interp_grid[:,0]*interp_grid[:,1]
    TH_fact = interp_grid[:,0]*interp_grid[:,1]
    
    #print BL_fact
    #print interp_grid
    
    out_data = []
    for x in xrange(len(in_lats)-1):

        print "Ray starting at %g degrees"%in_lats[x]
        # The four rays over which we'll interpolate: 
        # (Two adjacent latitudes, two adjacent frequencies)
        BL = ray_lf[in_lats[x]]   # bottom lat, low frequency
        BH = ray_hf[in_lats[x]]   # bottom lat, high frequency
        TL = ray_lf[in_lats[x+1]] # top lat, low frequency
        TH = ray_hf[in_lats[x+1]] # top lat, high frequency


        # Find peak power in the simulation -- assume 4kHz, 
        # slightly offset to avoid central null
        MAX_POWER = lightning_power(I0, center_lat, dlat, dfreq, 4000, center_lat, 0.7)

        # (initPwr)
        BL.power = BL.power*lightning_power(I0, center_lat, dlat, dfreq, BL.frequency, BL.lat, 0.7);
        BH.power = BH.power*lightning_power(I0, center_lat, dlat, dfreq, BH.frequency, BH.lat, 0.7);
        TL.power = TL.power*lightning_power(I0, center_lat, dlat, dfreq, TL.frequency, TL.lat, 0.7);
        TH.power = TH.power*lightning_power(I0, center_lat, dlat, dfreq, TH.frequency, TH.lat, 0.7);

        BL.scaled = True
        BH.scaled = True
        TL.scaled = True
        TH.scaled = True

        # Interpolate onto uniform time grid (first half of doInterp)
        t = np.linspace(0, sc.T_MAX, sc.NUM_STEPS)  # New sampling grid
        
        #BL_unif = interp_dataframe(BL, t, 'TG',['LAT','power','ELE'])
        BL_unif = interp_dataframe(BL, t, 'TG')
        BH_unif = interp_dataframe(BH, t, 'TG')
        TL_unif = interp_dataframe(TL, t, 'TG')
        TH_unif = interp_dataframe(TH, t, 'TG')

        # Next: Interpolate over frequency, check for crossings, etc
        # Mask off any values for which all four rays are more than L_MARGIN
        # away from the target field line
        # (Not outside --> coarse-grained points that may have crossings)
        
        mask = ~outside_L(BL_unif, BH_unif, TL_unif, TH_unif, L_shells)

        # Also get one timestep previous - so we can tell the direction of the segment
        mask_prev = np.roll(mask,-1)    
        
        # Interpolate yet again (just the hits only) over the fine-scale spatial grid:
        if any(mask):
            print "testing %g cases (coarse grid)"%(np.sum(mask))
            # Current timestep of ray, interpolated onto fine grid
            l_int       = (np.outer(BL_unif[mask].ELE, BL_fact) + 
                           np.outer(TL_unif[mask].ELE, TL_fact) + 
                           np.outer(BH_unif[mask].ELE, BH_fact) + 
                           np.outer(TH_unif[mask].ELE, TH_fact)).ravel()
            
            lat_int     = (np.outer(BL_unif[mask].LAT, BL_fact) +
                           np.outer(TL_unif[mask].LAT, TL_fact) + 
                           np.outer(BH_unif[mask].LAT, BH_fact) + 
                           np.outer(TH_unif[mask].LAT, TH_fact)).ravel()

            # Previous timestep of ray, interpolated onto fine grid
            l_int_prev  = (np.outer(BL_unif[mask_prev].ELE, BL_fact) +
                           np.outer(TL_unif[mask_prev].ELE, TL_fact) + 
                           np.outer(BH_unif[mask_prev].ELE, BH_fact) +
                           np.outer(TH_unif[mask_prev].ELE, TH_fact)).ravel()
            
            lat_int_prev = (np.outer(BL_unif[mask_prev].LAT, BL_fact) +
                            np.outer(TL_unif[mask_prev].LAT, TL_fact) + 
                            np.outer(BH_unif[mask_prev].LAT, BH_fact) + 
                            np.outer(TH_unif[mask_prev].LAT, TH_fact)).ravel()

            # Now we have lists of fine-scale interpolated L-shell and latitudes
            # for two timesteps, at each instance where we might have some crossings.
            # --> Go through the (l, lat, l_prev, lat_prev) lists and check for
            #     overlap with any of the simulation EA segments.
            
            t_inds, ea_inds = check_crossings(l_int, lat_int, l_int_prev, lat_int_prev, L_shells)
            

            #out_data.append(cross_inds)
    
    #return out_data

def check_crossings(l, lat, l_prev_vec, lat_prev_vec, L_target):
    ''' Finds line segments which cross L_target;
        returns indices of the entries to interpolate over.
        inputs: l, lat: L-shell and latitude arrays at current timestep
                l_prev, lat_prev: L-shell and latitude arrays at previous timestep
                L_target: L-shell to check crossings at
    '''

    # Generate the equal-area array
    EA_array = gen_EA_array(L_target)

    # print np.shape(l)
    # print np.shape(l_prev_vec)
    # print np.shape(lat)
    # print np.shape(lat_prev_vec)
    # Step through the l and lat vectors

    out_list = dict()


    # Vector l, no pre-masking
    r1ray = l_prev_vec * pow( np.cos(lat_prev_vec*sc.D2R) , 2 )
    x1ray = r1ray * np.cos(lat_prev_vec*sc.D2R)
    y1ray = r1ray * np.sin(lat_prev_vec*sc.D2R)

    r2ray = l * pow( np.cos(lat*sc.D2R) , 2 )
    x2ray = r1ray * np.cos(lat*sc.D2R)
    y2ray = r1ray * np.sin(lat*sc.D2R)

    Aray = y1ray - y2ray
    Bray = x2ray - x1ray
    Cray = x1ray*y2ray - y1ray*x2ray

    #print np.shape(Aray), np.shape(Bray), np.shape(Cray)

    #tmp = np.outer(Aray,EA_array['x1'])
    #print np.shape(tmp)
    val1 = np.outer(Aray,EA_array['x1']) + np.outer(Bray,EA_array['y1']) + Cray[:,np.newaxis]
    val2 = np.outer(Aray,EA_array['x2']) + np.outer(Bray,EA_array['y1']) + Cray[:,np.newaxis]

    val3 = np.outer(x1ray,EA_array['EAa']) + np.outer(y1ray,EA_array['EAb']) + EA_array['EAc'][np.newaxis,:]
    val4 = np.outer(x2ray,EA_array['EAa']) + np.outer(y2ray,EA_array['EAb']) + EA_array['EAc'][np.newaxis,:]

    mask = (val1*val2 <= 0) & (val3*val4 <= 0)

    # Rows and columns where mask is true:
    mr, mc = np.where(mask)

    # mr: index in time series
    # mc: index in EA array
    return mr, mc

    #outs = zip(mr,mc)
    # row: index in time series
    # col: index in EA

    # # Same format as before (with Cartesian coordinates for the ray segments)
    # out_inds = zip(mr,mc)
    # out_cartesian = dict()
    # for o in out_inds:
    #     #print o[0]
    #     out_cartesian[o[0]] = [(x1ray[o[0]], y1ray[o[0]]),(x2ray[o[0]],y2ray[o[0]])]

    # return out_cartesian

# ------------------------------
#     Loop version (slower, but easier to read)


    # for x in xrange(1, len(l)):
    #     #print x
    #     lat_prev = lat_prev_vec[x]
    #     lat_curr = lat[x]
    #     l_prev   = l_prev_vec[x]
    #     l_curr   = l[x]

    #     # # Order the limits (ray could be going either direction):
    #     # # if (lat_prev < lat_curr):
    #     # #     lat_low = lat_prev
    #     # #     lat_high= lat_curr
    #     # # else:
    #     # #     lat_low = lat_curr
    #     # #     lat_high= lat_prev
    #     # lat_low = min(lat_prev, lat_curr)
    #     # lat_high= max(lat_prev, lat_curr)

    #     # Find index of any EA segments the ray may have passed through:
    #     # Get index numbers into the EA array    
    #     # EA_iLow  = int(np.ceil(  (lat_low  - sc.EALimS) / sc.EAIncr ));
    #     # EA_iHigh = int(np.floor( (lat_high - sc.EALimS) / sc.EAIncr ));

    #     # if(   ( EA_iLow > EA_iHigh ) |
    #     #       ( EA_iLow < 0 )        | 
    #     #       ( EA_iHigh > (sc.EALimN-sc.EALimS)/sc.EAIncr )  ):
    #     #     pass
    #     # else:
    #         # Array of EA index numbers to check for crossings at:
    #         # EA_i = np.arange(EA_iLow, EA_iHigh + 1)

    #     # I'm just clipping this math straight from Jacob. Yikes!

    #     # Coordinates of the ray segment start and endpoints
    #     # (radial and Cartesian, in L-shells)
    #     r1ray = l_prev * pow( np.cos(lat_prev*sc.D2R) , 2 )
    #     x1ray = r1ray * np.cos(lat_prev*sc.D2R)
    #     y1ray = r1ray * np.sin(lat_prev*sc.D2R)

    #     r2ray = l_curr * pow( np.cos(lat_curr*sc.D2R) , 2 )
    #     x2ray = r1ray * np.cos(lat_curr*sc.D2R)
    #     y2ray = r1ray * np.sin(lat_curr*sc.D2R)

    #     Aray = y1ray - y2ray
    #     Bray = x2ray - x1ray
    #     Cray = x1ray*y2ray - y1ray*x2ray

    #     # Formal check for genuine crossings:        
    #     # val1 = Aray*EA_array.iloc[EA_i]['x1'] + Bray*EA_array.iloc[EA_i]['y1'] + Cray
    #     # val2 = Aray*EA_array.iloc[EA_i]['x2'] + Bray*EA_array.iloc[EA_i]['y1'] + Cray

    #     # val3 = x1ray*EA_array.iloc[EA_i]['EAa'] + y1ray*EA_array.iloc[EA_i]['EAb'] + EA_array.iloc[EA_i]['EAc']
    #     # val4 = x2ray*EA_array.iloc[EA_i]['EAa'] + y2ray*EA_array.iloc[EA_i]['EAb'] + EA_array.iloc[EA_i]['EAc']

    #     # Masking out a dataframe with .iloc takes twice as long as just doing all the multiplies!
    #     val1 = Aray*EA_array['x1'] + Bray*EA_array['y1'] + Cray
    #     val2 = Aray*EA_array['x2'] + Bray*EA_array['y1'] + Cray

    #     val3 = x1ray*EA_array['EAa'] + y1ray*EA_array['EAb'] + EA_array['EAc']
    #     val4 = x2ray*EA_array['EAa'] + y2ray*EA_array['EAb'] + EA_array['EAc']

    #     mask = (val1*val2 <= 0) & (val3*val4 <= 0)

    #     #print np.shape(mask)

    #     if any(mask):
    #         out_list[x] = [(x1ray,y1ray),(x2ray,y2ray)]
    #         #out_list[x] = np.where(mask)

    # return out_list


def outside_L(BL, BH, TL, TH, L_target):
    ''' Checks whether or not L_target is outside the volume bounded by
        the four rays. Returns true if outside, false if potentially inside.
        L_MARGIN is defined in sim_consts.py
    '''
    lowBound = L_target - sc.L_MARGIN
    upBound  = L_target + sc.L_MARGIN
    #print "lowBound: %g, upBound: %g"%(lowBound, upBound)

    # Bitwise operators (&, | vs and, or) yield boolean vectors
    return ((   (BL.ELE <= lowBound) &  
                (TL.ELE <= lowBound) &  
                (BH.ELE <= lowBound) &
                (TH.ELE <= lowBound) ) |
            (   (BL.ELE >= upBound)  &
                (TL.ELE >= upBound)  &
                (BH.ELE >= upBound)  &
                (TH.ELE >= upBound)  ))



def interp_dataframe(df, t_new, t_label, cols=None):
    ''' Interpolates the contents of a dataframe over a new axis. 
    Inputs:
        df:      The dataframe
        t_new:   The new axis
        t_label: The string label of the x-axis to interpolate over
        cols:    A list of columns to interpolate over.
    Outputs:
        A new dataframe
    '''
    tmp = pd.DataFrame()
    tmp[t_label] = t_new
 
    if not cols:
        cols = df.columns.values

    for c in cols:
        interpolater = interpolate.interp1d(df[t_label],df[c])
        tmp[c] = interpolater(t_new)
    
    return tmp



def lightning_power(I0, center_lat, dlat, dfreq, frequency, ray_lat, dist_long = 0.7):
    ''' Calculates initial power of the ray at top of ionosphere
        Accounts for ionospheric absorption (currently Helliwell model).
        Inputs:
        I0:          Flash peak current, in kA
        center_lat:  The flash latitude
        frequency:   The ray frequency
        ray_lat:     The latitude of the current ray
        dist_long:   Longitude distance, in degrees
        dlat:        Latitude separation between adjacent rays, in degrees (for geometric spreading)
        dfreq:       Frequency spread between adjacent rays, in Hz
    '''

    # Latitude distance, longitude distance between source and ray
    dist_lat = (sc.R_E + sc.H_IONO/2)*abs(center_lat - ray_lat)*sc.D2R    
    dist_long= (sc.R_E + sc.H_IONO/2)*abs(dist_long)*sc.D2R
    #print "dist_lat: ",dist_lat,"dist_long: ",dist_long

    dist_iono = np.hypot(dist_lat, dist_long)
    dist_tot =  np.hypot(dist_iono, sc.H_IONO)
    #print "Dist_iono: ", dist_iono,"dist_tot: ",dist_tot
    xi = np.arctan2(dist_iono, sc.H_IONO)
    #print "xi: ",xi
    # Power at sub-ionosphere point
    w = 2*np.pi*frequency
    w_squared = pow(w, 2)

    S = 1.0/(sc.Z0) * pow( sc.H_E * I0 * (2e-7)*(np.sin(xi)/dist_tot)* w * (sc.A - sc.B),2)/( (w_squared + pow(sc.A, 2))*(w_squared + pow(sc.B,2)))

    S_vert = S * np.cos(xi)   # Vertical component


    # Get ionosphere attenuation:
    atten_factor = pow(10, -ionoAbsorp(ray_lat, frequency)/10.0)

    # I don't know where Jacob got this attenuation factor from. Geometric spreading? Hmm.
    return S_vert * atten_factor * dlat * sc.D2R * (sc.R_E + sc.H_IONO) * dfreq * 0.87788331




def gen_EA_array(L):
    ''' Returns a dataframe of parameters for the EA (equal area) slices
        on which scattering is to be calculated. A clone of Jacob's initEA_Arr method.

        Inputs: simulation constants (maybe put this as a method input...)
                L_vector: List of L-shells we're computing on
        Outputs: a dataframe.
    '''


    # latitude vector (including endpoints)
    lam_vec = np.arange(sc.EALimS,sc.EALimN + sc.EAIncr, step=sc.EAIncr)
    
    # Trig functions
    clam = np.cos(lam_vec*sc.D2R)
    slam = np.sin(lam_vec*sc.D2R)
    clam2= pow(clam,2)
    slam2= pow(slam,2)
    rootTerm = np.sqrt(1 + 3*slam2)

    # dl_lam: still not sure where he gets this expression from
#    dl_lam = (clam2*clam / rootTerm)*sc.DL0
    dl_lam = (clam2*clam / rootTerm)*sc.L_MARGIN

    x_unit_vect = (3*clam2 - 2) / rootTerm
    y_unit_vect = (3*slam*clam) / rootTerm


    ptR = L*clam2  # Distance to field line (dipole model: r ~ L cos(latitude)^2 )
    ptX = ptR*clam      # Field line in Cartesian coordinates
    ptY = ptR*slam

    # Coordinates of the edges of the window
    x1 = ptX - x_unit_vect*dl_lam
    x2 = ptX + x_unit_vect*dl_lam
    y1 = ptY - y_unit_vect*dl_lam
    y2 = ptY + y_unit_vect*dl_lam

    # Output parameters:
    d = pd.DataFrame(columns=['lam','x1','x2','y1','y2','EAa','EAb','EAc','EA_length'])

    d['lam'] = lam_vec
    d['EAa'] = y1 - y2
    d['EAb'] = x2 - x1
    d['EAc'] = x1*y2 - y1*x2
    d['EA_length'] = np.sqrt(pow(x1-x2,2)+pow(y1-y2,2))*sc.R_E
    d['x1'] = x1
    d['x2'] = x2
    d['y1'] = y1
    d['y2'] = y2

    return d

if __name__ =='__main__':
    directory = '/shared/users/asousa/WIPP/WIPPy/python/'
    
    #calc_scattering(directory=directory, lower_freq=200, upper_freq=240,L_shells=None)
    directory='/shared/users/asousa/WIPP/WIPPy/python/'
    center_lat=35
    lower_freq=200
    upper_freq=240
    L_shells=2.5
    I0 = -100000
    dlat = 1
    calc_scattering(directory, I0, center_lat, lower_freq, upper_freq, L_shells, dlat)
