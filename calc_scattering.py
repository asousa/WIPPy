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
import time



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
    # BL_fact = 1 - interp_grid[:,0]-interp_grid[:,1] + interp_grid[:,0]*interp_grid[:,1]
    # TL_fact = interp_grid[:,0] - interp_grid[:,0]*interp_grid[:,1]
    # BH_fact = interp_grid[:,1] - interp_grid[:,0]*interp_grid[:,1]
    # TH_fact = interp_grid[:,0]*interp_grid[:,1]

    fine_grid_size = sc.DIV_LAT_NUM*sc.DIV_FREQ_NUM

    
    #print BL_fact
    #print interp_grid
    
    # A dataframe of interpolated ray parameters for every EA crossing
    out_data = pd.DataFrame()

    # Execution time
    tstart = time.time()



    # Loop over incident ray latitudes:
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
        BL.power = BL.power*lightning_power(I0, center_lat, dlat, dfreq, BL.frequency, BL.launch_lat, 0.7);
        BH.power = BH.power*lightning_power(I0, center_lat, dlat, dfreq, BH.frequency, BH.launch_lat, 0.7);
        TL.power = TL.power*lightning_power(I0, center_lat, dlat, dfreq, TL.frequency, TL.launch_lat, 0.7);
        TH.power = TH.power*lightning_power(I0, center_lat, dlat, dfreq, TH.frequency, TH.launch_lat, 0.7);

        BL.scaled = True
        BH.scaled = True
        TL.scaled = True
        TH.scaled = True

        # Interpolate onto uniform time grid (first half of doInterp)
        t = np.linspace(0, sc.T_MAX, sc.NUM_STEPS)  # New sampling grid
        
        #BL_unif = interp_dataframe(BL, t, 'TG',['LAT','power','ELE'])

        # Interpolating class -- initially interpolates onto a uniform time grid,
        # then can be called to get interpolated time vectors between the four corner vectors
        R = ray_interpolator(BL, TL, BH, TH, t)


        # Next: Interpolate over frequency, check for crossings, etc
        # Mask off any values for which all four rays are more than L_MARGIN
        # away from the target field line
        # (Not outside --> coarse-grained points that may have crossings)
        
        mask = ~outside_L(R.BL, R.BH, R.TL, R.TH, L_shells)

        # Also get one timestep previous - so we can tell the direction of the segment
        mask_prev = np.roll(mask,-1)    
        
        # Interpolate yet again (just the hits only) over the fine-scale spatial grid:
        if any(mask):
            print "testing %g cases (coarse grid)"%(np.sum(mask))

            l_int        = R.all_vals_at(interp_grid, 'l_sh', mask)
            lat_int      = R.all_vals_at(interp_grid, 'lat', mask)

            l_int_prev   = R.all_vals_at(interp_grid, 'l_sh', mask_prev)
            lat_int_prev = R.all_vals_at(interp_grid, 'lat', mask_prev)

            # Now we have lists of fine-scale interpolated L-shell and latitudes
            # for two timesteps, at each instance where we might have some crossings.
            # --> Go through the (l, lat, l_prev, lat_prev) lists and check for
            #     overlap with any of the simulation EA segments.
            
            cross_inds, ea_inds, cross_coords = check_crossings(l_int, lat_int, l_int_prev, lat_int_prev, L_shells)
            #print cross_inds
            
            # unravel the indices in the cross_inds list, to get time and weight indices
            mask_inds, fine_grid_inds =  np.unravel_index(cross_inds, (sum(mask), fine_grid_size))
            #print mask_inds
            #print fine_grid_inds
            # tmp = np.where(mask)[0]
            
            # print tmp
            # print tmp[mask_inds]

            # Time indexes in full ray -- index of confirmed crossings in the masked-off array
            time_inds = np.where(mask)[0][mask_inds]
            
            #print BL_fact[fine_grid_inds], TL_fact[fine_grid_inds], BH_fact[fine_grid_inds], TH_fact[fine_grid_inds]
            #print R.BL[time_inds]
            #print sum(R.BL['ELE'][mask] - R.BL[mask]['ELE'])
            #print max(time_inds), min(time_inds)
            # To confirm: Let's interpolate again and make sure we've got the right values for l and lat:
            #lcheck   = R.vals_at(interp_grid[fine_grid_inds,:],'ELE',mask=time_inds)
            #latcheck = R.vals_at(interp_grid[fine_grid_inds,:],'LAT',mask=time_inds)

            #print lat_int[cross_inds]
            #print latcheck
            # print lcheck
            # print latcheck

            df_fine = R.dataframe_at(interp_grid[fine_grid_inds,:],time_inds)
            df_fine['cross_coords'] = cross_coords  # Cartesian coordinates of start and endpoints (for plotting)
            df_fine['EA_index'] = ea_inds           # Index of which EA array it intercepts


            #print df_fine
            # print sum(abs(df_fine['ELE'] - l_int[cross_inds]))
            # print sum(abs(df_fine['LAT'] - lat_int[cross_inds]))
            #print df_fine['cross_coords']
            print np.shape(df_fine)

            out_data = pd.concat([out_data, df_fine])



    tstop = time.time()
    print "Elapsed time (Interpolation and crossing detection): %g seconds"%(tstop - tstart) 

    return out_data



def calc_resonant_pitchangle_change(crossing_df, L_target):

    # Start timer
    tstart = time.time()
    # Generate energy and velocity arrays
    E_tot_arr = pow(10,sc.E_EXP_BOT + sc.DE_EXP*np.arange(0,sc.NUM_E))
    v_tot_arr = sc.C*np.sqrt(1 - pow(sc.E_EL/(sc.E_EL + E_tot_arr),2))


    L = L_target

    epsm = (1/L)*(sc.R_E + sc.H_IONO)/sc.R_E
    alpha_eq = np.arcsin(np.sqrt( pow(epsm, 3)/np.sqrt(1 + 3*(1 - epsm))  ))


    # Gen EA array (maybe pass this as an arg)
    EA_array = gen_EA_array(L_target)


    tstop = time.time()

    #print EA_array
    # Loop through EA segments:
    #for EA_ind, EA_row in EA_array.iterrows():
    for EA_ind in np.unique(crossing_df['EA_index']):   
        #print EA_ind

        lat = EA_array['lam'][EA_ind]
        print 'EA segment at latitude = ',lat

        slat = np.sin(lat*sc.D2R)
        clat = np.cos(lat*sc.D2R)
        slat_term = np.sqrt(1 + 3*slat*slat)

        # Magnetic field at current location
        #B_local = sc.B_eq*slat_term/(pow(clat,6))


        # Lower hybrid frequency
        wh = (sc.Q_EL*sc.B0/sc.M_EL)*(1/pow(L,3))*slat_term/pow(clat,6)
        #wh = (2*np.pi*880000) / (pow(L,3)*slat_term/pow(clat,6))       # Jacob used a hard coded value for qB/m
        dwh_ds = (3*wh/(L*sc.R_E)) *(slat/slat_term) * (1/pow(slat_term,2) + 2/pow(clat,2))

        # flight time constants (divide by velocity to get time to ionosphere, in seconds)
        ftc_N, ftc_S = get_flight_time_constant(L, lat, alpha_eq)

        # Local loss cone angle
        alpha_lc = np.arcsin(np.sqrt(slat_term/pow(clat,6) )*np.sin(alpha_eq))
        salph = np.sin(alpha_lc)
        calph = np.cos(alpha_lc)

        ds = L*sc.R_E*slat_term*clat*sc.EAIncr*np.pi
        # # Sometimes Python is pretty and readable, but now is not one of those times
        # for cell_ind, cells in crossing_df[crossing_df['EA_index'] == EA_ind].iterrows():
        #     print np.shape(cell)

        # Mask out just the entries which cross the current EA segment:
        cells = crossing_df[crossing_df['EA_index'] == EA_ind]

        if cells.shape[0] > 0:
            print cells.shape
            #print cells.columns

            # where you left off: time and frequency (line 1800)
            t = cells.tg + sc.DT/2.0
            f = cells.frequency + sc.DF/2.0
            pwr = cells.power/sc.DT       # Jacob divides by num_rays too, but it looks like it's always 1
            psi = cells.psi*sc.D2R       

            # Maybe threshold power here?

            # Misc. parameters and trig
            mu = cells.mu
            stixP = cells.stixP
            stixR = cells.stixR
            stixL = cells.stixL

            spsi = np.sin(psi)
            cpsi = np.cos(psi)
            spsi_sq = pow(spsi,2)
            cpsi_sq = pow(cpsi,2)
            n_x  = mu*abs(spsi)
            n_z  = mu*cpsi
            mu_sq = mu*mu
            w = 2.0 * np.pi * f
            k = w*mu/sc.C
            kx = w*n_x/sc.C
            kz = w*n_z/sc.C
            Y = wh / w

            stixS = (stixR + stixL)/2.0
            stixD = (stixR - stixL)/2.0
            stixA = stixS + (stixP-stixS)*cpsi_sq
            stixB = stixP*stixS + stixR*stixL + (stixP*stixS - stixR*stixL)*cpsi_sq
            stixX = stixP/(stixP - mu_sq*spsi_sq)
          
            rho1=((mu_sq-stixS)*mu_sq*spsi*cpsi)/(stixD*(mu_sq*spsi_sq-stixP))
            rho2 = (mu_sq - stixS) / stixD

            Byw_sq = ( (2.0*sc.MU0/sc.C) * (pwr*stixX*stixX*rho2*rho2*mu*abs(cpsi)) /
                     np.sqrt( pow((np.tan(psi)-rho1*rho2*stixX),2) + pow((1+rho2*rho2*stixX),2)) )

            # RMS wave components
            Byw = np.sqrt(Byw_sq);
            Exw = abs(sc.C*Byw * (stixP - n_x*n_x)/(stixP*n_z))
            Eyw = abs(Exw * stixD/(stixS-mu_sq))
            Ezw = abs(Exw *n_x*n_z / (n_x*n_x - stixP))
            Bxw = abs((Exw *stixD*n_z/sc.C)/ (stixS - mu_sq))
            Bzw = abs((Exw *stixD *n_x) /(sc.C*(stixX - mu_sq)));

            # Oblique integration quantities
            R1 = (Exw + Eyw)/(Bxw+Byw)
            R2 = (Exw - Eyw)/(Bxw-Byw)
            w1 = (sc.Q_EL/(2*sc.M_EL))*(Bxw+Byw)
            w2 = (sc.Q_EL/(2*sc.M_EL))*(Bxw-Byw)
            alpha1 = w2/w1
          

            # Loop at line 1897 -- MRES loop
            # Since we're trying to keep it vectorized:
            # target dimensions are num(cell rows) rows x num(mres) columns.
            # (Sometimes Python is beautiful, but now is not one of those times)
            #  --> Note to future editors: the (vec)[:,np.newaxis] command is similar to repmat or tile:
            #  It will tile the vector along the new axis, with length dictated by whatever operation 
            #  follows it. In Python lingo this is "broadcasting"

            # Resonant modes to sum over: Integers, positive and negative
            mres = np.linspace(-5,5,11)
            
            # Parallel resonance velocity V_z^res  (pg 5 of Bortnik et al)
            t1 = w*w*kz*kz
            t2 = (mres*mres*wh*wh)[np.newaxis,:]-(w*w)[:,np.newaxis]
            t3 = (kz*kz)[:,np.newaxis] + ((mres*wh*mres*wh)[np.newaxis,:]/(pow(sc.C*np.cos(alpha_lc),2)))
            
            direction = np.outer(np.sign(kz),np.sign(mres))
            direction[:,mres==0] = -1*np.sign(kz)[:,np.newaxis] # Sign 

            # v_para_res = (ta - tb)
            v_para_res = ( direction*np.sqrt(abs(t1[:,np.newaxis] + t2*t3)) - (w*kz)[:,np.newaxis] ) / t3;

            # Getting some obnoxious numerical stuff in the case of m=0:
            v_para_res[:,mres==0] = (w/kz)[:,np.newaxis]
            v_tot_res = v_para_res / np.cos(alpha_lc)

            E_res = sc.E_EL*(1.0/np.sqrt( 1-(v_tot_res*v_tot_res/(sc.C*sc.C)) ) -1 )

            # Start and end indices in the energy array (+-20% energy band)
            e_starti = np.floor((np.log10(E_res)-sc.E_EXP_BOT-sc.E_BANDWIDTH)/sc.DE_EXP)
            e_endi   = np.ceil(( np.log10(E_res)-sc.E_EXP_BOT+sc.E_BANDWIDTH)/sc.DE_EXP)

            # Threshold to range of actual indexes (some energies may be outside our target range)
            np.clip(e_starti,0,sc.NUM_E, out=e_starti)
            np.clip(e_endi,  0,sc.NUM_E, out=e_endi)

            #inds = zip(e_starti.ravel(), e_endi.ravel())
            #print np.shape(inds)
            #v_tot = direction*v_tot_arr[e_]

            # Where you left off: Adding in the next loop (V_TOT). Need to iterate from estarti to eendi
            # for each cell in the 2d arrays... hmm.... how to vectorize nicely. Hmm.

    return v_tot_arr








def get_flight_time_constant(L, lat, alpha_eq):
 # * FUNCTION: getFltConst
 # * ---------------------
 # * This function calculates the 'flight-time-const' for a particular 
 # * latitude,  and returns the constant for a particle flying adiabatically
 # * from that latitude, to the Northern hemisphere, and a similar constant 
 # * to the Southern hemisphere.  This constant is then multiplied by 
 # * 1/v_tot for the particular particle, which gives the total flight
 # * time in seconds to the appropriate hemisphere.
 # *
    num_div = 10000   # number of steps for numeric evaluation of integral
    sin_alpha_eq_sq = pow(np.sin(alpha_eq),2)
    endx = np.sin(lat*sc.D2R)
    dx = endx/num_div

    # Evaluate Walt's flight-time constant, from equator to mirror pt.
    # (Equation 4.28 -- page 44 of Walt's book)
    walt_tau = (0.02925/sc.R_E)*sc.C*(1 - 0.4635*pow(sin_alpha_eq_sq,0.375))
    #walt_tau = (0.117/sc.R_E)*sc.C*(1 - 0.4635*pow(np.sin(alpha_eq),0.75))

    # Evaluate flight-time integral from equator to target latitude:
    # (Equation 4.27 -- page 44 of Walt's book)

    x = 0
    I = 0
    for n in xrange(0,num_div):
        xterm = np.sqrt(1 + 3*x*x)
        I += dx*xterm / np.sqrt(1 - (sin_alpha_eq_sq/pow(1-x*x,3))*xterm)
        x += dx

    # Northern flight time is the total flight time minus the offset
    # Soutehrn is total plus offset
    flight_time_N = L*sc.R_E*(walt_tau - I)
    flight_time_S = L*sc.R_E*(walt_tau + I)

    #print flight_time_N/sc.C
    #print flight_time_S/sc.C

    return flight_time_N, flight_time_S

def check_crossings(l, lat, l_prev_vec, lat_prev_vec, L_target):
    ''' Finds line segments which cross L_target;
        returns indices of the entries which cross an equal-area segment
        inputs: l, lat: L-shell and latitude arrays at current timestep
                l_prev, lat_prev: L-shell and latitude arrays at previous timestep
                L_target: L-shell to check crossings at
        outputs: 
        time_index, ea_index: vectors of index values 
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

    # mr: index in l and lat
    # mc: index in EA array
    #return mr, mc

    #outs = zip(mr,mc)
    # row: index in time series
    # col: index in EA

    # # Same format as before (with Cartesian coordinates for the ray segments)
    out_inds = zip(mr,mc)
    # print np.shape(out_inds)
    # print np.shape(x1ray)
    # print np.shape(mr)
    #out_cartesian = zip(x1ray[mr], y1ray[mr], x2ray[mr],y2ray[mr])#, zip(x2ray[mr],y2ray[mr])
    out_cartesian = []
    for o in out_inds:
        out_cartesian.append([(x1ray[o[0]], y1ray[o[0]]),(x2ray[o[0]],y2ray[o[0]])])

    # print x1ray[mr]
    # print x2ray[mr]
    #out_cartesian = x1ray[mr]
    

    return mr, mc, out_cartesian

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
    return ((   (BL['l_sh'] <= lowBound) &  
                (TL['l_sh'] <= lowBound) &  
                (BH['l_sh'] <= lowBound) &
                (TH['l_sh'] <= lowBound) ) |
            (   (BL['l_sh'] >= upBound)  &
                (TL['l_sh'] >= upBound)  &
                (BH['l_sh'] >= upBound)  &
                (TH['l_sh'] >= upBound)  ))



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


    ptR = L*clam2       # Distance to field line (dipole model: r ~ L cos(latitude)^2 )
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


class ray_interpolator(object):
    def __init__(self, BL, TL, BH, TH, t):
        '''Inputs: bottom low, bottom high, top low, top high interpolating factors'''
        # Interpolate dataframes onto a uniform time grid:
        self.BL = interp_dataframe(BL, t, 'tg')
        self.BH = interp_dataframe(BH, t, 'tg')
        self.TL = interp_dataframe(TL, t, 'tg')
        self.TH = interp_dataframe(TH, t, 'tg')

        # Ray metadata
        self.BL.frequency = BL.frequency
        self.BH.frequency = BH.frequency
        self.TL.frequency = TL.frequency
        self.TH.frequency = TH.frequency

    def interpolation_weights(self, interp_grid):
        BL_fact = 1 - interp_grid[:,0]-interp_grid[:,1] + interp_grid[:,0]*interp_grid[:,1]
        TL_fact = interp_grid[:,0] - interp_grid[:,0]*interp_grid[:,1]
        BH_fact = interp_grid[:,1] - interp_grid[:,0]*interp_grid[:,1]
        TH_fact = interp_grid[:,0]*interp_grid[:,1]

        return BL_fact, TL_fact, BH_fact, TH_fact
    def all_vals_at(self, interp_grid, val_name, mask):
        '''Returns interpolated values at each time step, and each factor.
           interp_grid: two-column array -- (t, f) interpolating values, 0 to 1.
           val_name: name of the column to interpolate
           mask: binary mask of time bins to look at '''
        # Interpolating weights:
        # BL_fact = 1 - interp_grid[:,0]-interp_grid[:,1] + interp_grid[:,0]*interp_grid[:,1]
        # TL_fact = interp_grid[:,0] - interp_grid[:,0]*interp_grid[:,1]
        # BH_fact = interp_grid[:,1] - interp_grid[:,0]*interp_grid[:,1]
        # TH_fact = interp_grid[:,0]*interp_grid[:,1]
        BL_fact, TL_fact, BH_fact, TH_fact = self.interpolation_weights(interp_grid)

        return (np.outer(self.BL[val_name][mask], BL_fact) + 
                    np.outer(self.BH[val_name][mask], TL_fact) + 
                    np.outer(self.TL[val_name][mask], BH_fact) + 
                    np.outer(self.TH[val_name][mask], TH_fact)).ravel()

    def vals_at(self, interp_grid, val_name, mask):
        ''' Same thing, but without the outer product. mask and factors must be same length.
           BL, TL, BH, TH: Interpolating weights (vectors).
           val_name: name of the column to interpolate
           mask: binary mask of time bins to look at '''
        # BL_fact = 1 - interp_grid[:,0]-interp_grid[:,1] + interp_grid[:,0]*interp_grid[:,1]
        # TL_fact = interp_grid[:,0] - interp_grid[:,0]*interp_grid[:,1]
        # BH_fact = interp_grid[:,1] - interp_grid[:,0]*interp_grid[:,1]
        # TH_fact = interp_grid[:,0]*interp_grid[:,1]

        BL_fact, TL_fact, BH_fact, TH_fact = self.interpolation_weights(interp_grid)

        return (self.BL[val_name][mask]*BL_fact + 
                self.BH[val_name][mask]*TL_fact + 
                self.TL[val_name][mask]*BH_fact + 
                self.TH[val_name][mask]*TH_fact)

    def dataframe_at(self, interp_grid, mask):

        cols = self.BL.columns.values
        df = pd.DataFrame()

        for c in cols:
            df[c] = self.vals_at(interp_grid, c, mask)

        BL_fact, TL_fact, BH_fact, TH_fact = self.interpolation_weights(interp_grid)

        df['frequency'] = (self.BL.frequency*BL_fact + 
                           self.BH.frequency*TL_fact + 
                           self.TL.frequency*BH_fact + 
                           self.TH.frequency*TH_fact)
        return df

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
