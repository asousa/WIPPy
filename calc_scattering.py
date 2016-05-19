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
from scipy.special import jn      # Bessel function of the 1st kind
from scipy.special import fresnel # Fresnel integrals
import itertools
import time
from line_intersection import check_intersection


# def calc_scattering(directory, I0, center_lat, lower_freq, upper_freq, L_shells, dlat):
def calc_scattering(directory, I0, center_lat, lower_freq, upper_freq, L_shells):

    ''' The main event! '''
    # Load rayfiles:
    ray_lf = load_rayfile(directory, frequency=lower_freq)
    ray_hf = load_rayfile(directory, frequency=upper_freq)

    # dfreq  = abs(upper_freq - lower_freq)
    # sc.DIV_FREQ_NUM = np.ceil(sc.DF/sc.F_STEP + 1)
    # print "DIV_FREQ_NUM:", sc.DIV_FREQ_NUM

    assert ray_lf.keys() == ray_hf.keys(), "Key mismatch between ray files"


    all_lats = np.array(sorted(ray_lf.keys()))

    in_lats = all_lats[(all_lats >= (center_lat - sc.LAT_SPREAD/2.0)) & (all_lats <= (center_lat + sc.LAT_SPREAD/2.0))]

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


        # generate interpolating grid (factors may be different depending on rays loaded)
        dlat = in_lats[x+1] - in_lats[x]
        dfreq  = abs(upper_freq - lower_freq)

        DIV_LAT_NUM = np.ceil(dlat/sc.LAT_STEP + 1)
        DIV_FREQ_NUM = np.ceil(dfreq/sc.F_STEP + 1)
        print "DIV_LAT_NUM:", DIV_LAT_NUM
        print "DIV_FREQ_NUM:", DIV_FREQ_NUM

        lat_fine_grid = np.linspace(0, 1, DIV_LAT_NUM)
        freq_fine_grid= np.linspace(0, 1, DIV_FREQ_NUM)

        print "Latitude interpolating steps:  ", lat_fine_grid
        print "Frequency interpolating steps: ", freq_fine_grid

        interp_grid = np.array([(l,f) for l in lat_fine_grid for f in freq_fine_grid])
        fine_grid_size = DIV_LAT_NUM*DIV_FREQ_NUM


        # Find peak power in the simulation -- assume 4kHz, 
        # slightly offset to avoid central null
        MAX_POWER = lightning_power(I0, center_lat, dlat, dfreq, 4000, center_lat, 0.7)
        print "center_lat: ",center_lat
        print "dlat:", dlat
        print "dfreq:", dfreq
        print "MAX_POWER:", MAX_POWER
        # (initPwr)
        BL.power *=lightning_power(I0, center_lat, dlat, dfreq, BL.frequency, BL.launch_lat, 0.7);
        BH.power *=lightning_power(I0, center_lat, dlat, dfreq, BH.frequency, BH.launch_lat, 0.7);
        TL.power *=lightning_power(I0, center_lat, dlat, dfreq, TL.frequency, TL.launch_lat, 0.7);
        TH.power *=lightning_power(I0, center_lat, dlat, dfreq, TH.frequency, TH.launch_lat, 0.7);

        # Interpolate onto uniform time grid (first half of doInterp)
        t = np.linspace(sc.T_STEP, sc.T_MAX, sc.NUM_STEPS)  # New sampling grid
        # print "T-grid:", t
        
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
            
            # unravel the indices in the cross_inds list, to get time and weight indices
            mask_inds, fine_grid_inds = np.unravel_index(cross_inds, (sum(mask), fine_grid_size))

            # Time indexes in full ray -- index of confirmed crossings in the masked-off array
            time_inds = np.where(mask)[0][mask_inds]

            # Interpolate everything else
            df_fine = R.dataframe_at(interp_grid=interp_grid[fine_grid_inds,:], mask=time_inds)
            df_fine['cross_coords'] = cross_coords  # Cartesian coordinates of start and endpoints (for plotting)
            df_fine['EA_index'] = ea_inds           # Index of which EA array it intercepts


            print np.shape(df_fine)

            out_data = pd.concat([out_data, df_fine])

    tstop = time.time()
    print "Elapsed time (Interpolation and crossing detection): %g seconds"%(tstop - tstart) 
    return out_data



def calc_resonant_pitchangle_change(crossing_df, L_target):

    # Start timer
    tstart = time.time()
    # # Generate energy and velocity arrays
    # E_tot_arr = pow(10,sc.E_EXP_BOT + sc.DE_EXP*np.arange(0,sc.NUM_E))
    # v_tot_arr = sc.C*np.sqrt(1 - pow(sc.E_EL/(sc.E_EL + E_tot_arr),2))


    L = L_target

    epsm = (1.0/L)*(sc.R_E + sc.H_IONO)/sc.R_E
    alpha_eq = np.arcsin(np.sqrt( pow(epsm, 3)/np.sqrt(1 + 3*(1 - epsm))  ))


    # Gen EA array (maybe pass this as an arg)
    EA_array = gen_EA_array(L_target)

    # Initialize alpha array

    DA_array_N = np.zeros((sc.NUM_E, sc.NUM_STEPS))
    DA_array_S = np.zeros((sc.NUM_E, sc.NUM_STEPS))


    tstop = time.time()

    # Loop through EA segments:
    #for EA_ind, EA_row in EA_array.iterrows():
    for EA_ind in np.unique(crossing_df['EA_index']):   
    # For each EA segment, CALCULATE:
    #     - wh
    #     - dwh/ds
    #     - flight-time constant
    #     - alpha_lc


        lat = EA_array['lam'][EA_ind]
        print 'EA segment at latitude = ',lat

        slat = np.sin(lat*sc.D2R)
        clat = np.cos(lat*sc.D2R)
        slat_term = np.sqrt(1 + 3*slat*slat)

        # Lower hybrid frequency
        wh = (sc.Q_EL*sc.B0/sc.M_EL)*(1.0/pow(L,3))*slat_term/pow(clat,6)
        #wh = (2*np.pi*880000) / (pow(L,3)*slat_term/pow(clat,6))       # Jacob used a hard coded value for qB/m
        dwh_ds = (3.0*wh/(L*sc.R_E)) *(slat/slat_term) * (1.0/(slat_term*slat_term) + 2.0/(clat*clat))

        # flight time constants (divide by velocity to get time to ionosphere, in seconds)
        ftc_N, ftc_S = get_flight_time_constant(L, lat, alpha_eq)

        # Local loss cone angle
        alpha_lc = np.arcsin(np.sqrt(slat_term/pow(clat,6) )*np.sin(alpha_eq))
        salph = np.sin(alpha_lc)
        calph = np.cos(alpha_lc)

        ds = L*sc.R_E*slat_term*clat*sc.EAIncr*np.pi/180.0
        dv_para_ds = -0.5*(salph*salph/(calph*wh))*dwh_ds


        # # Sometimes Python is pretty and readable, but now is not one of those times
        # for cell_ind, cells in crossing_df[crossing_df['EA_index'] == EA_ind].iterrows():
        #     print np.shape(cell)

        # Mask out just the entries which cross the current EA segment:
        cells = crossing_df[crossing_df['EA_index'] == EA_ind]
        #print cells
        for ir, c in cells.iterrows():
            #print c


            # where you left off: time and frequency (line 1800)
            t = c.tg + sc.T_STEP/2.0
            f = c.frequency + sc.F_STEP/2.0
            # Jacob divides by num_rays too -- since he stores total power in each cell. Average?
            # Do we need to do that here? I'm not doing that right now.
            pwr = c.power/(sc.T_STEP*EA_array['EA_length'][EA_ind])      
            psi = c.psi*sc.D2R

            # print "EA_length: ",EA_array['EA_length'][EA_ind]
            # print "Pwr: ",pwr       
            # print "DT: ",sc.DT
            #print pwr

            # Maybe threshold power here?

            # Misc. parameters and trig
            mu = c.mu
            stixP = c.stixP
            stixR = c.stixR
            stixL = c.stixL

            spsi = np.sin(psi)
            cpsi = np.cos(psi)
            spsi_sq = spsi*spsi
            cpsi_sq = cpsi*cpsi
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

            #print num
            Byw_sq = ( (2.0*sc.MU0/sc.C) * (pwr*stixX*stixX*rho2*rho2*mu*abs(cpsi)) /
                     np.sqrt( pow((np.tan(psi)-rho1*rho2*stixX),2) + pow((1+rho2*rho2*stixX),2)) )


            print ("t: %g, f: %g, pwr: %g, Byw_sq: %g"%(t,f,pwr,Byw_sq))
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
            for mres in np.linspace(-5,5,11):

                # Parallel resonance velocity V_z^res  (pg 5 of Bortnik et al)
                t1 = w*w*kz*kz
                t2 = pow(mres*wh, 2)-(w*w)
                t3 = kz*kz + pow((mres*wh),2)/(pow(sc.C*np.cos(alpha_lc),2))
                

                #print np.diff(t3[0,:])
                if mres==0:
                    direction = -1.0*np.sign(kz)
                else:
                    direction = np.sign(kz)*np.sign(mres)

                v_para_res = ( direction*np.sqrt(t1 + t2*t3) - w*kz) / t3

                # # Getting some obnoxious numerical stuff in the case of m=0:
                # if mres==0:
                #     v_para_res = (w/kz)
                # else: 
                #     v_tot_res = v_para_res / np.cos(alpha_lc)

                v_tot_res = v_para_res / np.cos(alpha_lc)
                E_res = sc.E_EL*(1.0/np.sqrt( 1-(v_tot_res*v_tot_res/(sc.C*sc.C)) ) - 1 )

                # Start and end indices in the energy array (+-20% energy band)
                e_starti = np.floor((np.log10(E_res)-sc.E_EXP_BOT-sc.E_BANDWIDTH)/sc.DE_EXP)
                e_endi   = np.ceil(( np.log10(E_res)-sc.E_EXP_BOT+sc.E_BANDWIDTH)/sc.DE_EXP)


                # Threshold to range of actual indexes (some energies may be outside our target range)
                e_starti = np.clip(e_starti,0,sc.NUM_E)
                e_endi   = np.clip(e_endi,  0,sc.NUM_E)

               

                # energy bins to calculate at:
                evec_inds = np.arange(e_starti,e_endi,dtype=int)

                v_tot = direction*sc.v_tot_arr[evec_inds]
                v_para = v_tot*calph
                v_perp = abs(v_tot*salph)

                # Relativistic factor
                gamma = 1.0/np.sqrt(1.0 - pow(v_tot/sc.C, 2))

                alpha2 = sc.Q_EL*Ezw /(sc.M_EL*gamma*w1*v_perp)
                
                beta = kx*v_perp/wh

                wtau_sq = (pow((-1),(mres-1)) * w1/gamma * 
                            ( jn( (mres-1), beta ) - 
                              alpha1*jn( (mres+1) , beta ) +
                              gamma*alpha2*jn( mres , beta ) ))

                T1 = -wtau_sq*(1 + calph*calph/(mres*Y - 1))

                # Analytical evaluation! (Line 1938)
                if (abs(lat) < 1e-3):        

                    eta_dot = mres*wh/gamma - w - kz*v_para

                    eta_mask = eta < 10
                    dalpha_eq[eta_mask]  = abs(T1[eta_mask] /v_para[eta_mask])*ds/np.sqrt(2)
                    dalpha_eq[~eta_mask] = abs(T1[~eta_mask]/eta_dot[~eta_mask])*np.sqrt(1-np.cos(ds*eta_dot[~eta_mask]/v_para[~eta_mask]))

                else:
                    v_para_star = v_para - dv_para_ds*ds/2.0
                    v_para_star_sq = v_para_star*v_para_star

                    AA =( (mres/(2.0*v_para_star*gamma))*dwh_ds* 
                          (1 + ds/(2.0*v_para_star)*dv_para_ds) - 
                           mres/(2.0*v_para_star_sq*gamma)*wh*dv_para_ds + 
                           w/(2.0*v_para_star_sq)*dv_para_ds )

                    BB =( mres/(gamma*v_para_star)*wh - 
                          mres/(gamma*v_para_star)*dwh_ds*(ds/2.0) -
                          w/v_para_star - kz )

                    Farg = (BB + 2.0*AA*ds) / np.sqrt(2.0*np.pi*abs(AA))
                    Farg0 = BB / np.sqrt(2*np.pi*abs(AA))

                    Fs,  Fc  = fresnel(Farg)
                    Fs0, Fc0 = fresnel(Farg0)

                    dFs_sq = pow(Fs - Fs0, 2)
                    dFc_sq = pow(Fc - Fc0, 2)

                    dalpha = np.sqrt( (np.pi/4.0)/abs(AA))*abs(T1/v_para)*np.sqrt(dFs_sq + dFc_sq)
                    alpha_eq_p = np.arcsin( np.sin(alpha_lc+dalpha)*pow(clat,3) / np.sqrt(slat_term) )
                    dalpha_eq  = alpha_eq_p - alpha_eq
                    #print "dalpha_eq:", dalpha_eq
                    #print "output vec shape:", np.shape(dalpha_eq)
                    #print "end - start:", endi - starti

                    #print "t offset:",t.iloc[index[0]]
                    # Add net change in alpha to total array:
                    if direction > 0: 
                        #print "Norf!"

                        tt = np.round( (t + ftc_N/v_para)/sc.T_STEP ).astype(int)
                        #print tt
                        tt = np.clip(tt,0,sc.NUM_STEPS - 1)
                        #tt.clip(0,sc.NUM_STEPS)
                        DA_array_N[evec_inds,tt] += dalpha_eq*dalpha_eq
                    else:
                        #print "Souf!" 
                        tt = np.round( (t + ftc_S/v_para)/sc.T_STEP ).astype(int)
                        tt = np.clip(tt,0,sc.NUM_STEPS - 1)
                        #print tt
                        DA_array_S[evec_inds,tt] += dalpha_eq*dalpha_eq






        # --------------------- Vector-ish version that isn't working ---------------
        # if cells.shape[0] > 0:
        #     print cells.shape
        #     #print cells.columns

        #     # where you left off: time and frequency (line 1800)
        #     t = cells.tg + sc.DT/2.0
        #     f = cells.frequency + sc.DF/2.0
        #     pwr = cells.power/sc.DT       # Jacob divides by num_rays too, but it looks like it's always 1
        #     psi = cells.psi*sc.D2R       

        #     #print pwr

        #     # Maybe threshold power here?

        #     # Misc. parameters and trig
        #     mu = cells.mu
        #     stixP = cells.stixP
        #     stixR = cells.stixR
        #     stixL = cells.stixL

        #     spsi = np.sin(psi)
        #     cpsi = np.cos(psi)
        #     spsi_sq = spsi*spsi
        #     cpsi_sq = cpsi*cpsi
        #     n_x  = mu*abs(spsi)
        #     n_z  = mu*cpsi
        #     mu_sq = mu*mu
        #     w = 2.0 * np.pi * f
        #     k = w*mu/sc.C
        #     kx = w*n_x/sc.C
        #     kz = w*n_z/sc.C
        #     Y = wh / w

        #     stixS = (stixR + stixL)/2.0
        #     stixD = (stixR - stixL)/2.0
        #     stixA = stixS + (stixP-stixS)*cpsi_sq
        #     stixB = stixP*stixS + stixR*stixL + (stixP*stixS - stixR*stixL)*cpsi_sq
        #     stixX = stixP/(stixP - mu_sq*spsi_sq)
          
        #     rho1=((mu_sq-stixS)*mu_sq*spsi*cpsi)/(stixD*(mu_sq*spsi_sq-stixP))
        #     rho2 = (mu_sq - stixS) / stixD

        #     #print num
        #     Byw_sq = ( (2.0*sc.MU0/sc.C) * (pwr*stixX*stixX*rho2*rho2*mu*abs(cpsi)) /
        #              np.sqrt( pow((np.tan(psi)-rho1*rho2*stixX),2) + pow((1+rho2*rho2*stixX),2)) )
        #     #print Byw_sq
        #     # RMS wave components
        #     Byw = np.sqrt(Byw_sq);
        #     Exw = abs(sc.C*Byw * (stixP - n_x*n_x)/(stixP*n_z))
        #     Eyw = abs(Exw * stixD/(stixS-mu_sq))
        #     Ezw = abs(Exw *n_x*n_z / (n_x*n_x - stixP))
        #     Bxw = abs((Exw *stixD*n_z/sc.C)/ (stixS - mu_sq))
        #     Bzw = abs((Exw *stixD *n_x) /(sc.C*(stixX - mu_sq)));

        #     # Oblique integration quantities
        #     R1 = (Exw + Eyw)/(Bxw+Byw)
        #     R2 = (Exw - Eyw)/(Bxw-Byw)
        #     w1 = (sc.Q_EL/(2*sc.M_EL))*(Bxw+Byw)
        #     w2 = (sc.Q_EL/(2*sc.M_EL))*(Bxw-Byw)
        #     alpha1 = w2/w1
          

        #     # Loop at line 1897 -- MRES loop
        #     # Since we're trying to keep it vectorized:
        #     # target dimensions are num(cell rows) rows x num(mres) columns.
        #     # (Sometimes Python is beautiful, but now is not one of those times)
        #     #  --> Note to future editors: the (vec)[:,np.newaxis] command is similar to repmat or tile:
        #     #  It will tile the vector along the new axis, with length dictated by whatever operation 
        #     #  follows it. In Python lingo this is "broadcasting"

        #     # Resonant modes to sum over: Integers, positive and negative
        #     mres = np.linspace(-5,5,11)
            
        #     # Parallel resonance velocity V_z^res  (pg 5 of Bortnik et al)
        #     t1 = w*w*kz*kz
        #     t2 = (mres*mres*wh*wh)[np.newaxis,:]-(w*w)[:,np.newaxis]
        #     t3 = (kz*kz)[:,np.newaxis] + ((mres*wh*mres*wh)[np.newaxis,:]/(pow(sc.C*np.cos(alpha_lc),2)))
            

        #     #print np.diff(t3[0,:])
        #     direction = np.outer(np.sign(kz),np.sign(mres))
        #     direction[:,mres==0] = -1*np.sign(kz)[:,np.newaxis] # Sign 

        #     # v_para_res = (ta - tb)
        #     v_para_res = ( direction*np.sqrt(abs(t1[:,np.newaxis] + t2*t3)) - (w*kz)[:,np.newaxis] ) / t3

        #     # Getting some obnoxious numerical stuff in the case of m=0:
        #     v_para_res[:,mres==0] = (w/kz)[:,np.newaxis]
        #     v_tot_res = v_para_res / np.cos(alpha_lc)

        #     E_res = sc.E_EL*(1.0/np.sqrt( 1-(v_tot_res*v_tot_res/(sc.C*sc.C)) ) - 1 )

        #     # Start and end indices in the energy array (+-20% energy band)
        #     e_starti = np.floor((np.log10(E_res)-sc.E_EXP_BOT-sc.E_BANDWIDTH)/sc.DE_EXP)
        #     e_endi   = np.ceil(( np.log10(E_res)-sc.E_EXP_BOT+sc.E_BANDWIDTH)/sc.DE_EXP)


        #     # Threshold to range of actual indexes (some energies may be outside our target range)
        #     np.clip(e_starti,0,sc.NUM_E, out=e_starti)
        #     np.clip(e_endi,  0,sc.NUM_E, out=e_endi)

        #     #print "Start index:",e_starti
        #     #print "Stop index: ",e_endi


        #     #inds = zip(e_starti.ravel(), e_endi.ravel())
        #     #print np.shape(inds)
        #     #v_tot = direction*v_tot_arr[e_]


        #     # emask = np.zeros((len(cells), len(mres), sc.NUM_E),dtype=bool)



        #     # for index, starti in np.ndenumerate(e_starti):
        #     #     emask[index[0], index[1],e_starti[index]:e_endi[index]] = True

        #     # print "total emask: ",np.sum(emask)

        #     # emask[:,:,e_starti:e_endi] = True
        #     #print "emask is: ", np.shape(emask)

        #     # Where you left off: Adding in the next loop (V_TOT). Need to iterate from estarti to eendi
        #     # for each cell in the 2d arrays... hmm.... how to vectorize nicely. Hmm.
        #     #print "e_starti is", np.shape(e_starti)
        #     for index, starti in np.ndenumerate(e_starti):
        #         endi = e_endi[index]

        #         #print starti, endi
        #         # Energies to do
        #         #print Ezw.index.values
                    
        #         #print Ezw[index(0)]

        #         if endi >= starti:
        #             # Dimensions of matrices: (num cells) rows x (num energies) columns
        #             #print "Index is", index

        #             v_tot = direction[index]*sc.v_tot_arr[starti:endi]
        #             v_para = v_tot*calph
        #             v_perp = abs(v_tot*salph)

        #             # Relativistic factor
        #             gamma = 1.0/np.sqrt(1.0 - pow(v_tot/sc.C, 2))

        #             alpha2 = (sc.Q_EL/sc.M_EL)*(Ezw.iloc[index[0]]/w1.iloc[index[0]])/(gamma*v_perp)
                    
        #             #print "gamma:", gamma
        #             #print "alpha2:",alpha2
        #             # beta = np.outer(kx/wh,v_perp)
        #             beta = (kx.iloc[index[0]]/wh)*v_perp
        #             #print "beta:", beta

        #             mr = mres[index[1]]
        #             # wtau_sq = (pow((-1),(mr-1))*np.outer(w1, 1.0/gamma) *
        #             #                     (jn(mr -1, beta) -
        #             #                      alpha1[:,np.newaxis]*jn(mr + 1, beta) +
        #             #                      gamma[np.newaxis,:]*alpha2*jn(mr, beta)))
        #             wtau_sq = (pow((-1),(mr-1))*(w1.iloc[index[0]]/gamma)*
        #                        jn(mr - 1, beta) -
        #                        alpha1.iloc[index[0]]*jn(mr + 1, beta) +
        #                        gamma*alpha2*jn(mr,beta))

        #             #T1 = -wtau_sq*(1 + calph*calph/(mr*Y[:,np.newaxis] - 1))
        #             T1 = -wtau_sq*(1 + calph*calph/(mr*Y.iloc[index[0]] - 1))

        #             # Analytical evaluation! (Line 1938)
        #             if (abs(lat) < 1e-3):
        #                 eta_dot = mr*wh/gamma - w.iloc[index[0]] - kz.iloc[index[0]]*v_para
        #                 dalpha_eq = abs(T1/eta_dot)*np.sqrt(1 - np.cos(ds*eta_dot/v_para))
        #             else:
        #                 v_para_star = v_para - dv_para_ds*ds/2.0
        #                 v_para_star_sq = v_para_star*v_para_star

        #                 AA = ((mr/(2.0*v_para_star*gamma))*dwh_ds * 
        #                       (1 + ds/(2.0*v_para_star)*dv_para_ds) - 
        #                       (mr/(2.0*v_para_star_sq*gamma))*wh*dv_para_ds + 
        #                       (w.iloc[index[0]]/(2.0*v_para_star_sq))*dv_para_ds)

        #                 BB = ((mr/(gamma*v_para_star))*wh - 
        #                       (mr/(gamma*v_para_star))*dwh_ds*(ds/2.0) -
        #                        w.iloc[index[0]]/v_para_star - kz.iloc[index[0]])

        #                 Farg = (BB + 2.0*AA*ds) / np.sqrt(2.0*np.pi*abs(AA))
        #                 Farg0 = BB / np.sqrt(2*np.pi*abs(AA))

        #                 Fs,  Fc  = fresnel(Farg)
        #                 Fs0, Fc0 = fresnel(Farg0)

        #                 dFs_sq = pow(Fs - Fs0, 2)
        #                 dFc_sq = pow(Fc - Fc0, 2)

        #                 dalpha = np.sqrt( (np.pi/4.0)*abs(AA))*abs(T1/v_para)*np.sqrt(dFs_sq + dFc_sq)
        #                 alpha_eq_p = np.arcsin( np.sin(alpha_lc+dalpha)*pow(clat,3) / np.sqrt(slat_term) )
        #                 dalpha_eq  = alpha_eq_p - alpha_eq
        #                 #print dalpha_eq
        #                 #print "output vec shape:", np.shape(dalpha_eq)
        #                 #print "end - start:", endi - starti

        #                 #print "t offset:",t.iloc[index[0]]
        #                 # Add net change in alpha to total array:
        #                 if direction[index] > 0:
        #                     #print "Norf!"

        #                     tt = np.round( (t.iloc[index[0]] + ftc_N/v_para)/sc.T_STEP ).astype(int)
        #                     #print tt
        #                     np.clip(tt,0,sc.NUM_STEPS - 1, out=tt)
        #                     #tt.clip(0,sc.NUM_STEPS)
        #                     DA_array_N[starti:endi,tt] += dalpha_eq*dalpha_eq
        #                 else:
        #                     #print "Souf!"
        #                     tt = np.round( (t.iloc[index[0]] + ftc_S/v_para)/sc.T_STEP ).astype(int)
        #                     np.clip(tt,0,sc.NUM_STEPS - 1, out=tt)
        #                     #print tt
        #                     DA_array_S[starti:endi,tt] += dalpha_eq*dalpha_eq
    tstop = time.time()

    print "Scattering took %g seconds"%(tstop - tstart)
    return DA_array_N, DA_array_S








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

    out_list = dict()


    ### --------- Jacob's Version ------
    # # Vector l, no pre-masking
    # r1ray = l_prev_vec * pow( np.cos(lat_prev_vec*sc.D2R) , 2 )
    # x1ray = r1ray * np.cos(lat_prev_vec*sc.D2R)
    # y1ray = r1ray * np.sin(lat_prev_vec*sc.D2R)

    # # Ugh. I hate that this works -- crossing detection alg uses r1 for x2 and y2. But Jacob's shit seems good.
    # r2ray = l * pow( np.cos(lat*sc.D2R) , 2 )
    # x2ray = r1ray * np.cos(lat*sc.D2R)
    # y2ray = r1ray * np.sin(lat*sc.D2R)

    # Aray = y1ray - y2ray
    # Bray = x2ray - x1ray
    # Cray = x1ray*y2ray - y1ray*x2ray

    # val1 = np.outer(Aray,EA_array['x1']) + np.outer(Bray,EA_array['y1']) + Cray[:,np.newaxis]
    # val2 = np.outer(Aray,EA_array['x2']) + np.outer(Bray,EA_array['y1']) + Cray[:,np.newaxis]

    # val3 = np.outer(x1ray,EA_array['EAa']) + np.outer(y1ray,EA_array['EAb']) + EA_array['EAc'][np.newaxis,:]
    # val4 = np.outer(x2ray,EA_array['EAa']) + np.outer(y2ray,EA_array['EAb']) + EA_array['EAc'][np.newaxis,:]

    # mask = np.array(val1*val2 <= 0) & (val3*val4 <= 0)

    # ---------- My Version -------------
    # (Does not check for parallel rays, doesn't watch for divide-by-zero cases... might be problematic?)

    r1ray = l_prev_vec * pow( np.cos(lat_prev_vec*sc.D2R) , 2 )
    x1ray = r1ray * np.cos(lat_prev_vec*sc.D2R)
    y1ray = r1ray * np.sin(lat_prev_vec*sc.D2R)

    r2ray = l * pow( np.cos(lat*sc.D2R) , 2 )
    x2ray = r2ray * np.cos(lat*sc.D2R)
    y2ray = r2ray * np.sin(lat*sc.D2R)

    A1 = y1ray - y2ray
    B1 = x2ray - x1ray
    C1 = x1ray*y2ray - x2ray*y1ray

    A2  = EA_array['EAa']
    B2  = EA_array['EAb']
    C2  = EA_array['EAc']

    x2EA= EA_array['x2']
    x1EA= EA_array['x1']

    # Xa is the X-value of intersection point
    Xa = (np.outer(B1,C2) - np.outer(C1,B2))/(np.outer(A1,B2) - np.outer(B1,A2))
    #print Xa

    # Check if within bounds of the segments:
    Imin = np.zeros_like(Xa)
    Imax = np.zeros_like(Xa)

    # minimum value: max(min(x2ray, x1ray), min(x2ea, x1ea))
    Imin = np.maximum(np.minimum(np.tile(x2ray, (len(x2EA),1)).transpose(), np.tile(x1ray, (len(x2EA),1)).transpose()),
                      np.minimum(np.tile(x2EA,  (len(x2ray),1)),            np.tile(x1EA, (len(x2ray),1))) )

    # maximum value: min(max(x2ray, x1ray), max(x2ea, x1ea))
    Imax = np.minimum(np.maximum(np.tile(x2ray, (len(x2EA),1)).transpose(), np.tile(x1ray, (len(x2EA),1)).transpose()),
                      np.maximum(np.tile(x2EA,  (len(x2ray),1)),            np.tile(x1EA, (len(x2ray),1))) )
    
    # for ind, ival in np.ndenumerate(Xa):
    #     ix = ind[0]
    #     iy = ind[1]
    #     #print ix, iy
    #     Imin[ix,iy] = max(min(x2ray[ix],x1ray[ix]),min(x2EA[iy],x1EA[iy]))
    #     Imax[ix,iy] = min(max(x2ray[ix],x1ray[ix]),max(x2EA[iy],x1EA[iy]))

    mask = (Xa >= Imin) & (Xa <= Imax)
    # print mask
    # print np.sum(mask)
    # print np.shape(mask)

    # Rows and columns where mask is true:
    mr, mc = np.where(mask)

    # mr: index in l and lat
    # mc: index in EA array
    #return mr, mc

    # row: index in time series
    # col: index in EA

    # # Same format as before (with Cartesian coordinates for the ray segments)
    out_inds = zip(mr,mc)


    x2ray = r2ray * np.cos(lat*sc.D2R)
    y2ray = r2ray * np.sin(lat*sc.D2R)

    out_cartesian = []
    for o in out_inds:
        out_cartesian.append([(x1ray[o[0]], y1ray[o[0]]),(x2ray[o[0]],y2ray[o[0]])])

    return mr, mc, out_cartesian


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
        interpolater = interpolate.interp1d(df[t_label],df[c],bounds_error=False)
        tmp[c] = interpolater(t_new)
    
    return tmp



def lightning_power(I0, center_lat, dlat, dfreq, frequency, ray_lat, dlong = 0.7):
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
    dist_lat = (sc.R_E + sc.H_IONO/2.0)*abs(center_lat - ray_lat)*sc.D2R    
    dist_long= (sc.R_E + sc.H_IONO/2.0)*abs(dlong)*sc.D2R
    #print "dist_lat: ",dist_lat,"dist_long: ",dist_long

    dist_iono = np.hypot(dist_lat, dist_long)
    dist_tot =  np.hypot(dist_iono, sc.H_IONO)
    #print "Dist_iono: ", dist_iono,"dist_tot: ",dist_tot
    xi = np.arctan2(dist_iono, sc.H_IONO)
    #print "xi: ",xi
    # Power at sub-ionosphere point
    w = 2*np.pi*frequency
    w_sq = pow(w, 2)

    S = ( 1.0/(sc.Z0) * pow( sc.H_E * I0 * (2e-7)*(np.sin(xi)/dist_tot)* w * (sc.A - sc.B),2)
        /( (w_sq + pow(sc.A, 2))*(w_sq + pow(sc.B,2)))  )

    # Get ionosphere attenuation:
    atten_factor = pow(10, -ionoAbsorp(ray_lat, frequency)/10.0)

    S_vert = S * np.cos(xi)*atten_factor   # Vertical component

    # print ("I0: %2.3f, center_lat: %2.3f, dlat: %2.3f, dfreq: %2.3f, f: %ld, lat: %2.3f, dlong: %g, S_vert: %e"
    #         %(I0, center_lat, dlat, dfreq, frequency, ray_lat, dlong, S_vert))

    print "S_vert:", S_vert
    # I don't know where Jacob got this attenuation factor from. Geometric spreading? Hmm.
    return S_vert * dlat * sc.D2R * (sc.R_E + sc.H_IONO) * dfreq * 0.87788331






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
    rootTerm = np.sqrt(1 + 3.0*slam2)

    # dl_lam: still not sure where he gets this expression from
#    dl_lam = (clam2*clam / rootTerm)*sc.DL0
    dl_lam = (clam2*clam / rootTerm)*sc.L_MARGIN

    x_unit_vect = (3.0*clam2 - 2) / rootTerm
    y_unit_vect = (3.0*slam*clam) / rootTerm


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
    d['EA_length'] = np.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))*sc.R_E
    d['x1'] = x1
    d['x2'] = x2
    d['y1'] = y1
    d['y2'] = y2

    # print "EA_length:", d['EA_length']

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

        # print self.BL['tg']

    def interpolation_weights(self, interp_grid):

        assert np.shape(interp_grid)[1] == 2, "Interpolation grid is not two columns!"

        BL_fact = 1 - interp_grid[:,0]-interp_grid[:,1] + interp_grid[:,0]*interp_grid[:,1]
        TL_fact = interp_grid[:,0] - interp_grid[:,0]*interp_grid[:,1]
        BH_fact = interp_grid[:,1] - interp_grid[:,0]*interp_grid[:,1]
        TH_fact = interp_grid[:,0]*interp_grid[:,1]

        return BL_fact, TL_fact, BH_fact, TH_fact
    def all_vals_at(self, interp_grid, val_name, mask):
        '''Returns interpolated values at each time step, and each factor.
           interp_grid: two-column array -- (l, f) interpolating values, 0 to 1.
           val_name: name of the column to interpolate
           mask: binary mask of time bins to look at '''
        # Interpolating weights:

        BL_fact, TL_fact, BH_fact, TH_fact = self.interpolation_weights(interp_grid)

        return (    np.outer(self.BL[val_name][mask], BL_fact) + 
                    np.outer(self.BH[val_name][mask], BH_fact) + 
                    np.outer(self.TL[val_name][mask], TL_fact) + 
                    np.outer(self.TH[val_name][mask], TH_fact)).ravel()

    def vals_at(self, interp_grid, val_name, mask):
        ''' Same thing, but without the outer product. mask and factors must be same length.
           BL, TL, BH, TH: Interpolating weights (vectors).
           val_name: name of the column to interpolate
           mask: binary mask of time bins to look at '''

        BL_fact, TL_fact, BH_fact, TH_fact = self.interpolation_weights(interp_grid)

        return (self.BL[val_name][mask]*BL_fact + 
                self.BH[val_name][mask]*BH_fact + 
                self.TL[val_name][mask]*TL_fact + 
                self.TH[val_name][mask]*TH_fact)

    def dataframe_at(self, interp_grid, mask):

        cols = self.BL.columns.values
        df = pd.DataFrame()

        for c in cols:
            df[c] = self.vals_at(interp_grid, c, mask)

        BL_fact, TL_fact, BH_fact, TH_fact = self.interpolation_weights(interp_grid)

        df['frequency'] = (self.BL.frequency*BL_fact + 
                           self.BH.frequency*BH_fact + 
                           self.TL.frequency*TL_fact + 
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
