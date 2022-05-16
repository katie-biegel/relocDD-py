import numpy as np
import datetime
import sys
import matplotlib.pyplot as plt
import os

# Import other hypoDD-py functions
from hypoDD.hypoDD_functions import cluster,dataprep,trialsrc,dtres,weighting,resstat
from hypoDD.hypoDD_functions import skip,sigcoherency,statres,eventstats
from hypoDD.inversion import svd,lsqr
from hypoDD.hypoDD_files import writeloc,writereloc,terminaloutputs
from utility.universal.geodetics import setorg,sdc2
from utility.universal.raytrace import partials

# Main HYPODD Function
def hypoDD(log,hinputs,outfol='outputs',reloctype=1,fileout=0,hypoDDdata=[],iboot=None):
    """
    Taken from hypoDD v1.3 11/2010
    HypoDD v1.3 Author: Felix Waldhauser
    Translated and edited: K. Biegel, katherine.biegel@ucalgary.ca
    :::
    Purpose:
    Program to determine high-resolution hypocenter locations using the
    double-difference algorithm. Residuals between observed and theoretical 
    travel time differences (or double-differences) are minimized for pairs
    of earthquakes at each station while linking together all observed
    event/station pairs. A least squares solution (SVD or LSQR) is found
    by iteratively adjusting the vector difference between hypocentral pairs.
    :::
    References:
    For a detailed description of the algorithm see:
    Waldhauser, F. and W.L. Ellsworth, A double-difference earthquake
        location algorithm: Method and application to the northern Hayward
        fault, Bull. Seismol. Soc. Am., 90, 1353-1368, 2000.

    For a user guide to hypoDD see USGS open-file report: 
        Waldhauser, F., HypoDD: A computer program to compute double-difference
        earthquake locations,  U.S. Geol. Surv. open-file report , 01-113,
        Menlo Park, California, 2001.

    Further publication citations will follow for the python
    specific version of this work.
    :::
    """

    # Start statement to log
    datet = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    log.write('\n\nstarting hypoDD (python version) %s \n' % datet)
    print('\n\nstarting hypoDD (python version) %s' % datet)
    
    # Make sub-output folder for textfiles
    txtfol = os.path.join(outfol,'txtoutputs')
    os.makedirs(txtfol,exist_ok='True')

    """
    Set up data and arrays
    """
    # Unpack input parameters
    [fn_cc,fn_ct,fn_sta,fn_eve,fn_loc,fn_reloc,fn_res,fn_stares,fn_srcpar,
     idata,iphase,minobs_cc,minobs_ct,amaxres_cross,amaxres_net,amaxdcc,amaxdct,maxdist,
     awt_ccp,awt_ccs,awt_ctp,awt_cts,adamp,istart,maxiter,isolv,niter,aiter,
     mod_nl,mod_ratio,mod_v,mod_top,iclust,ncusp,icusp] = hinputs

    # Get Data
    if fileout==0:
        """
        Read data
        """
        [ev_date,ev_time,ev_cusp,ev_lat,ev_lon,ev_dep,ev_mag,ev_herr,ev_zerr,ev_res,
         sta_lab,sta_lat,sta_lon,dt_sta1,dt_sta2,dt_dt,dt_qual,dt_c1,dt_c2,dt_idx,
         dt_ista1,dt_ista2,dt_ic1,dt_ic2,dt_offse,dt_offss,nev,nsta,ndt,ncc,nct,nccp,
         nccs,nctp,ncts]= dataprep(log,reloctype,fileout,fn_cc,fn_ct,fn_sta,fn_eve,idata,
                                   iphase,ncusp,icusp,maxdist,amaxdct[0],amaxdcc[0])
    elif fileout==1:
        """
        Unpack data if passed into function
        """
        [ev_date,ev_time,ev_cusp,ev_lat,ev_lon,ev_dep,ev_mag,ev_herr,ev_zerr,ev_res,
         sta_lab,sta_lat,sta_lon,dt_sta1,dt_sta2,dt_dt,dt_qual,dt_c1,dt_c2,dt_idx,
         dt_ista1,dt_ista2,dt_ic1,dt_ic2,dt_offse,dt_offss,nev,nsta,ndt,nccp,nccs,
         nctp,ncts] = hypoDDdata
        ncc = nccp+nccs
        nct = nctp+ncts
    else:
        raise Exception('(hypoDD) Fileout value error. HypoDD Function.')

    # if iboot:
    #     fn_loc = iboot + fn_loc
    #     fn_reloc = iboot + fn_reloc
    #     fn_res = iboot + fn_res 
    #     fn_stares = iboot + fn_stares
    #     fn_srcpar = iboot + fn_srcpar

    """
    Declare Open Variables
    """
    tmpr1 = float(0.)
    tmpr2 = float(0.)
    minwght = float(0.00001)
    rms_ccold = float(0.)
    rms_ctold = float(0.)
    rms_cc0old = float(0.)
    rms_ct0old = float(0.)
    ineg = int(0)
    rms_cc = float(0.)
    rms_ct = float(0.)
    rms_cc0 = float(0.)
    rms_ct0 = float(0.)
    xav = float(0.)
    yav = float(0.)
    zav = float(0.)
    tav = float(0.)
    sav = float(0.)
    alat = float(0.)
    alon = float(0.)
    adep = float(0.)
    
    """
    Reset .reloc and .loc file
    """
    loc = open(fn_loc,'w') 
    loc.close()
    reloc = open(fn_reloc,'w')
    reloc.close()

    """"
    Initial Clustering
    ------
    IF:
        Per data type if minobs_* is set to 0 in the input file this means no clustering.
        All events are then in one cluster:
            nclust therefore = 1
            and all event ids are including in clust array in first cluster
    ELSE:
        Function Cluster1 is called from script hypoDD_functions.py
        This function sorts events into clusters based on the minimum number of 
        shared data required. If an event-pair has > minobs_* both events are in
        the same cluster.
        * Note this is dependent on the minobs # in ph2dt so should be larger than
        that number if clustering is required.
    """
    if (idata==1 and minobs_cc==0) or (idata==2 and minobs_ct==0) or (idata==3 and minobs_ct+minobs_cc==0):
        # No clustering (one cluster includes all events)
        nclust = 1
        clust = np.zeros((1,nev+1))
        clust[0,0] = nev
        for i in range(0,nev):
            clust[0,i+1] = ev_cusp[i]
        log.write('(hypoDD) No clustering performed. \n')
        if fileout==0:
            print('(hypoDD) No clustering performed.')
    elif (minobs_ct+minobs_cc>0): # Clustering
        if reloctype==2:
            nclust = 1
            clust = np.zeros((1,nev+1))
            clust[0,0] = nev
            for i in range(0,nev):
                clust[0,i+1] = ev_cusp[i]

            # nclust = 100
            # clust = np.zeros((nclust,nev+1),dtype='int')
            # for i in range(100):
            #     try:
            #         evtoclust = open('txtoutputs/clust_%i.dat' % i,'r')
            #         clustlines = evtoclust.readlines()
            #         evtoclust.close()
            #     except:
            #         break

            #     count = 1
            #     for line in clustlines:
            #         line = line.strip()
            #         line = line.split(' ')
            #         line = list(filter(None,line))

            #         clust[i,count:count+len(line)] = line[0:len(line)]
            #         count += len(line)
            #     clust[i,0] = count-1

        else:
            [clust,noclust,nclust] = cluster(log,txtfol,nev,ndt,idata,minobs_cc,minobs_ct,
                                             dt_c1,dt_c2,ev_cusp)
    """
    Open Output Files if Writing to File
    """
    if fileout==0:
        if len(fn_stares)>1:
            stares = open(fn_stares,'w')
            stares.close()

    """
    Initalise loop over iterations
    """
    jiter=int(0)    # Counter for iter with no updating (air quakes)
                    # Total # bad iterations (not recorded) later subtracted from total
    
    # Big loop over clusters starts here:
    if iclust!=0:   # If specific clusters specified in input file
                    # Only allowed to specify one cluster
        if iclust<0 or iclust>nclust:   # Check if cluster ID is valid otherwise exit
            raise Exception('(hypoDD) Error: invalid cluster number %5i. Must be between 1 and nclust (%5i)' % (iclust,nclust))
        ibeg=iclust-1
        iend=iclust
    else:   # Otherwise loop over all clusters
        ibeg=int(0)
        iend=nclust+1     # nclust is the total number of clusters
        if nclust!=0:
            for i in range(nclust):
                if clust[i,0]<2:
                    print('Cluster %i has less than 2 events.  Only relocating first %i clusters.' % (i+1,i))
                    break
                else:
                    iend = i+1
    if iboot:
        ibeg = 0
        iend = 1

    # Initialise airquake array (cuspid's for airquakes later stored here)
    amcusp = np.zeros(1000,dtype='int')

    """
    Initialize variables
    """
    exav=int(0)
    eyav=int(0)
    ezav=int(0)
    etav=int(0)
    esav=int(0)
    dxav=int(0)
    dyav=int(0)
    dzav=int(0)
    dtav=int(0)
    dsav=int(0)

    """
    Big cluster loop
    ----
    All clusters are relocated separately.
    Relocation iterations called within clusters.
    """
    for icl in range(ibeg,iend):
        time_cluststart = datetime.datetime.now()
        datet = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        log.write('\n\n(hypoDD) RELOCATION OF CLUSTER: %2i     %s \n\n\n' % (icl+1,datet))
        print('\n(hypoDD) RELOCATION OF CLUSTER: %2i     %s' % (icl+1,datet))

        """
        Get data for each cluster if clustering invoked
        -----
        Reads in only correct data for each cluster
        """
        if (nclust!=1) and (minobs_cc>0 or minobs_ct>0):
            ncusp = int(clust[icl,0])
            icusp = clust[icl,1:ncusp+1] #np.zeros(ncusp,dtype='int')
            
            [ev_date,ev_time,ev_cusp,ev_lat,ev_lon,ev_dep,ev_mag,ev_herr,ev_zerr,ev_res,
             sta_lab,sta_lat,sta_lon,dt_sta1,dt_sta2,dt_dt,dt_qual,dt_c1,dt_c2,dt_idx,
             dt_ista1,dt_ista2,dt_ic1,dt_ic2,dt_offse,dt_offss,nev,nsta,ndt,ncc,nct,nccp,
             nccs,nctp,ncts]= dataprep(log,reloctype,fileout,fn_cc,fn_ct,fn_sta,fn_eve,
                                       idata,iphase,ncusp,icusp,maxdist,amaxdct[0],amaxdcc[0])
        # Recount all data (to update for cluster subset)
        nccold = ncc
        nctold = nct
        ncc = nccp+nccs
        nct = nctp+ncts
        nevold = nev

        # Initialize arrays for storing station residuals statistics
        sta_np = np.zeros(nsta,dtype='int')
        sta_ns = np.zeros(nsta,dtype='int')
        sta_nnp = np.zeros(nsta,dtype='int')
        sta_nns = np.zeros(nsta,dtype='int')
        sta_rmsc = np.zeros(nsta,dtype='float')
        sta_rmsn = np.zeros(nsta,dtype='float')

        # if iboot:
        #     if reloctype!=2:
        #         log.write('\nBootstrapping is running.  Randomly re-sampling data with replacement.\n')
        #         print('Bootstrapping is running.  Randomly re-sampling data with replacement.')
        #         dt_dt = np.random.choice(dt_dt,size=len(dt_dt),replace=True)
        #     elif reloctype==2:
        #         log.write('\nBootsrapping is running.  Adding noise to data to match event-pair noise level.\n')
        #         print('Bootstrapping is running.  Randomly re-sampling data with replacement.')
        #         dt_dt += np.random.normal(loc=0.,scale=0.005,size=len(dt_dt))
        #     log.write('Saving data.')
        #     np.savetxt('%sdt_postresample.txt',dt_dt)

        """
        Get cluster centroid
        -----
        From average of all events in cluster
        """
        sdc0_lat = 0.
        sdc0_lon = 0.
        sdc0_dep = 0.
        sdc0_lat = np.sum(ev_lat)
        sdc0_lon = np.sum(ev_lon)
        sdc0_dep = np.sum(ev_dep)
        sdc0_lat = sdc0_lat/nev
        sdc0_lon = sdc0_lon/nev
        sdc0_dep = sdc0_dep/nev
        log.write('Cluster centroid at: %10.6f  %11.6f  %9.6f \n' % (sdc0_lat,sdc0_lon,sdc0_dep))

        """"
        Set up cartesian coordinates from cluster centroid epicenter lat and lon
        ---
        This is important for coordinte conversions later on
        The cluster centroid lat and lon are set to the cartesian coordinate system origin point
        This means that relocation from cluster centroid (istart=1) all event locations start at (0,0,0)
        """
        setorg(sdc0_lat,sdc0_lon)

        """
        Convert all events to cartesian
        ---
        Get cartesian coordinates for epicenters (lat,lon,dep) to (x,y,z)
        """
        ev_x = np.zeros(nev)
        ev_y = np.zeros(nev)
        ev_z = np.zeros(nev)
        for i in range(0,nev):
            lat = ev_lat[i]
            lon = ev_lon[i]
            [x,y] = sdc2(lat,lon,-1)
            ev_x[i] = x*1000
            ev_y[i] = y*1000
            ev_z[i] = (ev_dep[i]-sdc0_dep)*1000
        log.write('# Events: %5i \n' % nev)

        """
        Save initial locations to file by cluster
        ---
        Write to output (mdat.loc)
        """
        if fileout==0:
            writeloc(fn_loc,nev,ev_cusp,ev_lat,ev_lon,ev_dep,ev_x,ev_y,ev_z,ev_herr,ev_zerr,
                     ev_date,ev_time,ev_mag)

        """
        Get initial trial sources
        """
        [nsrc,src_cusp,src_lat0,src_lon0,src_x0,src_y0,src_z0,src_t0,src_lat,
         src_lon,src_dep,src_x,src_y,src_z,src_t,src_xi,src_yi,src_zi,
         src_ti] = trialsrc(istart,sdc0_lat,sdc0_lon,sdc0_dep,nev,ev_cusp,ev_lat,ev_lon,ev_dep)
        src_s = np.zeros(nev,dtype='float')
        log.write('(hypoDD) # Initial trial sources: %6i \n' % nsrc)
        
        """
        Loop over iterations starts here:
        ------
        Define each iteration step at which re-weighting starts:
        this is dynam. since it depends on the number of neg 
        depths runs before (i.e. jiter must be subtracted out since
        these iterations are not recorded)
        """
        for i in range(0,niter):
            aiter[i] = aiter[i] - jiter     # Aiter from input file
                                            # Number of iterations per set
        maxiter = maxiter - jiter 
        kiter = 0   # counter for iter with data skipping
        jiter = 0   # counter for iter with no updating (air quakes)
        mbad = 0    # counter for air quakes
        normvar=np.zeros((int(maxiter),3))
        iteri = 0   # start on iteration 1

        # Output arrays
        calstart = np.empty(ndt)
        calend = np.empty(ndt)
        dtdtstart = np.empty(ndt)
        dtdtend = np.empty(ndt)
        loctru = np.empty((4,nev))
        locdel = np.empty((5,nev,int(maxiter)))
        locabs = np.empty((5,nev,int(maxiter)+1))
        time_preit = datetime.datetime.now()

        while iteri<maxiter:   # while iteri (current iteration) less than max number of iterations
            time_itstart = datetime.datetime.now()
            datet = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
            log.write('\n\n\n(hypoDD) \n=== ITERATION %3i (%3i) %s \n' % (iteri+1-jiter,iteri+1,datet))
            log.write('(hypoDD) Starting iteration %i with %i number of data points \n\n' % (iteri+1,ndt))

            """
            Get weighting parameters for this iteration:
            -----
            From the data weighting and re-weighting section of input file.
            """
            i=niter
            for i in range(0,niter):
                if iteri<aiter[i]:
                    break

            maxres_cross = amaxres_cross[i] # Max. residual threshold (corr)
            maxres_net = amaxres_net[i]     # Max. residual threshold (cat)
            maxdcc = amaxdcc[i]             # Max. interevent sep for linked pairs (corr)
            maxdct = amaxdct[i]             # Max. interevent sep for linked pairs (cat)
            wt_ccp = awt_ccp[i]             # Weight multiplier for cc_p (corr)
            wt_ccs = awt_ccs[i]             # Weight multiplier for cc_s (corr)
            wt_ctp = awt_ctp[i]             # Weight multiplier for ct_p (cat)
            wt_cts = awt_cts[i]             # Weight multiplier for ct_s (cat)
            damp = adamp[i]                 # Damping (only relevant for lsqr)

            # Write inversion parameters to log
            log.write('Weighting parameters for this iteration: \n')
            log.write('wt_ccp= %7.4f  wt_ccs= %7.4f  maxr_cc= %7.4f  maxd_cc= %7.4f \n' % 
                      (wt_ccp,wt_ccs,maxres_cross,maxdcc))
            log.write('wt_ctp= %7.4f  wt_cts= %7.4f  maxr_ct= %7.4f  maxd_ct= %7.4f  damp= %7.4f \n\n' % 
                      (wt_ctp,wt_cts,maxres_net,maxdct,damp))

            """
            Calculate travel times and slowness vectors
            -----
            Calls partials function from rt_functions.py
            """
            tempstart = datetime.datetime.now()
            [tmp_ttp,tmp_tts,tmp_xp,tmp_yp,tmp_zp,dist_all,az_all,
             ang_all] = partials(nsrc,src_cusp,src_lat,src_lon,src_dep,nsta,sta_lab,sta_lat,
                                 sta_lon,mod_nl,mod_ratio,mod_v,mod_top,fn_srcpar,return_all=True)
            """
            Calculate double difference vector and theor dt
            -----
            Calls dtres function from hypoDD_functions.py
            """
            [dt_cal,dt_res,tt] = dtres(reloctype,ndt,nsrc,dt_dt,dt_idx,dt_ista1,dt_ista2,dt_ic1,
                                       dt_ic2,src_t,tmp_ttp,tmp_tts)
            if fileout==0:
                """
                Save residual data for plotting
                ---
                Save CC and CT Data Separately
                Counts based on ncc and nct - will raise errors if ncc and nct
                are incorrect
                """
                if not iboot: # Only if not bootstrapping (too much memory otherwise)
                    """
                    Save CC Residuals
                    """
                    # Check to make sure ncc count is good
                    test = dt_idx[:ncc]
                    if not all(data<3 for data in test):
                        raise Exception('(hypoDD) NCC count incorrect. Fatal error.')
                    # Save dt_cal
                    fcal = 'CC_dtcal_%i_%i.txt' % (int(icl+1),int(iteri+1-jiter))
                    fcal = os.path.join(txtfol,fcal)
                    np.savetxt(fcal,dt_cal[:ncc])
                    # Save dt_dt
                    fdt = 'CC_dtdt_%i_%i.txt' % (int(icl+1),int(iteri+1-jiter))
                    fdt = os.path.join(txtfol,fdt)
                    np.savetxt(fdt,dt_dt[:ncc])
                    # Save traveltimes
                    ftt = 'CC_tt_%i_%i.txt' % (int(icl+1),int(iteri+1-jiter))
                    ftt = os.path.join(txtfol,ftt)
                    np.savetxt(ftt,tt[:ncc])

                    """
                    Save CT Residuals
                    """
                    # Check to make sure ncc count is good
                    test = dt_idx[ncc:ndt]
                    if not all(data>2 for data in test):
                        raise Exception('(hypoDD) NCT count incorrect. Fatal Error.')
                    if nct!=(ndt-ncc):
                        raise Exception('(hypoDD) NDT != NCC + NCT. Fatal Error.')
                    # Save dt_cal
                    fcal = 'CT_dtcal_%i_%i.txt' % (int(icl+1),int(iteri+1-jiter))
                    fcal = os.path.join(txtfol,fcal)
                    np.savetxt(fcal,dt_cal[ncc:ndt])
                    # Save dt_dt
                    fdt = 'CT_dtdt_%i_%i.txt' % (int(icl+1),int(iteri+1-jiter))
                    fdt = os.path.join(txtfol,fdt)
                    np.savetxt(fdt,dt_dt[ncc:ndt])
                    # Save traveltimes
                    ftt = 'CT_tt_%i_%i.txt' % (int(icl+1),int(iteri+1-jiter))
                    ftt = os.path.join(txtfol,ftt)
                    np.savetxt(ftt,tt[ncc:ndt])
            
            # Save first and last residuals for all runs
            # this is a bootstrapping output for evpair residuals
            if iteri==1:
                calstart = np.copy(dt_cal)
                dtdtstart = np.copy(dt_dt)
            elif iteri==maxiter:
                calend = np.copy(dt_cal)
                dtdtend = np.copy(dt_dt)
            """
            Get a priori weights and reweight residuals
            -----
            Calls weighting function from hypoDD_functions.py
            """
            [ineg,dt_wt,kiter] = weighting(log,reloctype,ndt,mbad,amcusp,idata,kiter,ineg,
                                           maxres_cross,maxres_net,maxdcc,maxdct,minwght,
                                           wt_ccp,wt_ccs,wt_ctp,wt_cts,dt_c1,dt_c2,dt_idx,
                                           dt_qual,dt_res,dt_offse,dt_offss)
            """
            Skip outliers and/or air quakes
            -----
            Data reorganization only
            Only triggered if events to be excluded or data to be excluded
            Calls skip function from hypoDD_functions.py
            """
            if ineg>0 or reloctype==2:
                log.write('(hypoDD) Skipping data triggered \n')
                [ndt,nev,nsrc,nsta,ev_cusp,ev_date,ev_time,ev_mag, ev_lat,ev_lon,ev_dep,
                 ev_x,ev_y,ev_z,ev_herr,ev_zerr,ev_res,src_cusp,src_lat,src_lon,src_dep,
                 src_lat0,src_lon0,src_x,src_y,src_z,src_t,src_x0,src_y0,src_z0,src_t0,
                 src_xi,src_yi,src_zi,src_ti,sta_lab,sta_lat,sta_lon,sta_rmsc,sta_rmsn,
                 sta_np,sta_ns,sta_nnp,sta_nns,dt_sta1,dt_sta2,dt_c1,dt_c2,dt_idx,dt_dt,
                 dt_qual,dt_cal,dt_ista1,dt_ista2,dt_ic1,dt_ic2,dt_res,dt_wt,dt_offse,
                 dt_offss,tmp_ttp,tmp_tts,tmp_xp,tmp_yp,tmp_zp,nct,
                 ncc] = skip(log,reloctype,kiter,minwght,ndt,nev,nsrc,nsta,ev_cusp,ev_date,
                             ev_time,ev_mag,ev_lat,ev_lon,ev_dep,ev_x,ev_y,ev_z,ev_herr,
                             ev_zerr,ev_res,src_cusp,src_lat,src_lon,src_dep,src_lat0,src_lon0,
                             src_x,src_y,src_z,src_t,src_x0,src_y0,src_z0,src_t0,src_xi,src_yi,
                             src_zi,src_ti,sta_lab,sta_lat,sta_lon,sta_rmsc,sta_rmsn,sta_np,
                             sta_ns,sta_nnp,sta_nns,dt_sta1,dt_sta2,dt_c1,dt_c2,dt_idx,dt_dt,
                             dt_qual,dt_cal,dt_ista1,dt_ista2,dt_ic1,dt_ic2,dt_res,dt_wt,dt_offse,
                             dt_offss,tmp_ttp,tmp_tts,tmp_xp,tmp_yp,tmp_zp,nct,ncc,amcusp,mbad)

                # Skip cluster if less than 2 events
                # i.e. if enough events are removed due to being airquakes
                if nev<2:
                    log.write('(hypoDD) Cluster has less than 2 events. \n')
                    print('(hypoDD) Cluster has less than 2 events.')
                    break
                log.write('(hypoDD) Number of data now %i \n\n' % ndt)
            else:
                log.write('(hypoDD) No data skipped. \n')

            """
            Get initial residual statistics (avrg,rms,var...)
            -----
            The inital run of resstat only triggered for iteration 1
            Else this function is called within inversion function
            Calls resstat from hypoDD_functions.py
            """
            if iteri==0:
                [rms_cc,rms_ct,rms_cc0,rms_ct0,rms_ccold,rms_ctold,rms_cc0old,rms_ct0old,
                 resvar1] = resstat(log,reloctype,idata,ndt,nev,dt_res,dt_wt,dt_idx,rms_cc,
                                    rms_ct,rms_cc0,rms_ct0,rms_ccold,rms_ctold,rms_cc0old,
                                    rms_ct0old)
            """
            Call inversion
            ---- 
            Functions location in hypoDD_inversion.py
            isolv indicated in input file
            ----
            SVD is an ok solution for small problems but we've noted some issues with shift
                biases and overfitting.
            We encourage all to use lsqr inversion regardless of problem size and also 
                incorporate a bootstrapping or other error analysis into uncertainty quantification.
            """
            #if iteri==5:
            #    import pdb; pdb.set_trace()
            log.write('\n(hypoDD) Pre-iversion iteration %i with %i number of data \n' % (iteri,ndt))
            if isolv==1:
                [src_cusp,src_dx,src_dy,src_dz,src_dt,src_ds,
                 src_ex,src_ey,src_ez,src_et,src_es,
                 exav,eyav,ezav,etav,esav,
                 dxav,dyav,dzav,dtav,dsav,
                 rms_cc,rms_ct,rms_cc0,rms_ct0,
                 rms_ccold,rms_ctold,rms_cc0old,rms_ct0old] = svd(log,reloctype,iteri,
                                                                  ndt,nev,nsrc,damp,mod_ratio,idata,
                                                                  ev_cusp,src_cusp,dt_res,dt_wt,
                                                                  dt_ista1,dt_ista2,dt_ic1,dt_ic2,
                                                                  exav,eyav,ezav,etav,esav,
                                                                  dxav,dyav,dzav,dtav,dsav,
                                                                  rms_cc,rms_ct,rms_cc0,rms_ct0,
                                                                  rms_ccold,rms_ctold,rms_cc0old,rms_ct0old,
                                                                  tmp_xp,tmp_yp,tmp_zp,dt_idx)
            elif isolv==2:
                #tempstart = datetime.datetime.now()
                [src_cusp,src_dx,src_dy,src_dz,src_dt,src_ds, 
                 src_ex,src_ey,src_ez,src_et,src_es,
                 exav,eyav,ezav,etav,esav,
                 dxav,dyav,dzav,dtav,dsav,
                 rms_cc,rms_ct,rms_cc0,rms_ct0,
                 rms_ccold,rms_ctold,rms_cc0old,rms_ct0old,
                 acond,normvar,sam] = lsqr(log,reloctype,iteri,ndt,nev,nsrc,damp,mod_ratio, 
                                             idata,ev_cusp,src_cusp,dt_res,dt_wt,
                                             dt_ista1,dt_ista2,dt_ic1,dt_ic2,
                                             exav,eyav,ezav,etav,esav,
                                             dxav,dyav,dzav,dtav,dsav,
                                             rms_cc,rms_ct,rms_cc0,rms_ct0,
                                             rms_ccold,rms_ctold,rms_cc0old,rms_ct0old,
                                             tmp_xp,tmp_yp,tmp_zp,dt_idx,normvar)
            else:
                raise Exception('(hypoDD) Error: invalid isolv value must be 1 or 2')
            log.write('\n(hypoDD) Post-iversion iteration %i with %i number of data \n' % (iteri,ndt))

            """
            Check for air quakes
            -----
            Air quakes are relocated events with a negative depth
            Negative depths indicate a depth above 0.0
            Air quake IDs are saved to amcusp and then removed before inversion
            in the next iteration
            """
            mbad = 0
            k = 0
            for i in range(nsrc):
                if (src_dep[i]+src_dz[i]/1000.)<0.: # If depth negative
                    log.write('(hypoDD) >>> Warning: negative depth - %12i \n' % ev_cusp[i])
                    amcusp[k]=ev_cusp[i]  # Save event ID
                    k += 1
                    if k>1000:  # Trigger exception if there are too many airquakes.
                                # Indicates poor data quality or poorly constrained inversion              
                        raise Exception('(hypoDD) >>> More than 1000 air quakes.')
            mbad = k # Number of neg depth events

            """
            Update iteration numbers:
            ----
            If there are airquakes the iteration DOES NOT COUNT
            If there are airquakes the event locations will not be updated (iskip)
            """
            iskip=0
            if mbad>0:
                for i in range(niter):
                    aiter[i] = aiter[i] + 1
                jiter+=1 # iteration with no update
                maxiter+=1
                normvar = np.append(normvar,np.zeros((1,3)),axis=0)

                log.write('(hypoDD) Number of air quakes (AQ) = %i \n' % mbad)
                if (nsrc-mbad)<=1:
                    log.write('(hypoDD) Warning: number of non-airquakes < 2.  Skipping cluster. \n')
                    print('(hypoDD) Warning: number of non-airquakes < 2.  Skipping cluster. \n')
                    continue
                iskip=1

            """
            Update source parameters:
            """
            if iskip==0:
                xav = 0 # mean centroid shift
                yav = 0
                zav = 0
                tav = 0
                sav = 0
                alon = 0
                alat = 0
                adep = 0

                if nsrc==1:
                    nsrc=nev

                # Save iteration updates to file
                # Open all files
                if fileout==0:
                    if not iboot:
                        fdstr = "delta_source_"+str(icl+1)+'_'+str(iteri+1-jiter)+".txt"
                        fds = os.path.join(txtfol,fdstr)
                        delta_source = open(fds,'w')
                        fastr = "abs_source_"+str(icl+1)+'_'+str(iteri+1-jiter)+".txt"
                        fas = os.path.join(txtfol,fastr)
                        abs_source = open(fas,'w')
                        finistr = "abs_source_"+str(icl+1)+".txt"
                        fini = os.path.join(txtfol,finistr)
                        abs_source0 = open(fini,'w')
                        ftstr = "tru_source_"+str(icl+1)+'_'+str(iteri+1-jiter)+".txt"
                        fts = os.path.join(txtfol,ftstr)
                        tru_source = open(fts,'w')

                # Loop over events
                for i in range(0,nsrc):
                    src_cusp[i] = ev_cusp[i]
                    # Update absolute source parameters
                    sensx = src_x
                    sensy = src_y
                    sensz = src_z
                    src_x[i] = src_x[i] + src_dx[i]
                    src_y[i] = src_y[i] + src_dy[i]
                    src_z[i] = src_z[i] + src_dz[i]
                    if reloctype==1:
                        src_t[i] = src_t[i] + src_dt[i]
                    if reloctype==2:
                        src_s[i] = src_s[i] + src_ds[i]
                    # Save updates and file locations to file
                    if fileout==0 and not iboot:
                        delta_source.write('%10.6f %10.6f %10.6f %10.6f %10.6f \n' % (src_dx[i],src_dy[i],src_dz[i],src_dt[i],src_ds[i]))
                        abs_source.write('%10.6f %10.6f %10.6f %10.6f %10.6f \n' % (src_x[i],src_y[i],src_z[i],src_t[i],src_s[i]))
                        abs_source0.write('%10.6f %10.6f %10.6f %10.6f %10.6f \n' % (src_x0[i],src_y0[i],src_z0[i],src_t0[i],src_t0[i]))
                        tru_source.write('%10.6f %10.6f %10.6f %10.6f \n' % (src_xi[i],src_yi[i],src_zi[i],src_ti[i]))
                        # Save delta and abs locations into arrays
                        locdel[:,i,iteri-jiter] = np.array([src_dx[i],src_dy[i],src_dz[i],src_dt[i],src_ds[i]])
                        locabs[:,i,iteri-jiter+1] = np.array([src_x[i],src_y[i],src_z[i],src_t[i],src_s[i]])
                        # Also save true locations
                        locabs[:,i,0] = np.array([src_x0[i],src_y0[i],src_z0[i],src_t0[i],src_t0[i]])
                        loctru[:,i] = np.array([src_xi[i],src_yi[i],src_zi[i],src_ti[i]])
                    # Update absolute source locations both in
                    # cartesian coordinates and in lat/lon
                    src_dep[i] = src_dep[i] + src_dz[i]/1000
                    [lat,lon] = sdc2(src_x[i]/1000,src_y[i]/1000,1)
                    src_lon[i] = lon
                    src_lat[i] = lat
                    alon = lon+alon
                    alat = lat+alat
                    adep = adep+src_dep[i]
                    # Calculate the centroid shift
                    ################## KB Note: Pull this out of loop maybe 
                    xav = xav + (src_x[i] - src_x0[i])
                    yav = yav + (src_y[i] - src_y0[i])
                    zav = zav + (src_z[i] - src_z0[i])
                    tav = tav + (src_t[i] - src_t0[i])
                    sav = sav + (src_s[i])
                #import pdb; pdb.set_trace()
                xav = xav/nsrc
                yav = yav/nsrc
                zav = zav/nsrc
                tav = tav/nsrc
                sav = sav/nsrc
                alon = alon/nsrc
                alat = alat/nsrc
                adep = adep/nsrc

                # Document centroid shift
                log.write('\n(hypoDD) Cluster centroid at: %10.6f  %11.6f  %9.6f \n' % (alat,alon,adep))
                log.write('(hypoDD) Mean centroid (origin) shift in x,y,z,t [m,ms]: %7.1f,%7.1f,%7.1f,%7.1f \n' % (xav,yav,zav,tav))
                # Close iteration files
                if fileout==0 and not iboot:
                    delta_source.close()
                    abs_source.close()
                    tru_source.close()

                # """
                # Parameter sensitivities
                # """
                # xchange = np.mean(np.abs((src_x-src_dx)/src_x))
                # ychange = np.mean(np.abs((src_y-src_dy)/src_y))
                # zchange = np.mean(np.abs((src_z-src_dz)/src_z))
                # print('\nParameter Sensitivities')
                # print('X: %f' % (1000*dt_change/xchange))
                # print('Y: %f' % (1000*dt_change/ychange))
                # print('Z: %f' % (1000*dt_change/zchange))
                # log.write('\nParameter Sensitivities\n')
                # log.write('X: %f \n' % (1000*dt_change/xchange))
                # log.write('Y: %f \n' % (1000*dt_change/ychange))
                # log.write('Z: %f \n' % (1000*dt_change/zchange))

                """
                Get interpair distance for each observation and average signal coherency
                """
                [ncc,nct,cohav,picav,
                 dt_offse] = sigcoherency(log,reloctype,nct,ncc,ndt,idata,src_x,src_y,src_z,
                                          dt_ic1,dt_ic2,dt_offse,dt_qual,dt_idx)

                """
                Get number of observations and mean residual at each station
                ----
                Resstat calculates residuals statistics per event
                This loop calculates residual statistics per station
                """
                # [sta_np,sta_ns,sta_nnp,sta_nns,
                #  sta_rmsc,sta_rmsn,tmpr1,tmpr2] = statres(log,nsta,ndt,idata,reloctype,
                #                                           sta_lab,dt_ista1,dt_ista2,
                #                                           dt_idx,dt_res)

                """
                Write output scratch mdat.reloc
                ----
                Separate files for each cluster and each iteration
                """
                if fileout==0:
                    i = iteri+1-jiter
                    str80 = '%s.%03i.%03i' % (fn_reloc,icl+1,i)
                    writereloc(str80,nev,src_cusp,src_lat,src_lon,src_dep,src_x,src_y,src_z,
                               src_ex,src_ey,src_ez,ev_date,ev_time,ev_mag,icl)
                    if (iteri+1)==maxiter:
                        str80 ='%s.%03i' % (fn_reloc,icl+1)
                        writereloc(str80,nev,src_cusp,src_lat,src_lon,src_dep,src_x,src_y,src_z,
                                   src_ex,src_ey,src_ez,ev_date,ev_time,ev_mag,icl)
                        str80 = fn_reloc
                        writereloc(str80,nev,src_cusp,src_lat,src_lon,src_dep,src_x,src_y,src_z,
                                   src_ex,src_ey,src_ez,ev_date,ev_time,ev_mag,icl,maxiter=True)   

            """
            Standard terminal outputs
            ------
            Only used if i/o turned on
            """
            if fileout==0:
                if isolv==1:
                    acond=-9
                terminaloutputs(reloctype,iteri,jiter,isolv,idata,
                                nev,nevold,nct,nctold,ncc,nccold,
                                rms_ct,rms_ctold,rms_cc,rms_ccold,tmpr1,tmpr2,
                                dxav,dyav,dzav,dtav,dsav,
                                xav,yav,zav,mbad,acond)

            datet = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
            log.write('\n\n(hypoDD) Iteration %2i finished with %i number of data at %s \n' % (iteri+1,ndt,datet))
            log.write('(hypoDD) Length of iteration: %s \n' % str(datetime.datetime.now()-time_itstart))
            log.write('(hypoDD) Time since start of cluster: %s \n' % str(datetime.datetime.now()-time_preit))
            #print('\n\n(hypoDD) Iteration ',iteri+1,' finished: ',datetime.datetime.now()-time_itstart)
            if iteri<maxiter:
                iteri+=1

        print('\n\n(hypoDD) Total number of iteration loops: ',iteri,' time spent: ',datetime.datetime.now()-time_preit)

        """
        Update origin times
        ----
        Only done on last iteration
        """
        if fileout==0:
            print('\n\n(hypoDD) Writing out results')
        for i in range(0,nev):
            src_t[i] = src_t[i]/1000. # change to s
            if src_t[i]>5:
                print('(hypoDD) Warning: org time diff > 5s for %i' % src_cusp[i])
            iyr = ev_date[i].year
            imo = ev_date[i].month
            idy = ev_date[i].day
            ihr = ev_time[i].hour
            imn = ev_time[i].minute
            isc = ev_time[i].second
            ims = ev_time[i].microsecond
            itf = datetime.datetime(iyr,imo,idy,ihr,imn,isc,ims)
            itf = itf - datetime.timedelta(seconds=src_t[i])
            ev_date[i] = datetime.date(itf.year,itf.month,itf.day)
            ev_time[i] = datetime.time(itf.hour,itf.minute,itf.second,itf.microsecond)

        """
        Get no. of obs per event and calculate rms statistics
        """
        [src_np,src_ns,src_nnp,src_nns,
         src_rmsc,src_rmsn] = eventstats(nev,ndt,dt_ic1,dt_ic2,dt_res,dt_idx)

        """
        Calculated RMS parameter sensitivities
        """
        mean_xp = np.mean(tmp_xp)
        mean_yp = np.mean(tmp_yp)
        mean_zp = np.mean(tmp_zp)
        rms_xp = np.sqrt(np.mean(tmp_xp*tmp_xp))
        rms_yp = np.sqrt(np.mean(tmp_yp*tmp_yp))
        rms_zp = np.sqrt(np.mean(tmp_zp*tmp_zp))

        """
        Output final residuals: mdat.res
        ---
        Includes final cluster only
        """
        if fileout==0:
            res = open(fn_res,'w')
            if reloctype==1:
                res.write('STA           DT        C1        C2    IDX     QUAL    RES [ms]   WT         OFFSE \n')
            if reloctype==2:
                res.write('STA1      STA2      DT        C1     IDX     QUAL    RES [ms]   WT         OFFSS \n')
            if reloctype==3:
                 res.write('STA1      STA2      DT        C1   C2      IDX     QUAL    RES [ms]   WT        OFFSE     OFFSS \n')               
            for j in range(0,ndt):
                if reloctype==1:
                    res.write('%7s %12.7f %9i %9i %1i %9.4f %12.6f %11.6f %8.1f \n' %
                              (dt_sta1[j],dt_dt[j],dt_c1[j],dt_c2[j],dt_idx[j],dt_qual[j],
                               dt_res[j]*1000,dt_wt[j],dt_offse[j]))
                if reloctype==2:
                    res.write('%7s %7s %12.7f %9i %1i %9.4f %12.6f %11.6f %8.1f \n' %
                              (dt_sta1[j],dt_sta2[j],dt_dt[j],dt_c1[j],dt_idx[j],dt_qual[j],
                               dt_res[j]*1000,dt_wt[j],dt_offss[j]))
                if reloctype==3:
                    res.write('%7s %7s %12.7f %9i %9i %1i %9.4f %12.6f %11.6f %8.1f %8.1f \n' %
                              (dt_sta1[j],dt_sta2[j],dt_dt[j],dt_c1[j],dt_c2[j],dt_idx[j],dt_qual[j],
                               dt_res[j]*1000,dt_wt[j],dt_offse[j],dt_offss[j]))
            res.close()
        #print('(hypoDD) Calc time cluster: ',datetime.datetime.now()-time_cluststart)


        """
        Output final locations: mdat.reloc
        """
        #if fileout==0:
        #    writereloc(fn_reloc,nev,src_cusp,src_lat,src_lon,src_dep,src_x,src_y,src_z,
        #               src_ex,src_ey,src_ez,ev_date,ev_time,ev_mag)

        """
        Output stations: mdat.station
        """
        if fileout==0:
            if len(fn_stares)>1:
                stares = open(fn_stares,'a')
                for i in range(0,nsta):
                    stares.write('%7s %9.4f %9.4f %7i %7i %7i %7i %9.4f %9.4f %3i \n' %
                                 (sta_lab[i],sta_lat[i],sta_lon[i],#sta_dist[i],sta_az[i],
                                 sta_np[i],sta_ns[i],sta_nnp[i],sta_nns[i],
                                 sta_rmsc[i],sta_rmsn[i],icl+1))
                stares.close()

    """
    Traveltime and Absolute Location Statistics
    -----
    ########## KB Note: Added to incorporate and understand
                        absolute location information.
    -----
    Calculates mean, med, RMS for location changes for events.
    Calculates mean, med, RMS for traveltime residuals of final event locations/times.
    """
    #if fileout==0:
    # Mean dx,dy,dz,dt
    src_dx_tot = np.abs(src_x-src_xi)
    src_dy_tot = np.abs(src_y-src_yi)
    src_dz_tot = np.abs(src_z-src_zi)
    src_dt_tot = np.abs(src_dt-src_ti)
    #mean_dx = np.mean(np.abs(src_dx_tot))
    #mean_dy = np.mean(np.abs(src_dy_tot))
    #mean_dz = np.mean(np.abs(src_dz_tot))
    #mean_dt = np.mean(np.abs(src_dt_tot))
    mean_dx = np.mean(src_dx_tot)
    mean_dy = np.mean(src_dy_tot)
    mean_dz = np.mean(src_dz_tot)
    mean_dt = np.mean(src_dt_tot)
    # Median
    med_dx = np.median(src_dx_tot)
    med_dy = np.median(src_dy_tot)
    med_dz = np.median(src_dz_tot)
    med_dt = np.median(src_dt_tot)
    # RMS
    rms_dx = np.sqrt(np.mean(src_dx_tot[:]**2))
    rms_dy = np.sqrt(np.mean(src_dy_tot[:]**2))
    rms_dz = np.sqrt(np.mean(src_dz_tot[:]**2))
    rms_dt = np.sqrt(np.mean(src_dt_tot[:]**2))
    #rms_dx = np.sqrt(np.sum(src_dx_tot[:]*src_dx_tot[:])/nsrc)
    #rms_dy = np.sqrt(np.sum(src_dy_tot[:]*src_dy_tot[:])/nsrc)
    #rms_dz = np.sqrt(np.sum(src_dz_tot[:]*src_dz_tot[:])/nsrc)
    #rms_dt = np.sqrt(np.sum(src_dt_tot[:]*src_dt_tot[:])/nsrc)
    # Print these:
    print('\n\n Relocation Stats for Cluster ',(icl+1))
    print('        DX (m)       DY (m)       DZ (m)       DT (ms)')
    print('Mean   %3.6f %3.6f %3.6f %3.6f' % (mean_dx,mean_dy,mean_dz,mean_dt))
    print('Median %3.6f %3.6f %3.6f %3.6f' % (med_dx,med_dy,med_dz,med_dt))
    print('RMS    %3.6f %3.6f %3.6f %3.6f' % (rms_dx,rms_dy,rms_dz,rms_dt))

    # Now do the traveltime residual statistics
    print('(hypoDD) writing update traveltime residuals to tt_res_final.txt')
    log.write('\n\n (hypoDD) writing update traveltime residuals to tt_res_final.txt \n')
    tmp_ttp_old = tmp_ttp
    tmp_tts_old = tmp_tts
    [tmp_ttp,tmp_tts,
     tmp_xp,tmp_yp,tmp_zp] = partials(nsrc,src_cusp,src_lat,src_lon,src_dep,
                                      nsta,sta_lab,sta_lat,sta_lon,
                                      mod_nl,mod_ratio,mod_v,mod_top,fn_srcpar)
    tmp_ttp_res = tmp_ttp_old-tmp_ttp
    tmp_tts_res = tmp_tts_old-tmp_tts

    # fig,ax = plt.subplots(ncols=2,nrows=1,figsize=(8,3),sharex=True,sharey=True)
    # ax[0,0].hist(tmp_ttp_res.ravel(),50,histtype='step')
    # ax[0,0].set_title('P TT Residuals')
    # ax[0,0].grid('True',linestyle='--',alpha=0.5)
    # ax[0,0].tick_params(right=True,top=True)
    # ax[0,0].set_ylabel('Counts')
    # ax[0,0].set_xlabel('Residual (s)')
    # ax[0,1].hist(tmp_tts_res.ravel(),50,histtype='step')
    # ax[0,1].set_title('S TT Residuals')
    # ax[0,1].grid('True',linestyle='--',alpha=0.5)
    # ax[0,1].tick_params(right=True,top=True)
    # ax[0,1].set_xlabel('Residual (s)')


    if not iboot:
        resstr = os.path.join(txtfol,'tt_res_final_icl.txt')
        res = open(resstr,'w')
        for i in range(0,nsrc):
            for j in range(0,nsta):
                res.write('%i %s P %3.8f %3.8f %3.8f' % (src_cusp[i],sta_lab[j],tmp_ttp_old[j,i],tmp_ttp[j,i],tmp_ttp_res[j,i]))
                res.write('%i %s S %3.8f %3.8f %3.8f' % (src_cusp[i],sta_lab[j],tmp_tts_old[j,i],tmp_tts[j,i],tmp_tts_res[j,i]))
        res.close()
    print('Traveltime Residuals P for Cluster ',(icl+1))
    print('Mean: %3.6f, Median: %3.6f, RMS: %3.6f' % (tmp_ttp_res.mean(axis=(0,1)),np.median(tmp_ttp_res,axis={0,1}),np.sqrt(np.sum(tmp_ttp_res[:,:]*tmp_ttp_res[:,:],axis=(0,1))/nsrc)))
    print('Traveltime Residuals S')
    print('Mean: %3.6f, Median: %3.6f, RMS: %3.6f' % (tmp_tts_res.mean(axis=(0,1)),np.median(tmp_tts_res,axis={0,1}),np.sqrt(np.sum(tmp_tts_res[:,:]*tmp_tts_res[:,:],axis=(0,1))/nsrc)))
    print('\n\n')
    log.write('\n\n Traveltime Residuals P for Cluster %i \n' % (icl+1))
    log.write('Mean: %3.6f, Median: %3.6f, RMS: %3.6f \n' % (tmp_ttp_res.mean(axis=(0,1)),np.median(tmp_ttp_res,axis={0,1}),np.sqrt(np.sum(tmp_ttp_res[:,:]*tmp_ttp_res[:,:],axis=(0,1))/nsrc)))
    log.write('Traveltime Residuals S \n')
    log.write('Mean: %3.6f, Median: %3.6f, RMS: %3.6f \n' % (tmp_tts_res.mean(axis=(0,1)),np.median(tmp_tts_res,axis={0,1}),np.sqrt(np.sum(tmp_tts_res[:,:]*tmp_tts_res[:,:],axis=(0,1))/nsrc)))

    """
    Parameter sensitivities
    """
    print('Parameter Sensitivities Cluster: ',(icl+1))
    print('Mean X: ',np.mean(dt_dt)/np.mean(src_dx_tot)) 
    print('Mean Y: ',np.mean(dt_dt)/np.mean(src_dy_tot)) 
    print('Mean Z: ',np.mean(dt_dt)/np.mean(src_dz_tot)) 

    datet = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    log.write('HypoDD Done %s \n' % datet)
    print('\nHypoDD Done %s \n\n' % datet)


    # If fileout is on (0) then outputs are located in files.
    # If fileout if off (1) then outputs are returned.
    if fileout==0:
        return None
    elif fileout==1:
        return calstart,calend,dtdtstart,dtdtend,locdel,locabs,loctru
