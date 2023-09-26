import numpy as np
import os
import sys
import datetime

# Import MethodType Data Reading Functions
from methodTypes.eventPair.hypoDD_ev import readccs_evpair,readcts_evpair
from methodTypes.stationPair.hypoDD_stat import readccs_stpair,readcts_stpair
from methodTypes.doublePair.hypoDD_doub import readccs_doubpair,readcts_doubpair

"""
This script contains hypoDD subfunctions for reading in files or
for writing to file if needed.
----
This script includes the following functions:
    hypoDD_input:   Read in to memory the hypoDD input file (hypoDD.inp)
    readevent:      Read in the event.dat file
    readcts:        Read in catalog data
    readccs:        Read in cc data
    readdata:       Call readcts and/or readccs and organize data arrays
    writeloc:       Write hypoDD.loc file
    writereloc:     Write hypoDD.reloc file
----
"""

def hypoDD_input(log='',fn_inp='inpfiles/hypoDD.inp',fileout=0):
    """
    Get input parameters from the hypoDD.inp file
    :::
    Inputs:
    log (file obj)  --- Opened file location for log file
    fn_inp (str)    --- File location for input file
    fileout (int)   --- Integer switch to indicate file i/o [default=0]
    :::
    Outputs:
    retlist (list)  --- List object containing all input parameters 
                        from hypoDD.inp file
    :::                
    """

    """
    Read in input file
    """
    # Open input file
    inputfile = open(fn_inp,'r')
    # Define counter variables
    ncusp = 0   # Number of events
    niter = 0   # Number of iteration blocks
    l = 0       # Line number
    ii = 0      # Number of listed events
    icusp = []  # Empty event list

    # Loop to read each parameter line, skipping comments
    inputs = inputfile.readlines()
    for line in inputs:
        # Skip comments
        if line[0:1] == '*' or line[1:2] == '*':
            continue

        line = line.strip()
        line = line.split()
        line = list(filter(None,line))

        """
        Input and Output File Locations
        """
        if l==0:
            fn_cc = str(line[0])        # CC Catalog
        if l==1:
            fn_ct = str(line[0])        # Dt Catalog
        if l==2:
            fn_eve = str(line[0])       # Event File
        if l==3:
            fn_sta = str(line[0])       # Station File
        if l==4:
            fn_loc = str(line[0])       # Location Output
        if l==5:    
            fn_reloc = str(line[0])     # Relocation Output
        if l==6:
            fn_stares = str(line[0])    # Station Residual Output
        if l==7:
            fn_res = str(line[0])       # Event-pair Residual Output
        if l==8:
            fn_srcpar = str(line[0])    # Source Pair Information

        """
        Data and Clustering Information
        """
        if l==9:
            idata = int(line[0])        # Type of Data
            iphase = int(line[1])       # Phases Used
            maxdist = float(line[2])    # Max. Event-Station Separation
        if l==10:
            minobs_cc = int(line[0])    # Min. # of Shared CC Measurements to Pair
            minobs_ct = int(line[1])    # Min. # of Shared CT Measurements to Pair

        """
        Iteration Variables
        """
        if l==11:
            istart = int(line[0])       # Start from Catalor or Cluster Centroid
            isolv = int(line[1])        # LSQR or SVD Inversion
            niter = int(line[2])        # Number of Inversion Iterations

            # Initialize Inversion Iteration Variable Arrays
            # Need to be the length of the number of sets of iterations
            aiter = np.zeros(niter)
            awt_ccp = np.zeros(niter)
            awt_ccs = np.zeros(niter)
            amaxres_cross = np.zeros(niter)
            amaxdcc = np.zeros(niter)
            awt_ctp = np.zeros(niter)
            awt_cts = np.zeros(niter)
            amaxres_net = np.zeros(niter)
            amaxdct = np.zeros(niter)
            adamp = np.zeros(niter)

        """
        Iteration Weighting Variables
        """
        if l >= 12 and l <= 11+niter:
            i = l-12
            aiter[i] = int(line[0])             # Number of iterations with these variables
            awt_ccp[i] = float(line[1])         # P CC Data Weighting
            awt_ccs[i] = float(line[2])         # S CC Data Weighting
            amaxres_cross[i] = float(line[3])   # Cross-Correlation Res. Cutoff
            amaxdcc[i] = float(line[4])         # Cross-Correlation Inter-event Dist. Cutoff
            awt_ctp[i] = float(line[5])         # P CT Data Weighting
            awt_cts[i] = float(line[6])         # S CT Data Weighting
            amaxres_net[i] = float(line[7])     # CT Res. Cutoff
            amaxdct[i] = float(line[8])         # CT Inter-event Dist. Cutoff
            adamp[i] = float(line[9])           # Damping Value for LSQR Inversion Only

        """
        Velocity Model Information
        """
        if l==12+niter:
            mod_nl = int(line[0])           # Number of Layers
            mod_ratio = float(line[1])      # VPVS Ratio
            mod_top = np.zeros(mod_nl)
            mod_v = np.zeros(mod_nl)
        if l==13+niter:
            for layer in range(0,mod_nl):
                mod_top[layer] = float(line[layer])     # Depth to Top of Layer in KM
        if l==14+niter:
            for layer in range(0,mod_nl):
                mod_v[layer] = float(line[layer])       # P Velocity of Layer in KM/S

        """
        Define Clusters
        ---
        If relocating certain cluster, this is used to define cluster ID
        Else iclust==0 for all clusters
        """
        if l==15+niter:
            iclust = int(line[0])

        """
        Define Event IDS
        ---
        If relocating certain events this is used
        Else this is an empty list
        """
        if l>=16+niter:
            line = [int(eid) for eid in line]
            icusp.extend(line)
        
        l = l+1

    inputfile.close()

    """
    Edit Variables as Needed
    """
    
    #If icusp is defined (i.e. only relocate certain events),
    #convert the list to numpy
    if len(icusp)>0:
        icusp = np.asarray(icusp,dtype='int')
    # Update the total number of events to the subsample of events number
    ncusp = len(icusp)

    # Rearange aiter numbers (i.e. aiter[i] as starting iteration with those 
    # variables rather than number of that iteration type)
    maxiter = np.sum(aiter)
    tmp = 0
    aitercopy = np.zeros(niter)
    for i in range(1,niter):
        tmp += aiter[i-1]
        aitercopy[i] = tmp
    aitercopy[0] = 0
    aiter = aitercopy

    # # Check event locations if file i/o needed
    # if fileout==0:
    #     if not os.path.exists(fn_sta):
    #         raise RuntimeError('File does not exist ' + fn_sta)
    #     if not os.path.exists(fn_eve):
    #         print('File does not exist ' + fn_eve)
    #     #if idata==1 or idata==3:
    #     #    if not os.path.exists(fn_cc):
    #     #        raise RuntimeError('File does not exist ' + fn_cc)
    #     if idata==2 or idata==3:
    #         if not os.path.exists(fn_inp):
    #             raise RuntimeError('File does not exist ' + fn_ct)  
    
    # Set max number of iterations
    #maxiter = aiter[niter-1]
    # Set first iteration start to 0
    #aiter[0] = 0 

    """
    Write to Log
    """
    if log!='':
        # Write files and variables to log file
        if fileout==0:
            log.write('INPUT FILES: \ncross dtime data: %s \ncatalog dtime data: %s \nevents: %s \nstations: %s \n\n'
                       % (fn_cc,fn_ct,fn_eve,fn_sta))
            log.write('OUTPUT FILES: \ninitial locations: %s \nrelocated events: %s \nevent pair residuals: %s \nstation residuals: %s \nsource parameters: %s \n\n'
                       % (fn_loc,fn_reloc,fn_res,fn_stares,fn_srcpar))
        log.write('\nINPUT PARAMETERS: \nIDATA: %i \nIPHASE: %i \nMAXDIST = %2.4f \nMINOBS_CC: %i \nMINOBS_CT = %i \nISTART: %i \nISOLV: %i \n\n'
                   % (idata,iphase,maxdist,minobs_cc,minobs_ct,istart,isolv))
        # Write inversion variables to log file
        if niter==1:
            log.write('ITER: 1 - %i \nDAMP: %i \nWT_CCP: %2.4f \nWT_CCS: %2.4f \nMAXR_CC: %2.4f \nMAXD_CC: %2.4f \nWT_CTP: %2.4f \nWT_CTS: %2.4f \nMAXR_CT %2.4f \nMAXD_CT: %2.4f \n\n' % (maxiter,adamp[i],awt_ccp[i],awt_ccs[i],amaxres_cross[i],amaxdcc[i],awt_ctp[i],awt_cts[i],amaxres_net[i],amaxdct[i]))
        else:
            for i in range(niter):
                if i==0:    # First set of iterations
                    log.write('ITER: 1 - %i \nDAMP: %i \nWT_CCP: %2.4f \nWT_CCS: %2.4f \nMAXR_CC: %2.4f \nMAXD_CC: %2.4f \nWT_CTP: %2.4f \nWT_CTS: %2.4f \nMAXR_CT %2.4f \nMAXD_CT: %2.4f \n\n' % (aiter[i+1],adamp[i],awt_ccp[i],awt_ccs[i],amaxres_cross[i],amaxdcc[i],awt_ctp[i],awt_cts[i],amaxres_net[i],amaxdct[i]))
                elif i==niter-1:    # Last set of iterations
                    log.write('ITER: %i - %i \nDAMP: %i \nWT_CCP: %2.4f \nWT_CCS: %2.4f \nMAXR_CC: %2.4f \nMAXD_CC: %2.4f \nWT_CTP: %2.4f \nWT_CTS: %2.4f \nMAXR_CT %2.4f \nMAXD_CT: %2.4f \n\n' % (aiter[i]+1,maxiter,adamp[i],awt_ccp[i],awt_ccs[i],amaxres_cross[i],amaxdcc[i],awt_ctp[i],awt_cts[i],amaxres_net[i],amaxdct[i]))
                else:
                    log.write('ITER: %i - %i \nDAMP: %i \nWT_CCP: %2.4f \nWT_CCS: %2.4f \nMAXR_CC: %2.4f \nMAXD_CC: %2.4f \nWT_CTP: %2.4f \nWT_CTS: %2.4f \nMAXR_CT %2.4f \nMAXD_CT: %2.4f \n\n' % (aiter[i]+1,aiter[i+1],adamp[i],awt_ccp[i],awt_ccs[i],amaxres_cross[i],amaxdcc[i],awt_ctp[i],awt_cts[i],amaxres_net[i],amaxdct[i]))
        # Write velocity model to log file
        log.write('MOD_NL: %i \nMOD_RATIO: %2.4f \n' % (mod_nl,mod_ratio))
        log.write('MOD_TOP       MOD_V \n')
        for i in range(0,mod_nl):
            log.write('%4.5f     %4.5f \n' % (mod_top[i],mod_v[i]))
        # If defined:
        # Repeat number of clusters/events to relocate
        if iclust==0:
            log.write('Relocate all cluster. \n\n')
            if fileout == 0:
                print('Relocate all clusters.')
        else:
            log.write('Relocate cluster number %i \n\n' % iclust)
            if fileout == 0:
                print('Relocate cluster number %i' % iclust)
        if ncusp==0:
            log.write('Relocate all events. \n\n')
            if fileout == 0:
                print('Relocate all events.')
        else:
            log.write('Relocate %i events \n\n' % ncusp)
            if fileout == 0:
                print('Relocate %i events' % ncusp)

    """
    Return input variables
    """
    retlist = [fn_cc,fn_ct,fn_sta,fn_eve,fn_loc,fn_reloc,fn_res,fn_stares,fn_srcpar,
               idata,iphase,minobs_cc,minobs_ct,amaxres_cross,amaxres_net,amaxdcc,amaxdct,
               maxdist,awt_ccp,awt_ccs,awt_ctp,awt_cts,adamp,istart,maxiter,
               isolv,niter,aiter,mod_nl,mod_ratio,mod_v,mod_top,iclust,ncusp,icusp]
    return retlist


def readevent(log,fn_eve,ncusp,icusp):
    """
    Read the events in the event.dat file into event arrays
    :::
    Parameters:
    log (file object) --- Log file
    fn_eve (str) --- File location for event.dat file
    :::
    Returns:
    retlist (list) --- List of event information and arrays
    :::
    """

    """
    Read event file
    """
    # Open and Read Event File
    evfile = open(fn_eve,'r')
    evs = evfile.readlines()

    # Define Loop Variables
    i = int(1)
    nev = int(len(evs))

    # Define Event Data Arrays
    ev_date = np.zeros(nev,dtype='object')  # Event Dates
    ev_time = np.zeros(nev,dtype='object')  # Event Times
    ev_lat = np.zeros(nev,dtype='float')    # Event Latitudes
    ev_lon = np.zeros(nev,dtype='float')    # Event Longitudes
    ev_dep = np.zeros(nev,dtype='float')    # Event Depths
    ev_mag = np.zeros(nev,dtype='float')    # Event Magnitudes
    ev_herr = np.zeros(nev,dtype='float')   # Event Horizontal Location Errors
    ev_zerr = np.zeros(nev,dtype='float')   # Event Vertical Location Errors
    ev_res = np.zeros(nev,dtype='float')    # Event Residuals
    ev_cusp = np.zeros(nev,dtype='int')     # Event Integer IDs

    """
    EQ read loop
    """
    for ind,ev in enumerate(evs):
        ev = ev.split()
        ev = list(filter(None,ev))
        ev_date[ind] = datetime.date(year=int(ev[0][0:4]),month=int(ev[0][4:6]),day=int(ev[0][6:]))
        ev_time[ind] = datetime.time(hour=int(ev[1][0:2]),minute=int(ev[1][2:4]),second=int(ev[1][4:6]),microsecond=int(ev[1][6:])*10000)
        ev_lat[ind] = float(ev[2])
        ev_lon[ind] = float(ev[3])
        ev_dep[ind] = float(ev[4])
        ev_mag[ind] = float(ev[5])
        ev_herr[ind] = float(ev[6])
        ev_zerr[ind] = float(ev[7])
        ev_res[ind] = float(ev[8])
        ev_cusp[ind] = int(ev[9])

        """
        Depth Correction
        ---
        If EQ is shallower than 10m, correct depth to 10m.
        The program will break for air quakes later so it is important for events
        to start below the surface.
        """
        if ev_dep[ind] < 0.01:
            ev_dep[ind] = 0.01

    evfile.close()

    """
    Check events to event list from hypodd.inp
    ---
    Only if hypoDD list is populated
    """
    count=int(0)    # New counter Variable
    if ncusp>0:
        for k in range(0,nev):          # Loop over all events in event file
            if ev_cusp[k] in icusp:     # Save i
                ev_date[count] = ev_date[k]
                ev_time[count] = ev_time[k]
                ev_lat[count] = ev_lat[k]
                ev_lon[count] = ev_lon[k]
                ev_dep[count] = ev_dep[k]
                ev_mag[count] = ev_mag[k]
                ev_herr[count] = ev_herr[k]
                ev_zerr[count] = ev_zerr[k]
                ev_res[count] = ev_res[k]
                ev_cusp[count] = ev_cusp[k]
                count = count+1
        # Cut Event Arrays to Number of Events Only
        ev_date = ev_date[0:ncusp]
        ev_time = ev_time[0:ncusp]
        ev_lat = ev_lat[0:ncusp]
        ev_lon = ev_lon[0:ncusp]
        ev_dep = ev_dep[0:ncusp]
        ev_mag = ev_mag[0:ncusp]
        ev_herr = ev_herr[0:ncusp]
        ev_zerr = ev_zerr[0:ncusp]
        ev_res = ev_res[0:ncusp]
        ev_cusp = ev_cusp[0:ncusp]
        nev = ncusp     # Reset number of events

    """
    Write to log
    """
    log.write('# events = %i \n' % nev)
    #print('# events = %i' % nev)

    """
    Raise exception if mismatch in arrays
    ---
    If there is a mismatch in the ncusp and nevs, there is an error.
    """
    if ncusp>0 and ncusp!=nev:
        log.write('>>> Events repeated in selection list or missing/repeated in event file.\n')
        print('>>> Events repeated in selection list or missing/repeated in event file.')
        for i in range(0,ncusp):
            k = 0
            for j in range(0,nev):
                if icusp[i]==ev_cusp[j]:
                    k+=1
            if k==0:    # Event is missing from event file.
                log.write('%i is missing \n' % icusp[i])
                print('%i is missing' % icusp[i])
            if k>=2:    # Event is repeated in event or selection list.
                log.write('%i is non-unique \n' % icusp[i])
                raise Exception('Event ID must be unique %i' % icusp[i])

    """
    Return Event Information
    """
    retlist = [nev,ev_date,ev_time,ev_lat,ev_lon,ev_dep,ev_mag,
               ev_herr,ev_zerr,ev_res,ev_cusp]
    return retlist


def readcts(log,reloctype,fn_ct,ncc,maxsep_ct,iphase,icusp,iicusp,
            ev_lat,ev_lon,ev_dep,sta_lab,sta_lat,sta_lon,
            dt_sta1,dt_sta2,dt_dt,dt_qual,dt_offse,dt_offss,
            dt_c1,dt_c2,dt_idx):
    """
    This function reads in the dt.ct file.
    :::
    Parameters:
    log (file object) --- Log file
    reloctype (int) --- Double-difference pairing type
    fn_ct (str) --- Catalog file location
    ncc (int) --- Number of cc values (for array indexing)
    maxsep_ct (float) --- Max. interevent sep. for ct data
    iphase (int) --- Phase switch
    icusp[nev] (int array) --- Sorted event indexes
    iicusp[nev] (int array) --- Sorted event IDs
    ev_lat[nev] (float array) --- Event latitudes
    ev_lon[nev] (float array) --- Event longitudes
    ev_dep[nev] (float array) --- Event depths
    sta_lab[nsta] (object array) --- Station codes
    sta_lat[nsta] (float array) --- Station latitudes
    sta_lon[nsta] (float array) --- Station longitudes
    dt_sta1[ndt] (object array) --- Data station 1 codes
    dt_sta2[ndt] (object array) --- Data station 2 codes
    dt_dt[ndt] (float array) --- Data diff. time array
    dt_qual[ndt] (float array) --- Data weight array
    dt_offse[ndt] (float array) --- Data interevent offsets
    dt_offss[ndt] (float array) --- Data interstation offsets
    dt_c1[ndt] (int array) --- Data event 1 IDs
    dt_c2[ndt] (int array) --- Data event 2 IDs
    dt_idx[ndt] (int array) --- Data type indexes
    :::
    Returns:
    retlist (list) --- List of updated data arrays
    :::
    """

    if reloctype==1:
        [nctp,ncts,dt_sta1,dt_dt,dt_qual,
         dt_offse,dt_c1,dt_c2,dt_idx] = readcts_evpair(log,fn_ct,ncc,maxsep_ct,
                                                       iphase,icusp,iicusp,
                                                       ev_lat,ev_lon,ev_dep,sta_lab,
                                                       dt_sta1,dt_dt,dt_qual,dt_offse,
                                                       dt_c1,dt_c2,dt_idx)
    elif reloctype==2:
        [nctp,ncts,dt_sta1,dt_sta2,dt_dt,
         dt_qual,dt_offss,dt_c1,dt_idx] = readcts_stpair(log,fn_ct,ncc,maxsep_ct,
                                                         iphase,icusp,iicusp,
                                                         ev_lat,ev_lon,ev_dep,
                                                         sta_lab,sta_lat,sta_lon,
                                                         dt_sta1,dt_sta2,dt_dt,dt_qual,
                                                         dt_offss,dt_c1,dt_idx)
    elif reloctype==3:
        [nctp,ncts,dt_sta1,dt_sta2,dt_dt,
         dt_qual,dt_offse,dt_offss,
         dt_c1,dt_c2,dt_idx] = readcts_doubpair(log,fn_ct,ncc,maxsep_ct,
                                                iphase,icusp,iicusp,
                                                ev_lat,ev_lon,ev_dep,
                                                sta_lab,sta_lat,sta_lon,
                                                dt_sta1,dt_sta2,dt_dt,dt_qual,
                                                dt_offse,dt_offss,dt_c1,dt_c2,dt_idx)

    """
    Record catalog dt measurements to log
    """
    if iphase!=2:
        log.write('# Catalog P dtimes = %7i \n' % (nctp))
        #print('# Catalog P dtimes = %7i' % (nctp))
    if iphase!=1:
        log.write('# Catalog S dtimes = %7i \n' % (ncts))
        #print('# Catalog S dtimes = %7i' % (ncts))
    """
    Return updated arrays
    """
    retlist = [nctp,ncts,dt_sta1,dt_sta2,dt_dt,dt_qual,dt_offse,dt_offss,dt_c1,dt_c2,dt_idx]
    return retlist


def readccs(log,reloctype,fn_cc,maxsep_cc,iphase,icusp,iicusp,
            ev_lat,ev_lon,ev_dep,sta_lab,sta_lat,sta_lon,
            dt_sta1,dt_sta2,dt_dt,dt_qual,dt_offse,dt_offss,dt_c1,dt_c2,dt_idx):
    """
    This function reads in the dt.cc file.
    :::
    Parameters:
    log (file object) --- Log file
    reloctype (int) --- Double-difference pairing type
    fn_cc (str) --- CC file location
    maxsep_cc (float) --- Max. interevent sep. for cc data
    iphase (int) --- Phase switch
    icusp[nev] (int array) --- Sorted event indexes
    iicusp[nev] (int array) --- Sorted event IDs
    ev_lat[nev] (float array) --- Event latitudes
    ev_lon[nev] (float array) --- Event longitudes
    ev_dep[nev] (float array) --- Event depths
    sta_lab[nsta] (object array) --- Station codes
    sta_lat[nsta] (float array) --- Station latitudes
    sta_lon[nsta] (float array) --- Station longitudes
    dt_sta1[ndt] (object array) --- Data station 1 codes
    dt_sta2[ndt] (object array) --- Data station 2 codes
    dt_dt[ndt] (float array) --- Data diff. time array
    dt_qual[ndt] (float array) --- Data weight array
    dt_offse[ndt] (float array) --- Data interevent offsets
    dt_offss[ndt] (float array) --- Data interstation offsets
    dt_c1[ndt] (int array) --- Data event 1 IDs
    dt_c2[ndt] (int array) --- Data event 2 IDs
    dt_idx[ndt] (int array) --- Data type indexes
    :::
    Returns:
    retlist (list) --- List of updated data arrays
    :::
    """

    if reloctype==1:
        [nccp,nccs,dt_sta1,dt_dt,dt_qual,
         dt_offse,dt_c1,dt_c2,
         dt_idx,iiotc] = readccs_evpair(log,fn_cc,maxsep_cc,iphase,icusp,iicusp,
                                        ev_lat,ev_lon,ev_dep,sta_lab,
                                        dt_sta1,dt_dt,dt_qual,dt_offse,
                                        dt_c1,dt_c2,dt_idx)
    elif reloctype==2:
        [nccp,nccs,dt_sta1,dt_sta2,dt_dt,
         dt_qual,dt_offss,dt_c1,
         dt_idx,iiotc] = readccs_stpair(log,fn_cc,maxsep_cc,iphase,icusp,iicusp,
                                        ev_lat,ev_lon,ev_dep,sta_lab,sta_lat,sta_lon,
                                        dt_sta1,dt_sta2,dt_dt,dt_qual,dt_offss,dt_c1,dt_idx)
    elif reloctype==3:
        [nccp,nccs,dt_sta1,dt_sta2,dt_dt,
         dt_qual,dt_offse,dt_offss,dt_c1,dt_c2,
         dt_idx,iiotc] = readccs_doubpair(log,fn_cc,maxsep_cc,iphase,icusp,iicusp,
                                          ev_lat,ev_lon,ev_dep,sta_lab,sta_lat,sta_lon,
                                          dt_sta1,dt_sta2,dt_dt,dt_qual,
                                          dt_offse,dt_offss,dt_c1,dt_c2,dt_idx)

    """
    Record to log
    ---
    Record the number of cc values read into memory
    """
    if iphase!=2:
        log.write('# Cross-correlation P dtimes = %7i (no OTC for: %7i) \n' % (nccp,iiotc))
        #print('# Cross-correlation P dtimes = %7i (no OTC for: %7i)' % (nccp,iiotc))
    if iphase!=1:
        log.write('# Cross-correlation S dtimes = %7i (no OTC forL %7i) \n' % (nccs,iiotc))
        #print('# Cross-correlation S dtimes = %7i (no OTC for: %7i)' % (nccs,iiotc))

    """
    Return the read in data
    """
    retlist = [nccp,nccs,dt_sta1,dt_sta2,dt_dt,dt_qual,
               dt_offse,dt_offss,dt_c1,dt_c2,dt_idx]
    return retlist 


def readdata(log,reloctype,idata,iphase,fn_cc,fn_ct,nev,ev_cusp,ev_lat,ev_lon,ev_dep,
             nsta,sta_lab,sta_lat,sta_lon,maxsep_ct,maxsep_cc):
    """
    The read data function read in and returns data arrays
    from the dt.ct and dt.cc files.
    :::
    Parameters:
    log (file object) --- Log file
    idata (int) --- Data type switch
    iphase (int) --- Phase switch
    fn_cc (str) --- Cross-correlation file location
    fn_ct (str) --- Catalog file location
    nev (int) --- No. of events
    ev_cusp[nev] (int array) --- Event ID array
    ev_lat[nev] (float array) --- Event latitude
    ev_lon[nev] (float array) --- Event longitude
    ev_dep[nev] (float array) --- Event depths
    nsta (int) --- No. of stations
    sta_lab[nsta] (object array) --- Station codes
    sta_lat[nsta] (float array) --- Station latitudes
    sta_lon[nsta] (float array) --- Station longitudes
    maxsep_ct
    maxsep_cc
    :::
    Returns:
    retlist (list) --- Return list of data arrays
    :::
    """

    """
    Declare counter arrays
    """
    ncc = int(0)
    nct = int(0)
    nctp = int(0)
    ncts = int(0)
    nccp = int(0)
    nccs = int(0)

    """
    Declare Full Data Arrays
    """
    nevp = 500000
    maxdat = 100
    maxdata = int(nevp*maxdat)          # Max. possible no. of data

    dt_sta1 = np.empty(maxdata,dtype='U8')      # Data station 1 codes
    dt_sta2 = np.empty(maxdata,dtype='U8')      # Data station 2 codes (stpair + doubpair)
    dt_dt = np.zeros(maxdata,dtype='float64')     # Data diff. times (double for doubpair)
    dt_qual = np.zeros(maxdata,dtype='float16')   # Data initial weights
    dt_offse = np.zeros(maxdata,dtype='float64')  # Data interevent offsets
    dt_offss = np.zeros(maxdata,dtype='float64')  # Data interstation offsets (stpair + doubpair)
    dt_c1 = np.zeros(maxdata,dtype='uint64')       # Data eventID 1
    dt_c2 = np.zeros(maxdata,dtype='uint64')       # Data eventID 2
    dt_idx = np.zeros(maxdata,dtype='uint16')      # Data type indexes
    dt_ista1 = np.zeros(maxdata,dtype='uint16')    # Data station 1 indexes
    dt_ista2 = np.zeros(maxdata,dtype='uint16')    # Data station 1 indexes
    dt_ic1 = np.zeros(maxdata,dtype='uint16')      # Data event 1 indexes
    dt_ic2 = np.zeros(maxdata,dtype='uint16')      # Data event 2 indexes

    """
    Sort event IDs
    """
    iicusp = np.argsort(ev_cusp)
    icusp = np.zeros(nev)
    icusp = ev_cusp[iicusp[0:nev]]

    if idata==1 or idata==3:
        """
        First read in CC values
        """
        [nccp,nccs,dt_sta1,dt_sta2,dt_dt,dt_qual,dt_offse,dt_offss,dt_c1,dt_c2,
         dt_idx] = readccs(log,reloctype,fn_cc,maxsep_cc,iphase,icusp,iicusp,ev_lat,
                           ev_lon,ev_dep,sta_lab,sta_lat,sta_lon,dt_sta1,dt_sta2,dt_dt,
                           dt_qual,dt_offse,dt_offss,dt_c1,dt_c2,dt_idx)
        ncc = nccp+nccs

    if idata==2 or idata==3:
        """
        Now read in catalog values and define catalog dts
        """
        [nctp,ncts,dt_sta1,dt_sta2,dt_dt,dt_qual,dt_offse,dt_offss,dt_c1,dt_c2,
         dt_idx] = readcts(log,reloctype,fn_ct,ncc,maxsep_ct,iphase,icusp,iicusp,ev_lat,
                           ev_lon,ev_dep,sta_lab,sta_lat,sta_lon,dt_sta1,dt_sta2,dt_dt,
                           dt_qual,dt_offse,dt_offss,dt_c1,dt_c2,dt_idx)
        nct = nctp+ncts

    """
    Trim data arrays to needed size
    """
    ndt = ncc+nct
    dt_sta1 = dt_sta1[:ndt]
    dt_sta2 = dt_sta2[:ndt]
    dt_dt = dt_dt[:ndt]
    dt_qual = dt_qual[:ndt]
    dt_offse = dt_offse[:ndt]
    dt_offss = dt_offss[:ndt]
    dt_c1 = dt_c1[:ndt]
    dt_c2 = dt_c2[:ndt]
    dt_idx = dt_idx[:ndt]
    if reloctype==1:
        dt_sta2 = []
        dt_offss = []
    if reloctype==2:
        dt_c2 = []
        dt_offse = []

    """
    Return List
    """
    retlist = [ndt,nccp,nccs,nctp,ncts,dt_sta1,dt_sta2,dt_dt,dt_qual,dt_offse,dt_offss,
               dt_c1,dt_c2,dt_idx]
    return retlist

def writeloc(fn_loc,nev,ev_cusp,ev_lat,ev_lon,ev_dep,ev_x,ev_y,ev_z,ev_herr,ev_zerr,
             ev_date,ev_time,ev_mag):
    """
    Write the hypoDD.loc file.
    :::
    Parameters:
    fn_loc (str) --- hypoDD.loc file location (from hypoDD.inp)
    nev (int) --- No. of events
    ev_cusp[nev] (int array) --- Event IDs
    ev_lat[nev] (float array) --- Event latitudes
    ev_lon[nev] (float array) --- Event longitudes
    ev_dep[nev] (float array) --- Event depths
    ev_x[nev] (float array) --- Event cartesian x coordinates
    ev_y[nev] (float array) --- Event cartesian y coordinates
    ev_z[nev] (float array) --- Event cartesian z coordinates
    ev_herr[nev] (float array) --- Event horizontal error
    ev_zerr[nev] (float array) --- Event vertical error
    ev_date[nev] (object array) --- Event dates
    ev_time[nev] (object array) --- Event times
    ev_mag[nev] (float array) --- Event magnitudes
    iclust (int) --- Current cluster ID
    :::
    Returns:
    NONE
    :::
    """
    loc = open(fn_loc,'a')  
    for i in range(0,nev):
        """
        For each event,
        write to file the event information
        """
        loc.write('%9i %10.8f %11.8f %9.6f %10.1f %10.1f %10.1f %8.1f %8.1f %8.1f %4i %2i %2i %2i %2i %5.2f %4.1f \n' %
                  (ev_cusp[i],ev_lat[i],ev_lon[i],ev_dep[i],ev_x[i],ev_y[i],ev_z[i],
                   ev_herr[i]*1000,ev_herr[i]*1000,ev_zerr[i]*1000,ev_date[i].year,
                   ev_date[i].month,ev_date[i].day,
                   ev_time[i].hour,ev_time[i].minute,
                   ev_time[i].second + ev_time[i].microsecond/1000000.,ev_mag[i]))
    loc.close()
    return None

def writereloc(fn_reloc,nev,src_cusp,src_lat,src_lon,src_dep,src_x,src_y,src_z,src_ex,
               src_ey,src_ez,ev_date,ev_time,ev_mag,icl,maxiter=False):
    """
    Write the hypoDD.loc file.
    :::
    Parameters:
    reloctype (int) --- Type of double-difference pairing
    fn_reloc (str) --- hypoDD.reloc file location (from hypoDD.inp)
    nev (int) --- No. of events
    src_cusp[nev] (int array) --- Event IDs
    src_lat[nev] (float array) --- Event latitudes
    src_lon[nev] (float array) --- Event longitudes
    src_dep[nev] (float array) --- Event depths
    src_x[nev] (float array) --- Event cartesian x coordinates
    src_y[nev] (float array) --- Event cartesian y coordinates
    src_z[nev] (float array) --- Event cartesian z coordinates
    src_ex[nev] (float array) --- Error x coordinates
    src_ey[nev] (float array) --- Error y coordinates
    src_ez[nev] (float array) --- Error z coordinates
    ev_date[nev] (object array) --- Event dates
    ev_time[nev] (object array) --- Event times
    ev_mag[nev] (float array) --- Event magnitudes
    icl (int) --- Cluster ID
    :::
    Returns:
    NONE
    :::
    """
    if maxiter:
        relocs = open(fn_reloc,'a')
    else:
        relocs = open(fn_reloc,'w')
    for i in range(0,nev):
        relocs.write('%9i %10.8f %11.8f %9.6f %10.1f %10.1f %10.1f %8.1f %8.1f %8.1f %4i %2i %2i %2i %2i %6.3f %4.1f %3i \n' %
                     (src_cusp[i],src_lat[i],src_lon[i],src_dep[i],src_x[i],src_y[i],src_z[i],src_ex[i],src_ey[i],src_ez[i],
                      ev_date[i].year,ev_date[i].month,ev_date[i].day,ev_time[i].hour,
                      ev_time[i].minute,ev_time[i].second + ev_time[i].microsecond/1000000.,ev_mag[i],int(icl)+1))
    relocs.close()
    return None

def terminaloutputs(reloctype,iteri,jiter,isolv,idata,nev,nevold,nct,nctold,ncc,nccold,
                    rms_ct,rms_ctold,rms_cc,rms_ccold,tmpr1,tmpr2,dxav,dyav,dzav,dtav,dsav,
                    xav,yav,zav,mbad,acond):
    """
    Write terminal outputs
    """
    #if jiter==0:
    #    str3 = '    '
    #else:
    n = iteri+1-jiter
    if n<1000:
        str3 = ' %3i' % n
    if n<100:
        str3 = '  %2i' % n
    if n<10:
        str3 = '   %1i' % n

    """
    Output for SVD and both CT and CC
    """
    if isolv==1 and idata==3: 
        if reloctype==1:
            if iteri==0:
                print('IT\tEV\tCT\tCC\t\t\tRMSCT\t\t\tRMSCC\t\tRMSST\t\tDX\t\tDY\t\tDZ\t\tDT\t\tOS\t\tAQ')
                print('\t%\t%\t%\t\tms\t%\t\tms\t%\t\tms\t\tm\t\tm\t\tm\t\tms\t\tm')
            print('%2i%s\t%3i\t%3i\t%3i\t\t%5i\t%5s\t\t%5i\t%5s\t\t%5i\t\t%4i\t\t%4i\t\t%4i\t\t%4i\t\t%4i\t\t%4i' %
                  (iteri+1,str3,np.rint(nev*100./nevold),np.rint(nct*100./nctold),
                   np.rint(ncc*100./nccold),
                   np.rint(rms_ct*1000),str(round((rms_ct-rms_ctold)*100./rms_ctold,1)).rjust(5),
                   np.rint(rms_cc*1000.),str(round((rms_cc-rms_ccold)*100./rms_ccold,1)).rjust(5),
                   np.rint(np.maximum(tmpr1,tmpr2)),np.rint(dxav),np.rint(dyav),
                   np.rint(dzav),np.rint(dtav),
                   np.rint(np.maximum(np.abs(xav),np.maximum(np.abs(yav),np.abs(zav)))),mbad))
        elif reloctype==2:
            if iteri==0:
                print(' IT    EV   CT   CC      RMSCT        RMSCC   RMSST   DX   DY   DZ   DS   OS   AQ')
                print('        %    %    %    ms      %    ms      %    ms    m    m    m   ms    m')
            print('%2i%s %3i  %3i  %3i %5i  %5s %5i  %5s %5i %4i %4i %4i %4i %4i %4i' %
                  (iteri+1,str3,np.rint(nev*100./nevold),np.rint(nct*100./nctold),
                   np.rint(ncc*100./nccold),
                   np.rint(rms_ct*1000),str(round((rms_ct-rms_ctold)*100./rms_ctold,1)).rjust(5),
                   np.rint(rms_cc*1000.),str(round((rms_cc-rms_ccold)*100./rms_ccold,1)).rjust(5),
                   np.rint(np.maximum(tmpr1,tmpr2)),np.rint(dxav),np.rint(dyav),
                   np.rint(dzav),np.rint(dsav),
                   np.rint(np.maximum(np.abs(xav),np.maximum(np.abs(yav),np.abs(zav)))),mbad))
        elif reloctype==3:
            if iteri==0:
                print(' IT    EV   CT   CC      RMSCT        RMSCC   RMSST   DX   DY   DZ   OS   AQ')
                print('        %    %    %    ms      %    ms      %    ms    m    m    m    m')
            print('%2i%s %3i  %3i  %3i %5i  %5s %5i  %5s %5i %4i %4i %4i %4i %4i' %
                  (iteri+1,str3,np.rint(nev*100./nevold),np.rint(nct*100./nctold),
                   np.rint(ncc*100./nccold),
                   np.rint(rms_ct*1000),str(round((rms_ct-rms_ctold)*100./rms_ctold,1)).rjust(5),
                   np.rint(rms_cc*1000.),str(round((rms_cc-rms_ccold)*100./rms_ccold,1)).rjust(5),
                   np.rint(np.maximum(tmpr1,tmpr2)),np.rint(dxav),np.rint(dyav),np.rint(dzav),
                   np.rint(np.maximum(np.abs(xav),np.maximum(np.abs(yav),np.abs(zav)))),mbad))

    """
    Output for SVD and CT only
    """
    if isolv==1 and idata==2: 
        if reloctype==1:
            if iteri==0:
                print(' IT    EV   CT      RMSCT   RMSST   DX   DY   DZ   DT   OS   AQ')
                print('        %    %    ms      %    ms    m    m    m   ms    m')
            print('%2i%s %3i  %3i %5i  %5s %5i %4i %4i %4i %4i %4i %4i' %
                  (iteri+1,str3,np.rint(nev*100./nevold),np.rint(nct*100./nctold),
                   np.rint(rms_ct*1000),str(round((rms_ct-rms_ctold)*100./rms_ctold,1)).rjust(5),
                   np.rint(np.maximum(tmpr1,tmpr2)),np.rint(dxav),np.rint(dyav),
                   np.rint(dzav),np.rint(dtav),
                   np.rint(np.maximum(np.abs(xav),np.maximum(np.abs(yav),np.abs(zav)))),mbad))
        elif reloctype==2:
            if iteri==0:
                print(' IT    EV   CT      RMSCT   RMSST   DX   DY   DZ   DS   OS   AQ')
                print('        %    %    ms      %    ms    m    m    m   ms    m')
            print('%2i%s %3i  %3i %5i  %5s %5i %4i %4i %4i %4i %4i %4i' %
                  (iteri+1,str3,np.rint(nev*100./nevold),np.rint(nct*100./nctold),
                   np.rint(rms_ct*1000),str(round((rms_ct-rms_ctold)*100./rms_ctold,1)).rjust(5),
                   np.rint(np.maximum(tmpr1,tmpr2)),np.rint(dxav),np.rint(dyav),
                   np.rint(dzav),np.rint(dsav),
                   np.rint(np.maximum(np.abs(xav),np.maximum(np.abs(yav),np.abs(zav)))),mbad))
        elif reloctype==3:
            if iteri==0:
                print(' IT    EV   CT      RMSCT   RMSST   DX   DY   DZ   OS   AQ')
                print('        %    %    ms      %    ms    m    m    m    m')
            print('%2i%s %3i  %3i %5i  %5s %5i %4i %4i %4i %4i %4i' %
                  (iteri+1,str3,np.rint(nev*100./nevold),np.rint(nct*100./nctold),
                   np.rint(rms_ct*1000),str(round((rms_ct-rms_ctold)*100./rms_ctold,1)).rjust(5),
                   np.rint(np.maximum(tmpr1,tmpr2)),np.rint(dxav),np.rint(dyav),np.rint(dzav),
                   np.rint(np.maximum(np.abs(xav),np.maximum(np.abs(yav),np.abs(zav)))),mbad))

    """
    Output for SVD and CC Only
    """
    if isolv==1 and idata==1:
        if reloctype==1:
            if iteri==0:
                print(' IT    EV   CC      RMSCC   RMSST   DX   DY   DZ   DT   OS   AQ')
                print('        %    %    ms      %    ms    m    m    m   ms    m')
            print('%2i%s %3i  %3i  %3i %5i  %5s %5i %4i %4i %4i %4i %4i %4i' %
                  (iteri+1,str3,np.rint(nev*100./nevold),
                   np.rint(ncc*100./nccold),
                   np.rint(rms_cc*1000.),str(round((rms_cc-rms_ccold)*100./rms_ccold,1)).rjust(5),
                   np.rint(np.maximum(tmpr1,tmpr2)),np.rint(dxav),np.rint(dyav),
                   np.rint(dzav),np.rint(dtav),
                   np.rint(np.maximum(np.abs(xav),np.maximum(np.abs(yav),np.abs(zav)))),mbad))
        elif reloctype==2:
            if iteri==0:
                print(' IT    EV   CC      RMSCC   RMSST   DX   DY   DZ   DS   OS   AQ')
                print('        %    %    ms      %    ms    m    m    m   ms    m')
            print('%2i%s %3i  %3i  %3i %5i  %5s %5i %4i %4i %4i %4i %4i %4i' %
                  (iteri+1,str3,np.rint(nev*100./nevold),
                   np.rint(ncc*100./nccold),
                   np.rint(rms_cc*1000.),str(round((rms_cc-rms_ccold)*100./rms_ccold,1)).rjust(5),
                   np.rint(np.maximum(tmpr1,tmpr2)),np.rint(dxav),np.rint(dyav),
                   np.rint(dzav),np.rint(dsav),
                   np.rint(np.maximum(np.abs(xav),np.maximum(np.abs(yav),np.abs(zav)))),mbad))
        elif reloctype==3:
            if iteri==0:
                print(' IT    EV   CC      RMSCC   RMSST   DX   DY   DZ   OS   AQ')
                print('        %    %    ms      %    ms    m    m    m    m')
            print('%2i%s %3i  %3i  %3i %5i  %5s %5i %4i %4i %4i %4i %4i' %
                  (iteri+1,str3,np.rint(nev*100./nevold),
                   np.rint(ncc*100./nccold),
                   np.rint(rms_cc*1000.),str(round((rms_cc-rms_ccold)*100./rms_ccold,1)).rjust(5),
                   np.rint(np.maximum(tmpr1,tmpr2)),np.rint(dxav),np.rint(dyav),np.rint(dzav),
                   np.rint(np.maximum(np.abs(xav),np.maximum(np.abs(yav),np.abs(zav)))),mbad))

    """
    Output for LSQR and both CC and CT
    """
    if isolv==2 and idata==3:
        if reloctype==1:
            if iteri==0:
                print('IT\tEV\tCT\tCC\tRMSCT\t\tRMSCC\t\tRMSST\tDX\tDY\tDZ\tDT\tOS\tAQ\tCND')
                print('\t%\t%\t%\tms\t%\t\tms\t%\tms\tm\tm\tm\tms\tm')
            print('%2i%s\t%3i\t%3i\t%3i\t%5i\t%5s\t%5i\t%5s\t%5i\t%4i\t%4i\t%4i\t%4i\t%4i\t%4i\t%5i' %
                  (iteri+1,str3,np.rint(nev*100./nevold),np.rint(nct*100./nctold),
                   np.rint(ncc*100./nccold),
                   np.rint(rms_ct*1000),str(round((rms_ct-rms_ctold)*100./rms_ctold,1)).rjust(5),
                   np.rint(rms_cc*1000.),str(round((rms_cc-rms_ccold)*100./rms_ccold,1)).rjust(5),
                   np.rint(np.maximum(tmpr1,tmpr2)),np.rint(dxav),np.rint(dyav),
                   np.rint(dzav),np.rint(dtav),
                   np.rint(np.maximum(np.abs(xav),np.maximum(np.abs(yav),np.abs(zav)))),mbad,np.rint(acond)))
        elif reloctype==2:
            if iteri==0:
                print(' IT    EV   CT   CC      RMSCT        RMSCC   RMSST   DX   DY   DZ   DS   OS   AQ   CND')
                print('        %    %    %    ms      %    ms      %    ms    m    m    m   ms    m')
            print('%2i%s %3i  %3i  %3i %5i  %5s %5i  %5s %5i %4i %4i %4i %4i %4i %4i %5i' %
                  (iteri+1,str3,np.rint(nev*100./nevold),np.rint(nct*100./nctold),
                   np.rint(ncc*100./nccold),
                   np.rint(rms_ct*1000),str(round((rms_ct-rms_ctold)*100./rms_ctold,1)).rjust(5),
                   np.rint(rms_cc*1000.),str(round((rms_cc-rms_ccold)*100./rms_ccold,1)).rjust(5),
                   np.rint(np.maximum(tmpr1,tmpr2)),np.rint(dxav),np.rint(dyav),
                   np.rint(dzav),np.rint(dsav),
                   np.rint(np.maximum(np.abs(xav),np.maximum(np.abs(yav),np.abs(zav)))),mbad,np.rint(acond)))
        elif reloctype==3:
            if iteri==0:
                print(' IT    EV   CT   CC      RMSCT        RMSCC   RMSST   DX   DY   DZ   OS   AQ   CND')
                print('        %    %    %    ms      %    ms      %    ms    m    m    m    m')
            print('%2i%s %3i  %3i  %3i %5i  %5s %5i  %5s %5i %4i %4i %4i %4i %4i %5i' %
                  (iteri+1,str3,np.rint(nev*100./nevold),np.rint(nct*100./nctold),
                   np.rint(ncc*100./nccold),
                   np.rint(rms_ct*1000),str(round((rms_ct-rms_ctold)*100./rms_ctold,1)).rjust(5),
                   np.rint(rms_cc*1000.),str(round((rms_cc-rms_ccold)*100./rms_ccold,1)).rjust(5),
                   np.rint(np.maximum(tmpr1,tmpr2)),np.rint(dxav),np.rint(dyav),np.rint(dzav),
                   np.rint(np.maximum(np.abs(xav),np.maximum(np.abs(yav),np.abs(zav)))),mbad,np.rint(acond)))

    """
    Output for LSQR and CT Only
    """
    if isolv==2 and idata==2:
        if reloctype==1:
            if iteri==0:
                print(' IT    EV   CT      RMSCT   RMSST   DX   DY   DZ   DT   OS   AQ   CND')
                print('        %    %    ms      %    ms    m    m    m   ms    m')
            print('%2i%s %3i  %3i %5i  %5s %5i %4i %4i %4i %4i %4i %4i %5i' %
                  (iteri+1,str3,np.rint(nev*100./nevold),np.rint(nct*100./nctold),
                   np.rint(rms_ct*1000),str(round((rms_ct-rms_ctold)*100./rms_ctold,1)).rjust(5),
                   np.rint(np.maximum(tmpr1,tmpr2)),np.rint(dxav),np.rint(dyav),
                   np.rint(dzav),np.rint(dtav),
                   np.rint(np.maximum(np.abs(xav),np.maximum(np.abs(yav),np.abs(zav)))),mbad,np.rint(acond)))
        elif reloctype==2:
            if iteri==0:
                print(' IT    EV   CT      RMSCT   RMSST   DX   DY   DZ   DS   OS   AQ   CND')
                print('        %    %    ms      %    ms    m    m    m   ms    m')
            print('%2i%s %3i  %3i %5i  %5s %5i %4i %4i %4i %4i %4i %4i %5i' %
                  (iteri+1,str3,np.rint(nev*100./nevold),np.rint(nct*100./nctold),
                   np.rint(rms_ct*1000),str(round((rms_ct-rms_ctold)*100./rms_ctold,1)).rjust(5),
                   np.rint(np.maximum(tmpr1,tmpr2)),np.rint(dxav),np.rint(dyav),
                   np.rint(dzav),np.rint(dsav),
                   np.rint(np.maximum(np.abs(xav),np.maximum(np.abs(yav),np.abs(zav)))),mbad,np.rint(acond)))
        elif reloctype==3:
            if iteri==0:
                print(' IT    EV   CT      RMSCT   RMSST   DX   DY   DZ   OS   AQ   CND')
                print('        %    %    ms      %    ms    m    m    m    m')
            print('%2i%s %3i  %3i %5i  %5s %5i %4i %4i %4i %4i %4i %5i' %
                  (iteri+1,str3,np.rint(nev*100./nevold),np.rint(nct*100./nctold),
                   np.rint(rms_ct*1000),str(round((rms_ct-rms_ctold)*100./rms_ctold,1)).rjust(5),
                   np.rint(np.maximum(tmpr1,tmpr2)),np.rint(dxav),np.rint(dyav),np.rint(dzav),
                   np.rint(np.maximum(np.abs(xav),np.maximum(np.abs(yav),np.abs(zav)))),mbad,np.rint(acond)))

    """
    Output for LSQR and CC only
    """
    if isolv==2 and idata==1: 
        if reloctype==1:
            if iteri==0:
                print(' IT    EV   CC      RMSCC   RMSST   DX   DY   DZ   DT   OS   AQ   CND')
                print('        %    %    ms      %    ms    m    m    m   ms    m')
            print('%2i%s %3i  %3i %5i  %5s %5i %4i %4i %4i %4i %4i %4i %5i' %
                  (iteri+1,str3,np.rint(nev*100./nevold),
                   np.rint(ncc*100./nccold),
                   np.rint(rms_cc*1000.),str(round((rms_cc-rms_ccold)*100./rms_ccold,1)).rjust(5),
                   np.rint(np.maximum(tmpr1,tmpr2)),np.rint(dxav),np.rint(dyav),
                   np.rint(dzav),np.rint(dtav),
                   np.rint(np.maximum(np.abs(xav),np.maximum(np.abs(yav),np.abs(zav)))),mbad,np.rint(acond)))
        elif reloctype==2:
            if iteri==0:
                print(' IT    EV   CC      RMSCC   RMSST   DX   DY   DZ   DS   OS   AQ   CND')
                print('        %    %    ms      %    ms    m    m    m   ms    m')
            print('%2i%s %3i  %3i %5i  %5s %5i %4i %4i %4i %4i %4i %4i %5i' %
                  (iteri+1,str3,np.rint(nev*100./nevold),
                   np.rint(ncc*100./nccold),
                   np.rint(rms_cc*1000.),str(round((rms_cc-rms_ccold)*100./rms_ccold,1)).rjust(5),
                   np.rint(np.maximum(tmpr1,tmpr2)),np.rint(dxav),np.rint(dyav),
                   np.rint(dzav),np.rint(dsav),
                   np.rint(np.maximum(np.abs(xav),np.maximum(np.abs(yav),np.abs(zav)))),mbad,np.rint(acond)))
        elif reloctype==3:
            if iteri==0:
                print(' IT    EV   CC      RMSCC   RMSST   DX   DY   DZ   OS   AQ   CND')
                print('        %    %    ms      %    ms    m    m    m    m')
            print('%2i%s %3i  %3i %5i  %5s %5i %4i %4i %4i %4i %4i %5i' %
                  (iteri+1,str3,np.rint(nev*100./nevold),
                   np.rint(ncc*100./nccold),
                   np.rint(rms_cc*1000.),str(round((rms_cc-rms_ccold)*100./rms_ccold,1)).rjust(5),
                   np.rint(np.maximum(tmpr1,tmpr2)),np.rint(dxav),np.rint(dyav),np.rint(dzav),
                   np.rint(np.maximum(np.abs(xav),np.maximum(np.abs(yav),np.abs(zav)))),mbad,np.rint(acond)))

    return None








