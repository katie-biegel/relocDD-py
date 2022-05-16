import numpy as np
"""
This script contains subfuntions from hypoDD that are used in the station-pair
double-difference method.
---
This script contains the functions:
    readcts_stpair      Read in the dt.ct catalog
    readccs_stpair      Read in the dt.cc catalog
    dtcal_st            Returns dt_cal array
    svdg_st             Returns the svd g kernel matrix
    svdunpack_st        Unpacks and returns svd inversion results
    lsqrprep_st         Preps arrays for lsqr inversion
---
"""

def readcts_stpair(log,fn_ct,ncc,maxsep_ct,iphase,icusp,iicusp,
                   ev_lat,ev_lon,ev_dep,sta_lab,sta_lat,sta_lon,
                   dt_sta1,dt_sta2,dt_dt,dt_qual,dt_offss,dt_c1,dt_idx):
    """
    Reads in the dt.ct file in the station-pair format.
    :::
    Parameters:
    log (file object) --- Log file
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
    dt_offss[ndt] (float array) --- Data interstation offsets
    dt_c1[ndt] (int array) --- Data event 1 IDs
    dt_idx[ndt] (int array) --- Data type indexes
    :::
    Returns:
    nctp (int) --- No. of P-phase cat data
    ncts (int) --- No. of S-phase cat data
    dt_sta1[ndt] (object array) --- Data station 1 codes
    dt_sta2[ndt] (object array) --- Data station 2 codes
    dt_dt[ndt] (float array) --- Data diff. time array
    dt_qual[ndt] (float array) --- Data weight array
    dt_offss[ndt] (float array) --- Data interstation offsets
    dt_c1[ndt] (int array) --- Data event 1 IDs
    dt_idx[ndt] (int array) --- Data type indexes
    :::
    """

    # Define PI for consistency
    PI = float(3.141593)

    """
    Counter Variables
    """
    nctp=int(0)
    ncts=int(0)
    i=ncc
    nsta=len(sta_lab)

    """
    Read in the full catalog file
    """
    ctfile = open(fn_ct,'r')
    cts = ctfile.readlines()
    for line in cts:
        """
        Read in Catalog Header Line
        ---
        Catalog header data follows the format:
        # EV1ID
        where
        EV1ID is the event 1 integer ID
        """

        if line[0:1]=='#':
            line = line.split()
            line = list(filter(None,line))

            """
            Read in Event IDs
            """
            ic1 = int(line[1])
            iskip = 0

            """
            Skip event pairs with events not in event list
            """
            try:
                k1 = int(tuple(np.argwhere(icusp==ic1))[0])
            except:
                iskip=1
                continue
        else:
            """
            Read in catalog line
            ---
            The catalog line format follows:
            STA1CODE STA2CODE TT1 TT2 WEIGHT PHA
            where
            STACODE is the string station code
            TT1 is the phase traveltime for event 1
            TT2 is the phase traveltime for event 2
            WEIGHT is the data quality (usually defined in the event pick catalog)
            PHA is a 'P' or 'S' character string
            """

            if iskip==1:    # If bad event pair, skip
                continue

            line = line.split()
            line = list(filter(None,line))

            """
            Read in catalog line
            """
            dt_sta1[i] = str(line[0])
            dt_sta2[i] = str(line[1])
            dt1 = float(line[2])
            dt2 = float(line[3])
            dt_qual[i] = float(line[4])
            pha = str(line[5])
            dt_c1[i] = ic1

            """
            Store time difference of the two traveltimes
            """
            dt_dt[i] = dt1-dt2

            """
            Skip data is station not saved
            """
            if dt_sta1[i] in sta_lab and dt_sta2[i] in sta_lab:
                """
                Check station-pair offset
                """
                i1 = np.argwhere(sta_lab==dt_sta1[i])[0][0]
                slat1=sta_lat[i1]
                slon1=sta_lon[i1]
                i2 = np.argwhere(sta_lab==dt_sta2[i])[0][0]
                slat2=sta_lat[i2]
                slon2=sta_lon[i2]

                slat=slat1-slat2
                slon=slon1-slon2
                offs=np.sqrt((slat*111)**2 + 
                             (slon*np.cos(slat1*PI/180)*111)**2)
                #if maxsep_ct>0 and offs>maxsep_ct:
                #    continue

                """
                Store data correctly based on phase
                """
                if pha=='P' and iphase!=2:
                    dt_idx[i]=3
                    nctp+=1
                elif pha=='S' and iphase!=1:
                    dt_idx[i]=4
                    ncts+=1
                else:
                    continue
            
                dt_offss[i] = offs
                i+=1

    ctfile.close()

    return nctp,ncts,dt_sta1,dt_sta2,dt_dt,dt_qual,dt_offss,dt_c1,dt_idx

#@profile(stream=open('mem_logs/readccs_stpair.mem','w+'))
def readccs_stpair(log,fn_cc,maxsep_cc,iphase,icusp,iicusp,ev_lat,ev_lon,ev_dep,
                   sta_lab,sta_lat,sta_lon,
                   dt_sta1,dt_sta2,dt_dt,dt_qual,dt_offss,dt_c1,dt_idx):
    """
    Reads in the dt.cc file in the event-pair format.
    :::
    Parameters:
    log (file object) --- Log file
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
    dt_offss[ndt] (float array) --- Data interevent offsets
    dt_c1[ndt] (int array) --- Data event 1 IDs
    dt_idx[ndt] (int array) --- Data type indexes
    :::
    Returns:
    nccp (int) --- No. of P-phase cc data
    nccs (int) --- No. of S-phase cc data
    dt_sta1[ndt] (object array) --- Data station 1 codes
    dt_sta2[ndt] (object array) --- Data station 2 codes
    dt_dt[ndt] (float array) --- Data diff. time array
    dt_qual[ndt] (float array) --- Data weight array
    dt_offss[ndt] (float array) --- Data interevent offsets
    dt_c1[ndt] (int array) --- Data event 1 IDs
    dt_idx[ndt] (int array) --- Data type indexes
    :::
    """

    # Define PI for consistency
    PI = float(3.141593)
    """
    Counter Variables
    """
    nccp = int(0)
    nccs = int(0)
    i=int(0)
    iiotc=int(0)
    nsta=len(sta_lab)

    """
    Read CC File Loop
    """
    # Read in the file
    ccfile = open(fn_cc,'r')
    ccs = ccfile.readlines()
    for line in ccs:
        """
        Read in the Event Header Line
        ---
        Event Headers in the CC Format MUST Contain
        # EV1ID OTC
        where 
        EV1ID is the event 1 integer ID
        OTC is the origin time correction
                    **seems defined in early examples from
                        Waldhauser and Ellsworth but not used
                        in many examples.  Set to 0.0 if unsure.
        All values are required.
        """
        if line[0:1]=='#':
            line = line.split()
            line = list(filter(None,line))

            """
            Read in Values
            """
            ic1 = int(line[1])
            try:
                otc = float(line[2])
            except:
                """
                Skip events with no otc (origin time correction)
                ---
                OTC can be defined as 0.0 (should be in most cases)
                """
                log.write('No otc for %i. Event Skipped. \n' % (ic1))
                print('No otc for %i. Event skipped.' % (ic1))
                iiotc+=1
                iskip=1
                continue
            iskip = int(0)

            """
            Skip event if not in event list
            """
            try:
                k1 = int(tuple(np.argwhere(icusp==ic1))[0])
            except:
                iskip=1
                continue

        else:
            """
            Read in cc lines
            ---
            The data format for cc lines:
            STACODE1 STACODE2 DT_CC WEIGHT PHASE
            where
            STACODE is the station string code
            DT_CC is the diff. time defined by the CC
            WEIGHT is the initial data weight (I usually define this as the CC coefficent)
            PHASE is a 'P' or 'S' character string
            """
            if iskip==1:    # Skip if bad event pair
                continue

            line = line.split()
            line = list(filter(None,line))

            """
            Read in cc lines
            """
            dt_sta1[i] = str(line[0])
            dt_sta2[i] = str(line[1])
            dt_dt[i] = float(line[2])
            dt_qual[i] = float(line[3])
            pha = str(line[4])
            dt_c1[i] = ic1
            dt_dt[i] = dt_dt[i]-otc # OTC term applied directly to the dt

            """
            Skip if bad station
            """
            if not dt_sta1[i] in sta_lab:
                continue
            if not dt_sta2[i] in sta_lab:
                continue

            """
            Check station-pair offset
            """
            for j in range(0,nsta):
                if dt_sta1[i]==sta_lab[j]:
                    slat1=sta_lat[j]
                    slon1=sta_lon[j]
                elif dt_sta2[i]==sta_lab[j]:
                    slat2=sta_lat[j]
                    slon2=sta_lon[j]
            slat=slat1-slat2
            slon=slon1-slon2
            offs=np.sqrt((slat*111)**2 + 
                         (slon*np.cos(slat1*PI/180)*111)**2)
            if maxsep_cc>0 and offs>maxsep_cc:
                continue

            """
            Store data correctly based on phase
            """
            if pha=='P' and iphase!=2:
                dt_idx[i]=1
                nccp+=1
            elif pha=='S' and iphase!=1:
                dt_idx[i]=2
                nccs+=1
            else:
                raise Exception('>>> Phase identifier format error')
            
            dt_offss[i] = offs
            i += 1
    ccfile.close()

    return nccp,nccs,dt_sta1,dt_sta2,dt_dt,dt_qual,dt_offss,dt_c1,dt_idx,iiotc

#@profile(stream=open('mem_logs/dtcal_stpair.mem','w+'))
def dtcal_st(ndt,dt_cal,dt_dt,dt_idx,dt_ista1,dt_ista2,dt_ic1,src_t,tmp_ttp,tmp_tts):
    """
    This function returns the dt_cal array (calculated diff. time)
    :::
    Parameters:
    ndt (int) --- No. of data
    dt_dt[ndt] (float array) --- Data diff. time
    dt_idx[ndt] (int array) --- Data type indexes
    dt_ista1[ndt] (int array) --- Data station indexes
    dt_ista2[ndt] (int array) --- Data station indexes
    dt_ic1[ndt] (int array) --- Data event 1 indexes
    src_t[nev] (float array) --- Source times
    tmp_ttp[nsta,nev] (float array) --- P raytracing times
    tmp_tts[nsta,nev] (float array) --- S raytracing times
    :::
    Returns:
    dt_cal[ndt] --- Calcualted diff. times from raytracing
    :::
    """

    tt1 = 0.
    tt2 = 0.
    for i in range(ndt):
        """
        Pull correct traveltimes for pair
        """
        if dt_idx[i]==1 or dt_idx[i]==3: # P-phase
            tt1 = tmp_ttp[dt_ista1[i],dt_ic1[i]] - src_t[dt_ic1[i]]/1000.
            tt2 = tmp_ttp[dt_ista2[i],dt_ic1[i]] - src_t[dt_ic1[i]]/1000.
        elif dt_idx[i]==2 or dt_idx[i]==4: # S phase
            tt1 = tmp_tts[dt_ista1[i],dt_ic1[i]] - src_t[dt_ic1[i]]/1000.
            tt2 = tmp_tts[dt_ista2[i],dt_ic1[i]] - src_t[dt_ic1[i]]/1000.
        if tt1==0. or tt2==0.:
            raise Exception('Fatal Error (theor tt)')

        """
        Calculated diff. time
        """
        dt_cal[i] = tt1-tt2

    return dt_cal,tt1,tt2


def svdg_st(ndt,nev,nsrc,dt_ic1,dt_ista1,dt_ista2,tmp_xp,tmp_yp,tmp_zp,dt_idx,mod_ratio):
    """
    This function calculates and returns the g matrix for SVD inversion
    :::
    Parameters:
    ndt (int) ---- No. of data
    nev (int) ---- No. of events
    nsrc (int) ---- No. of sources
    dt_ic1[ndt] (int array) ---- Array holding ev1 indexes
    dt_ista1[ndt] (int array) ---- Array holding station indexes
    dt_ista2[ndt] (int array) ---- Array holding station indexes
    tmp_xp[nsta,nev] (float array) --- X partial derivatives from ray tracing
    tmp_yp[nsta,nev] (float array) --- Y partial derivatives from ray tracing
    tmp_zp[nsta,nev] (float array) --- Z partial derivatives from ray tracing
    dt_idx[ndt] (float array) ---- Array holding data type indexes
    mod_ratio (float) ---- VPVS ratio
    :::
    Returns:
    g[ndt+4,4*nev] (float array) --- SVD inversion g kernel matrix
    :::
    """

    # Declare G (extra 4 rows for symmetry)
    g = np.zeros((ndt+4,4*nev))

    # If all events start from same source:
    if nsrc==1:
        for i in range(0,ndt):
            k3 = 4*dt_ic1[i]

            g[i,k3] = (tmp_xp[dt_ista1[i],0]-tmp_xp[dt_ista2[i],0])
            g[i,k3+1] = (tmp_yp[dt_ista1[i],0]-tmp_yp[dt_ista2[i],0])
            g[i,k3+2] = (tmp_zp[dt_ista1[i],0]-tmp_zp[dt_ista2[i],0])
            g[i,k3+3] = 1.
    else:        
        k1 = dt_ic1[:ndt]
        k3 = 4*k1
        for i in range(0,ndt):
            g[i,k3[i]] = (tmp_xp[dt_ista1[i],k1[i]]-tmp_xp[dt_ista2[i],k1[i]])
            g[i,k3[i]+1] = (tmp_yp[dt_ista1[i],k1[i]]-tmp_yp[dt_ista2[i],k1[i]])
            g[i,k3[i]+2] = (tmp_zp[dt_ista1[i],k1[i]]-tmp_zp[dt_ista2[i],k1[i]])
            g[i,k3[i]+3] = 1.

    # Multiply by the VPVS ratio where data describes S waves
    g[0:ndt,:] = np.where(np.logical_or(dt_idx[0:ndt].reshape((ndt,1))==2,dt_idx[0:ndt].reshape((ndt,1))==4),g[0:ndt,:]*mod_ratio,g[0:ndt,:])

    return g

#@profile(stream=open('mem_logs/svdunpack_st.mem','w+'))
def svdunpack_st(log,nev,iteri,cvm,factor,norm,x,ev_cusp,exav,eyav,ezav,esav,
                 dxav,dyav,dzav,dsav):
    """
    Unpack the svd results
    :::
    Parameters:
    log (file object) --- Log file
    nev (int) --- No. of events
    iteri (int) --- Iteration number
    cvm[4*nev,4*nev] (float array) --- Covariance matrix
    factor (float) --- 95% confidence level
    norm[4*nev] (float array) --- G matrix column norms
    x[4*nev] (float array) --- Model matrix
    ev_cusp[nev] (int array) --- Event IDs
    exav (float) ---- Avg. x error
    eyav (float) ---- Avg. y error
    ezav (float) ---- Avg. z error
    esav (float) ---- Avg. t error
    dxav (float) ---- Change in x value
    dyav (float) ---- Change in y value
    dzav (float) ---- Change in z value
    dsav (float) ---- Change in t value
    :::
    Returns:
    retlist (list) ---- List of returned variables
    :::
    """

    """
    Error matrix
    """
    se = np.zeros(4*nev)
    for i in range(0,4*nev):
        se[i] = np.sqrt(cvm[i,i])*factor
    # Rescale Errors
    se = se/norm

    # Store solution and errors
    src_dx = np.zeros(nev)
    src_dy = np.zeros(nev)
    src_dz = np.zeros(nev)
    src_ds = np.zeros(nev)
    src_ex = np.zeros(nev)
    src_ey = np.zeros(nev)
    src_ez = np.zeros(nev)
    src_es = np.zeros(nev)
    src_dx = -x[0:4*nev-3:4]
    src_dy = -x[1:4*nev-2:4]
    src_dz = -x[2:4*nev-1:4]
    src_ds = -x[3:4*nev:4]
    src_ex = se[0:4*nev-3:4]
    src_ey = se[1:4*nev-2:4]
    src_ez = se[2:4*nev-1:4]
    src_es = se[3:4*nev:4]
    src_cusp = np.copy(ev_cusp)

    # Output statistics
    exavold = exav
    eyavold = eyav
    ezavold = ezav 
    esavold = esav 
    dxavold = dxav 
    dyavold = dyav 
    dzavold = dzav 
    dsavold = dsav 
    exav = np.sum(src_ex)/nev
    eyav = np.sum(src_ey)/nev
    ezav = np.sum(src_ez)/nev
    esav = np.sum(src_es)/nev
    dxav = np.sum(np.abs(src_dx))/nev
    dyav = np.sum(np.abs(src_dy))/nev
    dzav = np.sum(np.abs(src_dz))/nev
    dsav = np.sum(np.abs(src_ds))/nev
    if iteri==1:
        exavold = exav
        eyavold = eyav
        ezavold = ezav
        esavold = esav
        dxavold = dxav
        dyavold = dyav
        dzavold = dzav
        dsavold = dsav

    log.write('Location summary: \n')
    log.write('mean 2sig-error (x,y,z,st) [m,ms]: %7.1f %7.1f %7.1f %7.1f \n' % (exav,eyav,ezav,esav))
    log.write('  ( %7.1f %7.1f %7.1f %7.1f )  \n' % (exav-exavold,eyav-eyavold,ezav-ezavold,esav-esavold))
    log.write('mean shift (x,y,z,st) [m,ms] (DX,DY,DZ,DT): %7.1f %7.1f %7.1f %7.1f \n' % (dxav,dyav,dzav,dsav))
    log.write('  ( %7.1f %7.1f %7.1f %7.1f )  \n' % (dxav-dxavold,dyav-dyavold,dzav-dzavold,dsav-dsavold))

    retlist = [src_dx,src_dy,src_dz,src_ds,src_ex,src_ey,src_ez,src_es,src_cusp,
               exavold,eyavold,ezavold,esavold,dxavold,dyavold,dzavold,dsavold,
               exav,eyav,ezav,esav,dxav,dyav,dzav,dsav]
    return retlist

#@profile(stream=open('mem_logs/lsqrprep_st.mem','w+'))
def lsqrprep_st(log,ndt,nsrc,dt_idx,tmp_xp,tmp_yp,tmp_zp,dt_ista1,dt_ista2,dt_ic1,mod_ratio,wt):
    """
    Prepare kernal matrices for LSQR inversion
    :::
    Parameters:
    log (file obj) ---- Log file
    ndt (int) ---- No. of data
    nsrc (int) ---- No. of events
    dt_idx[ndt] (float array) ---- Array holding data type indexes
    tmp_xp[NSTA,NEV] (float array) --- X partial derivatives from ray tracing
    tmp_yp[NSTA,NEV] (float array) --- Y partial derivatives from ray tracing
    tmp_zp[NSTA,NEV] (float array) --- Z partial derivatives from ray tracing
    dt_ista1[ndt] (int array) ---- Array holding station indexes
    dt_ista2[ndt] (int array) ---- Array holding station indexes
    dt_ic1[ndt] (int array) ---- Array holding ev1 indexes
    mod_ratio (float) ---- VPVS ratio
    wt[ndt] (float array) ---- Weights
    :::
    Returns:
    nar (int) --- Number of g matrix values (4 per data)
    nndt (int) --- Number of data 
    row_i[nar] (int array) --- G matrix row indexes
    col_i[nar] (int array) --- G matrix column indexes
    rw[nar] (int array) --- G matrix values
    :::
    """
    log.write('~ Setting up G matrix: \n')

    nar = 4*ndt
    nndt = ndt
    row_i = np.empty((nar),dtype='uint32')   # Row indexes for sparse matrix
    col_i = np.empty((nar),dtype='uint32')   # Column indexes for sparse matrix
    rw = np.zeros(nar,dtype='float64')       # Values for sparse G matrix
    
    # Set up row and column indexes and data input for sparse a matrix
    row_i[0:ndt] = np.arange(ndt)
    row_i[ndt:2*ndt] = np.arange(ndt)
    row_i[2*ndt:3*ndt] = np.arange(ndt)
    row_i[3*ndt:4*ndt] = np.arange(ndt)

    if nsrc==1:
        k1 = 0
        #rw[0:ndt]       = np.where((dt_idx[0:ndt]==2)|(dt_idx[0:ndt]==4),(tmp_xp[dt_ista1[0:ndt],k1]-tmp_xp[dt_ista2[0:ndt],k1])*wt[0:ndt]*mod_ratio,(tmp_xp[dt_ista1[0:ndt],k1]-tmp_xp[dt_ista2[0:ndt],k1])*wt[0:ndt])
        #rw[ndt:2*ndt]   = np.where((dt_idx[0:ndt]==2)|(dt_idx[0:ndt]==4),(tmp_yp[dt_ista1[0:ndt],k1]-tmp_yp[dt_ista2[0:ndt],k1])*wt[0:ndt]*mod_ratio,(tmp_yp[dt_ista1[0:ndt],k1]-tmp_yp[dt_ista2[0:ndt],k1])*wt[0:ndt])
        #rw[2*ndt:3*ndt] = np.where((dt_idx[0:ndt]==2)|(dt_idx[0:ndt]==4),(tmp_zp[dt_ista1[0:ndt],k1]-tmp_zp[dt_ista2[0:ndt],k1])*wt[0:ndt]*mod_ratio,(tmp_zp[dt_ista1[0:ndt],k1]-tmp_zp[dt_ista2[0:ndt],k1])*wt[0:ndt])
        #rw[3*ndt:4*ndt] = wt[0:ndt]
        rw[0:ndt]       = np.where((dt_idx[0:ndt]==2)|(dt_idx[0:ndt]==4),(tmp_xp[dt_ista1[0:ndt],k1])*wt[0:ndt]*mod_ratio,(tmp_xp[dt_ista1[0:ndt],k1])*wt[0:ndt])
        rw[ndt:2*ndt]   = np.where((dt_idx[0:ndt]==2)|(dt_idx[0:ndt]==4),(tmp_yp[dt_ista1[0:ndt],k1])*wt[0:ndt]*mod_ratio,(tmp_yp[dt_ista1[0:ndt],k1])*wt[0:ndt])
        rw[2*ndt:3*ndt] = np.where((dt_idx[0:ndt]==2)|(dt_idx[0:ndt]==4),(tmp_zp[dt_ista1[0:ndt],k1])*wt[0:ndt]*mod_ratio,(tmp_zp[dt_ista1[0:ndt],k1])*wt[0:ndt])
        rw[3*ndt:4*ndt] = wt[0:ndt]
        rw[0:ndt]      += np.where((dt_idx[0:ndt]==2)|(dt_idx[0:ndt]==4),(-tmp_xp[dt_ista2[0:ndt],k1])*wt[0:ndt]*mod_ratio,(-tmp_xp[dt_ista2[0:ndt],k1])*wt[0:ndt])
        rw[ndt:2*ndt]  += np.where((dt_idx[0:ndt]==2)|(dt_idx[0:ndt]==4),(-tmp_yp[dt_ista2[0:ndt],k1])*wt[0:ndt]*mod_ratio,(-tmp_yp[dt_ista2[0:ndt],k1])*wt[0:ndt])
        rw[2*ndt:3*ndt]+= np.where((dt_idx[0:ndt]==2)|(dt_idx[0:ndt]==4),(-tmp_zp[dt_ista2[0:ndt],k1])*wt[0:ndt]*mod_ratio,(-tmp_zp[dt_ista2[0:ndt],k1])*wt[0:ndt])
        rw[3*ndt:4*ndt]+= -wt[0:ndt]


        # Set up column indexes with non-zero elements
        col_i[0:ndt] = 4*dt_ic1[0:ndt]
        col_i[ndt:2*ndt] = 4*dt_ic1[0:ndt]+1
        col_i[2*ndt:3*ndt] = 4*dt_ic1[0:ndt]+2
        col_i[3*ndt:4*ndt] = 4*dt_ic1[0:ndt]+3

    else:
        k1 = dt_ic1[0:ndt].astype(int)

        rw[0:ndt]       = np.where((dt_idx[0:ndt]==2)|(dt_idx[0:ndt]==4),(tmp_xp[dt_ista1[0:ndt],k1[0:ndt]]-tmp_xp[dt_ista2[0:ndt],k1[0:ndt]])*wt[0:ndt]*mod_ratio,(tmp_xp[dt_ista1[0:ndt],k1[0:ndt]]-tmp_xp[dt_ista2[0:ndt],k1[0:ndt]])*wt[0:ndt])
        rw[ndt:2*ndt]   = np.where((dt_idx[0:ndt]==2)|(dt_idx[0:ndt]==4),(tmp_yp[dt_ista1[0:ndt],k1[0:ndt]]-tmp_yp[dt_ista2[0:ndt],k1[0:ndt]])*wt[0:ndt]*mod_ratio,(tmp_yp[dt_ista1[0:ndt],k1[0:ndt]]-tmp_yp[dt_ista2[0:ndt],k1[0:ndt]])*wt[0:ndt])
        rw[2*ndt:3*ndt] = np.where((dt_idx[0:ndt]==2)|(dt_idx[0:ndt]==4),(tmp_zp[dt_ista1[0:ndt],k1[0:ndt]]-tmp_zp[dt_ista2[0:ndt],k1[0:ndt]])*wt[0:ndt]*mod_ratio,(tmp_zp[dt_ista1[0:ndt],k1[0:ndt]]-tmp_zp[dt_ista2[0:ndt],k1[0:ndt]])*wt[0:ndt])
        rw[3*ndt:4*ndt] = wt[0:ndt]
        #rw[0:ndt]       = np.where((dt_idx[0:ndt]==2)|(dt_idx[0:ndt]==4),(tmp_xp[dt_ista1[0:ndt],k1[0:ndt]])*wt[0:ndt]*mod_ratio,(tmp_xp[dt_ista1[0:ndt],k1[0:ndt]])*wt[0:ndt])
        #rw[ndt:2*ndt]   = np.where((dt_idx[0:ndt]==2)|(dt_idx[0:ndt]==4),(tmp_yp[dt_ista1[0:ndt],k1[0:ndt]])*wt[0:ndt]*mod_ratio,(tmp_yp[dt_ista1[0:ndt],k1[0:ndt]])*wt[0:ndt])
        #rw[2*ndt:3*ndt] = np.where((dt_idx[0:ndt]==2)|(dt_idx[0:ndt]==4),(tmp_zp[dt_ista1[0:ndt],k1[0:ndt]])*wt[0:ndt]*mod_ratio,(tmp_zp[dt_ista1[0:ndt],k1[0:ndt]])*wt[0:ndt])
        #rw[3*ndt:4*ndt] = wt[0:ndt]
        #rw[0:ndt]      += np.where((dt_idx[0:ndt]==2)|(dt_idx[0:ndt]==4),(-tmp_xp[dt_ista2[0:ndt],k1[0:ndt]])*wt[0:ndt]*mod_ratio,(-tmp_xp[dt_ista2[0:ndt],k1[0:ndt]])*wt[0:ndt])
        #rw[ndt:2*ndt]  += np.where((dt_idx[0:ndt]==2)|(dt_idx[0:ndt]==4),(-tmp_yp[dt_ista2[0:ndt],k1[0:ndt]])*wt[0:ndt]*mod_ratio,(-tmp_yp[dt_ista2[0:ndt],k1[0:ndt]])*wt[0:ndt])
        #rw[2*ndt:3*ndt]+= np.where((dt_idx[0:ndt]==2)|(dt_idx[0:ndt]==4),(-tmp_zp[dt_ista2[0:ndt],k1[0:ndt]])*wt[0:ndt]*mod_ratio,(-tmp_zp[dt_ista2[0:ndt],k1[0:ndt]])*wt[0:ndt])
        #rw[3*ndt:4*ndt]+= -wt[0:ndt]

        # Set up column indexes with non-zero elements
        col_i[0:ndt] = 4*k1
        col_i[ndt:2*ndt] = 4*k1+1
        col_i[2*ndt:3*ndt] = 4*k1 + 2
        col_i[3*ndt:4*ndt] = 4*k1 + 3

    return nar,nndt,row_i,col_i,rw

#@profile(stream=open('mem_logs/lsqrunpack_st.mem','w+'))
def lsqrunpack_st(log,nev,iteri,x,se,resvar1,factor,ev_cusp,
                  exav,eyav,ezav,esav,dxav,dyav,dzav,dsav):
    """
    Unpack the lsqr inversion results.
    :::
    Parameters:
    log (file object) --- Log file
    nev (int) --- No. of events
    iteri (int) --- Iteration number
    x[4*nev] (float array) --- Model parameters
    se[4*nev] (float array) --- Model error parameters
    resvar1 (float) --- Resstat avg. residual
    factor (float) --- 95% confidence interval
    ev_cusp[nev] (int array) --- Event IDs
    exav (float) ---- Avg. x error
    eyav (float) ---- Avg. y error
    ezav (float) ---- Avg. z error
    esav (float) ---- Avg. st. corr. term error
    dxav (float) ---- Change in x value
    dyav (float) ---- Change in y value
    dzav (float) ---- Change in z value
    dsav (float) ---- Change in station correction value
    :::
    Returns:
    retlist (list) --- List of variables to return
    :::
    """

    # Store solution and errors
    src_dx = np.zeros(nev)
    src_dy = np.zeros(nev)
    src_dz = np.zeros(nev)
    src_ds = np.zeros(nev)

    src_ex = np.zeros(nev)
    src_ey = np.zeros(nev)
    src_ez = np.zeros(nev)
    src_es = np.ones(nev)

    src_dx = -x[0:4*nev-3:4]
    src_dy = -x[1:4*nev-2:4]
    src_dz = -x[2:4*nev-1:4]
    src_ds = -x[3:4*nev:4]

    src_ex = np.sqrt(se[0:4*nev-3:4])*np.sqrt(resvar1)*factor
    src_ey = np.sqrt(se[1:4*nev-2:4])*np.sqrt(resvar1)*factor
    src_ez = np.sqrt(se[2:4*nev-1:4])*np.sqrt(resvar1)*factor
    src_es = np.sqrt(se[3:4*nev:4])*np.sqrt(resvar1)*factor
    src_cusp = np.copy(ev_cusp)

    #Get average errors and vector changes
    exavold = exav
    eyavold = eyav
    ezavold = ezav
    esavold = esav

    dxavold = dxav
    dyavold = dyav
    dzavold = dzav
    dsavold = dsav

    exav = np.sum(src_ex)/nev
    eyav = np.sum(src_ey)/nev
    ezav = np.sum(src_ez)/nev
    esav = np.sum(src_es)/nev

    dxav = np.sum(np.abs(src_dx))/nev
    dyav = np.sum(np.abs(src_dy))/nev
    dzav = np.sum(np.abs(src_dz))/nev
    dsav = np.sum(np.abs(src_ds))/nev

    if iter==1:
        exavold = exav
        eyavold = eyav
        ezavold = ezav
        esavold = esav

        dxavold = dxav
        dyavold = dyav
        dzavold = dzav
        dsavold = dsav

    # Output location statistics
    log.write('Location summary: \n')
    log.write(' mean 2sig-error (x,y,z,t) [m,ms]: %7.1f, %7.1f, %7.1f, %7.1f, ( %7.1f, %7.1f, %7.1f, %7.1f), \n' %
              (exav, eyav, ezav, esav, exav-exavold, eyav-eyavold, ezav-ezavold, esav-esavold))
    log.write(' mean shift (x,y,z,t) [m,ms] (DX,DY,DZ,DT): %7.1f, %7.1f, %7.1f, %7.1f, ( %7.1f, %7.1f, %7.1f, %7.1f), \n' %
              (dxav, dyav, dzav, dsav, dxav-dxavold, dyav-dyavold, dzav-dzavold, dsav-dsavold))

    retlist = [src_dx,src_dy,src_dz,src_ds,src_ex,src_ey,src_ez,src_es,src_cusp,
               exavold,eyavold,ezavold,esavold,dxavold,dyavold,dzavold,dsavold,
               exav,eyav,ezav,esav,dxav,dyav,dzav,dsav]
    return retlist


