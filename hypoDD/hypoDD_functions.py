import numpy as np
import datetime
import os
# Other HypoDD functions
from utility.universal.geodetics import sdc2
from hypoDD.hypoDD_files import readdata,readevent
from utility.universal.misc import readstat
from methodTypes.eventPair.hypoDD_ev import dtcal_ev
from methodTypes.stationPair.hypoDD_stat import dtcal_st
from methodTypes.doublePair.hypoDD_doub import dtcal_doub

"""
This script contains assorted supplementary functions used in the
preparation or running of the hypoDD command.
---
This script includes the functions:
    cleanarrays:        Function cleaning the event and station arrays
                        and updating the data indexing arrays
    cluster:            Clusters based on nearest neighbor
    trialsrc:           Sets up the initial trial locations
    dtres:              Calculate the double-difference values
    weighting:          Reweights data
    skip:               Cleans and skips data
    resstat:            Calculates residual statistics
    dataprep:           Read in and prep data for hypoDD
    sigcoherency:       Recalculate inter-pair offsets and signal/pick coherency
    statres:            Station residuals
"""

def cleanarrays(log,reloctype,nev,ev_cusp,ev_date,ev_time,ev_lat,ev_lon,ev_dep,ev_mag,
                ev_herr,ev_zerr,ev_res,nsta,sta_lab,sta_lat,sta_lon,ndt,dt_c1,dt_c2,
                dt_sta1,dt_sta2,dt_ic1,dt_ic2,dt_ista1,dt_ista2,amcusp=[]):
    """
    This functions cleans the event, station, and indexing
    data arrays.
    :::
    Parameters:
    log (file object) --- Log file
    reloctype (int) --- Double-difference type
    nev (int) --- No. of events
    ev_cusp[nev] (int array) --- Event IDs
    ev_date[nev] (object array) --- Event dates (datetime.date objects)
    ev_time[nev] (object array) --- Event times (datetime.time objects)
    ev_lat[nev] (float array) --- Event latitudes
    ev_lon[nev] (float array) --- Event longitudes
    ev_dep[nev] (float array) --- Event depths
    ev_mag[nev] (float array) --- Event mag
    ev_herr[nev] (float array) --- Event horizontal errors
    ev_zerr[nev] (float array) --- Event vertical errors
    ev_res[nev] (float array) --- Event residuals
    nsta (int) --- No. of stationss
    sta_lab[nsta] (object array) --- Station codes
    sta_lat[nsta] (float array) --- Station latitudes
    sta_lon[nsta] (float array) --- Station longitudes
    ndt (int) --- No. of data
    dt_c1[ndt] (int array) --- Data event 1 IDs
    dt_c2[ndt] (int array) --- Data event 2 IDs
    dt_sta1[ndt] (object array) --- Data station codes
    dt_sta2[ndt] (object array) --- Data station codes
    dt_ic1[ndt] (int array) --- Data event 1 indexes
    dt_ic2[ndt] (int array) --- Data event 2 indexes
    dt_ista1[ndt] (int array) --- Data station indexes
    dt_ista2[ndt] (int array) --- Data station indexes
    :::
    Returns:
    retlist (list) --- Updated arrays
    :::
    """
    log.write('\n\nCleanarrays ...\n')

    """
    Clean events
    ---
    This loop cleans out any events that do not have
    associated data.
    """
    #dt_ic1 = np.sort(dt_c1[0:ndt])
    #if reloctype==1 or reloctype==3:
    #    dt_ic2 = np.sort(dt_c2[0:ndt])
    k=int(0)
    ################ KB Note: Entirely Possible to remove loop revisit
    #print('Events before: ',nev)
    for i in range(0,nev):
        if reloctype==2:
            if len(amcusp)>0 and ev_cusp[i] in amcusp:
                #print('Airquake: ',ev_cusp[i])
                continue
            if ev_cusp[i] in dt_c1:
                ev_date[k] = ev_date[i]
                ev_time[k] = ev_time[i]
                ev_lat[k] = ev_lat[i]
                ev_lon[k] = ev_lon[i]
                ev_dep[k] = ev_dep[i]
                ev_mag[k] = ev_mag[i]
                ev_herr[k] = ev_herr[i]
                ev_zerr[k] = ev_zerr[i]
                ev_res[k] = ev_res[i]
                ev_cusp[k] = ev_cusp[i]
                k = k+1
            #else:
            #    #print('Event no data: ',ev_cusp[i])
        else:
            if len(amcusp)>0 and ev_cusp[i] in amcusp:
                #print('Airquake: ',ev_cusp[i])
                continue
            if ev_cusp[i] in dt_c1 or ev_cusp[i] in dt_c2:
                ev_date[k] = ev_date[i]
                ev_time[k] = ev_time[i]
                ev_lat[k] = ev_lat[i]
                ev_lon[k] = ev_lon[i]
                ev_dep[k] = ev_dep[i]
                ev_mag[k] = ev_mag[i]
                ev_herr[k] = ev_herr[i]
                ev_zerr[k] = ev_zerr[i]
                ev_res[k] = ev_res[i]
                ev_cusp[k] = ev_cusp[i]
                k = k+1
            #else:
            #    #print('Event no data: ',ev_cusp[i])
    nev=k
    ev_date = ev_date[:nev]
    ev_time = ev_time[:nev]
    ev_lat = ev_lat[:nev]
    ev_lon = ev_lon[:nev]
    ev_dep = ev_dep[:nev]
    ev_mag = ev_mag[:nev]
    ev_herr = ev_herr[:nev]
    ev_zerr = ev_zerr[:nev]
    ev_res = ev_res[:nev]
    ev_cusp = ev_cusp[:nev]
    #print('# Post clean event loop (nev): %10i' % nev)
    log.write('# Post clean event loop (nev): %10i \n' % nev)

    """
    Clean stations
    ---
    Remove any stations without associated data.
    """
    #print('Stations before: ',nsta)
    #sta_itmp1 = np.zeros(nsta)
    #sta_itmp1 = np.where(np.in1d(sta_lab,dt_sta1),1,0)
    #if reloctype==2 or reloctype==3:
    #   sta_itmp2 = np.zeros(nsta)
    #   sta_itmp2 = np.where(np.in1d(sta_lab,dt_sta2),1,0)
    k=int(0)
    #################### KB Note: Again I think it's entirely possible to remove loop. Revisit.
    for i in range(0,nsta):
        if reloctype==1:
            if sta_lab[i] in dt_sta1:
                sta_lab[k] = sta_lab[i]
                sta_lat[k] = sta_lat[i]
                sta_lon[k] = sta_lon[i]
                k=k+1
            #else:
            #    print('Bad station: ',sta_lab[i])
        else:
            if sta_lab[i] in dt_sta1 or sta_lab[i] in dt_sta2:
                sta_lab[k] = sta_lab[i]
                sta_lat[k] = sta_lat[i]
                sta_lon[k] = sta_lon[i]
                k=k+1
            #else:
            #    print('Bad station: ',sta_lab[i])

    nsta=k
    sta_lab = sta_lab[:nsta]
    sta_lat = sta_lat[:nsta]
    sta_lon = sta_lon[:nsta]
    log.write('# Post clean stations loop (nsta) = %6i \n' % nsta)
    #print('# Post clean stations loop (nsta) = %6i' % nsta)

    """
    Fast indexing arrays as implemented by Waldhauser 2001.
    ---
    This loop shouldn't break.  An error here signifies bad
    data values or organization or a bug in the code.
    """
    for i in range(ndt):
        """
        Double check each data has a station and event pair that
        are defined.
        ---
        If this breaks it means the event and station arrays do not
        match the data arrays.  This causes a fatal error.
        """
        try:
            dt_ista1[i] = np.where(sta_lab==dt_sta1[i])[0][0]
            if reloctype==2 or reloctype==3:
                dt_ista2[i] = np.where(sta_lab==dt_sta2[i])[0][0]
            dt_ic1[i] = np.where(ev_cusp==dt_c1[i])[0][0]
            if reloctype==1 or reloctype==3:
                dt_ic2[i] = np.where(ev_cusp==dt_c2[i])[0][0]
        except:
            raise Exception('FATAL ERROR INDEXING. Clean arrays.')
    """
    Return all the updated arrays
    """
    retlist = [nev,ev_cusp[:nev],ev_date[:nev],ev_time[:nev],ev_lat[:nev],ev_lon[:nev], 
               ev_dep[:nev],ev_mag[:nev],ev_herr[:nev],ev_zerr[:nev],ev_res[:nev],nsta,
               sta_lab[:nsta],sta_lat[:nsta],sta_lon[:nsta],dt_ic1[:ndt],dt_ic2[:ndt],
               dt_ista1[:ndt],dt_ista2[:ndt]]
    return retlist


def cluster(log,fol,nev,ndt,idata,minobs_cc,minobs_ct,dt_c1,dt_c2,ev_cusp):
    """
    Cluster Events
    :::
    Only called if minobs_cc and minobs_ct are > 0 depending
    on data type.
    :::
    Parameters:
    log (file obj) ---- Log file
    nev (int) ---- Number of events
    ndt (int) ---- Number of differential times
    idata (int) ---- Data indicator
                     0 - Synthetic Data
                     1 - Cross Correlations
                     2 - Catalog
                     3 - Cross and Catalog
                     From input file
    minobs_cc (int) ---- Min. # of CC obsverations to 
                         consider an event-pair
                         From Input File
    minobs_ct (int) ---- Min. # of catalog obsverations to 
                         consider an event-pair
                         From Input File
    dt_ic1[ndt] (int array) ---- Event 1 indexes in ev_cusp 
                                 for each data
    dt_ic2[ndt] (int array) ---- Event 2 indexes in ev_cusp 
                                 for each data
    ev_cusp[nev] (int array) ---- Event IDs
    :::
    Returns:
    clust[nclust,nev] (int array) ---- Cluster holds the number of events
                                       in each cluster and a list of events
                                       in each cluster
                                       ---
                                       Each column indicates a cluster of events
                                       - row 0 = # of events in cluster
                                       - row 1:end = list of event IDs in cluster
    noclust[nev] (int array) --- Holds IDs of events that are not clustered
    nclust (int) --- Number of clusters
    :::
    """
    log.write('\n\nClustering....\n')

    """
    Set up event-pair arrays
    ---
    Symmetical matrix (no duplicate pairs)
    """
    # apair_n holds the number of measurements per event pair
    apair_n = np.zeros((nev,nev),dtype='int8') 
    icusp = np.copy(ev_cusp[:nev])
    icusp = np.sort(icusp) # Sorted event IDS
    # Loop over measurements          
    for i in range(ndt):
        # Pull event IDs for each data
        j = np.argwhere(icusp==dt_c1[i])[0][0]
        k = np.argwhere(icusp==dt_c2[i])[0][0]
        # Add plus one to the counter
        # If the indexes are switched this is added to the lower matrix
        # Else added to the upper matrix
        # Only lower matrix will be counted
        apair_n[j,k] += 1
        apair_n[k,j] += 1

    # If synthetic data, only have 1 cluster.
    # If only one type of data, the other type will not impede clustering.
    if idata==0 or idata==1:
        minobs_ct=0
    if idata==0 or idata==2:
        minobs_cc=0
    
    """
    Cluster Events
    ---
    If event pairs share more than the minimum number of obs.
    then both events are clustered.
    All common pairs are then clustered.
    ---
    Clustering based solely on number of measurements.
    """
    # Initialize array acl to store cluster index for each event
    acl = nev*np.ones(nev,dtype='int') # Initialise all events to higher 
                                       # than largest possible cluster #
                                       # (Max. # of clusters == Nev)
    k=0 # Number of event-pairs < min # of obs
    n=0 # Cluster number
    # Loop over lower matrix event pairs
    for i in range(1,nev):
        for j in range(0,i):
            if apair_n[i,j]>=(minobs_cc+minobs_ct):
                # If the number of shared obs is over minimum
                # Then both events are in the same cluster
                if acl[i]<acl[j]:
                    # If event 1 is in a smaller number cluster 
                    # then move event 2 to new cluster.
                    if acl[j]==nev:
                        # If event 2 not clustered than add to
                        # event 1 cluster
                        acl[j] = acl[i]
                    else:
                        ia = acl[j] # Event 2 Cluster ID
                        #print('Ia: ',ia,' for Event ',j)
                        # Loop over previous events and move
                        # to event 1 cluster if in event 2
                        # cluster
                        for ii in range(nev):
                            if acl[ii]==ia:
                                acl[ii]=acl[i]
                elif acl[j]<acl[i]:
                    # If event 2 is in smaller number cluster
                    # then move event 2 to new cluster
                    if acl[i]==nev:
                        # If event 1 not clustered than add to
                        # event 2 cluster
                        acl[i] = acl[j]
                    else:
                        ia=acl[i] # Event 1 Cluster ID
                        #print('Ia: ',ia,' for Event ',i)
                        # Loop over previous events and move
                        # to event 2 cluster if in event 1
                        # cluster
                        for ii in range(nev):
                            if acl[ii]==ia:
                                acl[ii]=acl[j]
                elif acl[i]==nev and acl[j]==nev:
                    # If event 1 isn't clustered
                    # **(event 2 wouldn't be cluster either)
                    # Put both events in same cluster.
                    acl[i] = n
                    acl[j] = n
                    n=n+1
    
    """
    Store event keys for each cluster in clust
    """
    acli = np.argsort(acl) # Sorted array indices 
                           # Largest cluster first
    # Initialize both clust and noclust to largest
    # number of possible clusters (nev)
    # Clust holds IDs so it is 2D
    clust = np.zeros((nev,nev+1),dtype='uint32')
    noclust = np.zeros(nev,dtype='uint32')

    n = 0 # Cluster
    nn = 1 # Event in cluster
    iskip=0
    #Check if clustering occurred
    if acl[acli[0]]==nev:
        # If no events clustered
        i=1
        iskip=1
    else:
        clust[n,nn] = icusp[acli[0]]

    nclust=0
    if iskip==0:
        """ 
        If there was clustering 
        """
        # Loop over all events
        for i in range(1,nev):
            # Check if new cluster
            if acl[acli[i]]>acl[acli[i-1]]:
                clust[n,0] = nn # Save number of events in prev. cluster
                n = n+1 # Move to next cluster
                nn = 0 # Reset event list
            if acl[acli[i]]==nev: # If event not clustered
                iskip=1    # Indicate skipped event
                break       # Break out of for loop
            # Otherwise same cluster
            nn = nn+1 # Add to event count
            clust[n,nn] = icusp[acli[i]] # Save event name to cluster list 
        clust[n,0] = nn  # Save last cluster
        nclust = n #+1 # Save number of clusters
        noclust[0] = 0
    if iskip==1:
        """
        Skip events
        """
        for j in range(i,nev): # Loop over events
            noclust[j-i+1] = icusp[acli[j]] # Find event ID
        noclust[0] = nev-i # Save number of skipped events
    
    """
    Sort-biggest cluster first
    """
    if nclust>1:
        for i in range(nclust-1): # Loop over clusters i
            for j in range(i,nclust): # Loop over all other clusters j
                # Check if cluster i has fewer events than j to switch order
                # otherwise order is correct
                if clust[i,0] <= clust[j,0]:
                    # Move cluster i to end of cluster array
                    for k in range(clust[i,0]+1):
                        clust[nev-1,k] = clust[i,k]
                    # Move cluster j to cluster i
                    for k in range(0,clust[j,0]+1):
                        clust[i,k] = clust[j,k]
                    # Move cluster i to cluster j location
                    for k in range(clust[nev-1,0]+1):
                        clust[j,k] = clust[nev-1,k]
    """
    Count No. of Clustered Events
    """
    k = 0
    for i in range(nclust):
        k += clust[i,0]

    """
    Record Clusters
    """
    #print('> Clustered events: %5i' % k)
    #print('> Isolated events: %5i' % noclust[0])
    #print('> No. of clusters: %5i' % nclust)
    for i in range(0,nclust):
        print('Cluster %4i: %5i events' % (i+1,clust[i,0]))
    log.write('# Clustered events: %5i \n' % k)
    log.write('# Isolated events: %5i \n' % noclust[0])
    for i in range(0,int(noclust[0])):
        log.write(' %15s ' % noclust[i])
    log.write('# Clusters= %5i, for min. number of links set to %5i \n' % (nclust,minobs_ct+minobs_cc))
    for i in range(0,nclust):
        log.write('Cluster %5i: %5i events \n' % (i,clust[i,0]))
        for j in range(1,int(clust[i,0])+1):
            log.write('%15s ' % clust[i,j])
    log.write('\n\n')

    if nclust!=0:
        clust = clust[:nclust,:(clust[0,0]+1)]
        noclust = noclust[:(noclust[0]+1)]
    
    if nclust==0:
        # No clustering (one cluster includes all events)
        nclust = 1
        clust = np.zeros((1,nev+1))
        clust[0,0] = nev
        for i in range(0,nev):
            clust[0,i+1] = ev_cusp[i]
        #print('Nclust: ',nclust)
        #print('Noclust: ',len(noclust))
        #print('Clust[0]: ',clust[0,:])
        #raise Exception('No clusters.')

    # Save clusters
    for cl in range(nclust):
        count = 0
        cfilestr = os.path.join(fol,'clust_%i.dat' % cl)
        cfile = open(cfilestr,'w')
        cstring = ''
        for i in range(1,int(clust[cl,0])):
            if count==8:
                cstring += '\n'
                cfile.write(cstring)
                cstring = ''
                count=0
            cstring += ('%i ' % clust[cl,i])
            count+=1
        cfile.close()
    # Return the cluster array, the non-clustered events,
    # and the number of clusters
    return clust,noclust,nclust


def trialsrc(istart,sdc0_lat,sdc0_lon,sdc0_dep,nev,ev_cusp,ev_lat,ev_lon,ev_dep):
    """
    Set up source locations 
    :::
    If istart=1, all events start at cluster centroid
    If istart=2, all events start at catalog locations
    Lat/Lons are converted to cartesian x/y
    Cartesian coordinate system initialised from cluster centroid epicenter
    :::
    Parameters:
    istart (int) ---- Integer indicating initial source locations
    sdc0_lat float) ---- Starting cluster centroid latitude
    sdc0_lon (float) ---- Starting cluster centroid longitude
    sdc0_dep (float) ---- Starting cluster centroid depth
    nev (int) ---- Number of events
    ev_cusp[nev] (int array) ---- Numpy array of event IDs
    ev_lat[nev] (float array) ---- Numpy array of event latitudes
    ev_lon[nev] (float array) ---- Numpy array of event longitudes
    ev_dep[nev] (float array) ---- Numpy array of event depths
    :::
    Returns:
    nsrc (int) ---- Number of events
    src_cusp (int array) ---- Numpy array of event IDs
    src_lat0[nev] (float array) ---- Initial source latitudes
                                     if istart=1, this is cluster centroid latitude
                                     if istart=2, this is event latitudes
    src_lon0[nev] (float array) ---- Initial source longitudes
                                        ""
                                        ""
    src_x0[nev] (float array) ---- Initial source x
                                   if istart=1, this is zero (coordinate system origin)
                                   if istart=2, this is the event x
    src_y0[nev] (float array) ---- Initial source y
                                      ""
                                      ""
    src_z0[nev] (float array) ---- Initial source z
                                      ""
                                      ""
    src_t0[nev] (float array) ---- Initial source times
                                      ""
                                      ""
    src_lat[nev] (float array) ---- Current source lat (updated later)
    src_lon[nev] (float array) ---- Current source lon
    src_dep[nev] (float array) ---- Current source dep
    src_x[nev] (float array) ----  Current x (updated later)
    src_y[nev] (float array) ----  Current y
    src_z[nev] (float array) ----  Current z
    src_t[nev] (float array) ----  Current t
    src_xi[nev] (float array) ---- * ORIGINAL SOURCE LOCATIONS
    src_yi[nev] (float array) ---- * Event locations both times 
    src_zi[nev] (float array) ---- * Regardless of istart
    src_ti[nev] (float array) ---- * To compare absolute loc shifts
                                   * For synthetic model testing
    :::
    """

    """
    Copy common arrays
    """
    src_x = np.zeros(nev,dtype='float')
    src_y = np.zeros(nev,dtype='float')
    src_t = np.zeros(nev,dtype='float')
    src_ti = np.zeros(nev,dtype='float')
    src_t0 = np.zeros(nev,dtype='float')
    src_cusp = np.copy(ev_cusp)

    """
    Set up parameters for initial inversion
    """

    if istart==1:
        """
        Starting from centroid of cluster
        """
        """
        Cluster center as initial trial source
        """
        nsrc = 1
        src_lon = np.full((nev),sdc0_lon)
        src_lat = np.full((nev),sdc0_lat)
        src_dep = np.full((nev),sdc0_dep)
        src_lon0 = np.full((nev),sdc0_lon)
        src_lat0 = np.full((nev),sdc0_lat)

        """
        Set to cartesian origin (cluster centroid)
        """
        src_z = np.zeros(nev,dtype='float')
        src_x0 = np.zeros(nev,dtype='float')
        src_y0 = np.zeros(nev,dtype='float')
        src_z0 = np.zeros(nev,dtype='float')

        """
        Initial cartesian coordinates
        """
        try:
            true = np.loadtxt('origXY.loc')
            src_xi = np.copy(true[:,0]) 
            src_yi = np.copy(true[:,1])
            src_zi = np.copy(true[:,2])
        except:
            # Convert to cartesian
            for i in range(0,nev):
                [x,y] = sdc2(ev_lat[i],ev_lon[i],-1)
                src_xi[i] = x*1000.
                src_yi[i] = y*1000.
            src_z = np.full((nev),(ev_dep-sdc0_dep)*1000.)

    else:
        """
        Starting model are catalog source locations
        """
        nsrc=nev

        """
        Copy catalog locations
        """
        src_lon = np.copy(ev_lon)
        src_lat = np.copy(ev_lat)
        src_dep = np.copy(ev_dep)

        # Convert to cartesian
        for i in range(0,nev):
            [x,y] = sdc2(src_lat[i],src_lon[i],-1)
            src_x[i] = x*1000.
            src_y[i] = y*1000.
        src_z = np.full((nev),(ev_dep-sdc0_dep)*1000.)
        src_lon0 = np.copy(ev_lon)
        src_lat0 = np.copy(ev_lat)

        """
        Copy catalog initial locations
        """
        src_x0 = np.copy(src_x)
        src_y0 = np.copy(src_y)
        src_z0 = np.copy(src_z)

        try:
            true = np.loadtxt('origXY.loc')
            src_xi = np.copy(true[0]) 
            src_yi = np.copy(true[1])
            src_zi = np.copy(true[2])
        except:
            src_xi = np.copy(src_x) 
            src_yi = np.copy(src_y)
            src_zi = np.copy(src_z)
            
    return [nsrc,src_cusp,src_lat0,src_lon0,
            src_x0,src_y0,src_z0,src_t0,
            src_lat,src_lon,src_dep,
            src_x,src_y,src_z,src_t,
            src_xi,src_yi,src_zi,src_ti]


def dtres(reloctype,ndt,nsrc,dt_dt,dt_idx,dt_ista1,dt_ista2,dt_ic1,dt_ic2,src_t,
          tmp_ttp,tmp_tts):
    """
    Calculates difference vector dt_cal and double
    difference vector dt_dt
    :::
    Parameters:
    reloctype (int) --- Double-difference method type
    ndt (int) --- No. of data
    nsrc (int) --- No. of sources
    dt_dt[ndt] (float array) --- Data diff. time
    dt_idx[ndt] (int array) --- Data type indexes
    dt_ista1[ndt] (int array) --- Data station indexes
    dt_ista2[ndt] (int array) --- Data station indexes
    dt_ic1[ndt] (int array) --- Data event 1 indexes
    dt_ic2[ndt] (int array) --- Data event 2 indexes
    src_t[nev] (float array) --- Source times
    tmp_ttp[nsta,nev] (float array) --- P raytracing times
    tmp_tts[nsta,nev] (float array) --- S raytracing times
    :::
    Returns:
    dt_cal[ndt] (float array) --- Data cal. diff times
    dt_res[ndt] (float arrya) --- Data double-diff. times
    :::
    """

    """
    Declare Arrays
    """
    dt_res = np.zeros(ndt,dtype='float')
    dt_cal = np.zeros(ndt,dtype='float')
    tt = []

    if nsrc==1:
        """
        Single source
        """
        dt_res = np.copy(dt_dt)
    else:
        """
        Multiple sources
        """

        if reloctype==1:
            [dt_cal,tt1,tt2] = dtcal_ev(ndt,dt_cal,dt_dt,dt_idx,dt_ista1,dt_ic1,dt_ic2,
                                        src_t,tmp_ttp,tmp_tts)
            tt.append(tt1)
            tt.append(tt2)
        elif reloctype==2:
            [dt_cal,tt1,tt2] = dtcal_st(ndt,dt_cal,dt_dt,dt_idx,dt_ista1,dt_ista2,dt_ic1,
                                        src_t,tmp_ttp,tmp_tts)
            tt.append(tt1)
            tt.append(tt2)
        elif reloctype==3:  
            [dt_cal,tt1,tt2,tt3,tt4] = dtcal_doub(ndt,dt_cal,dt_dt,dt_idx,dt_ista1,dt_ista2,
                                                  dt_ic1,dt_ic2,src_t,tmp_ttp,tmp_tts)
            tt.append(tt1) 
            tt.append(tt2)
            tt.append(tt3)
            tt.append(tt4)

        """
        Calculate double difference
        """
        dt_res = dt_dt - dt_cal
        tt = np.asarray(tt,dtype='float')

    return dt_cal,dt_res,tt


def weighting(log,reloctype,ndt,mbad,amcusp,idata,kiter,ineg,maxres_cross,maxres_net,
              maxdcc,maxdct,minwght,wt_ccp,wt_ccs,wt_ctp,wt_cts,dt_c1,dt_c2,dt_idx,
              dt_qual,dt_res,dt_offse,dt_offss):
    """
    Weighting/re-weighting function.
    :::
    PARAMETERS:
    log (file obj) ---- Log file
    reloctype (int) ---- Double-difference pairing type
    ndt (int) ---- Number of data
    mbad (int) --- Airquake counter
    amcusp (int array) --- Array of airquake IDS
    idata (int) --- Data type switch
    kiter (int) --- Skipped data counter
    ineg (int) --- If skip function is triggered
    maxres_cross (float) --- Max. res. allowed for cc pairs
    maxres_net (float) --- Max. res. allowed for ct pairs
    maxdcc (float) --- Max. interevent dist. for cc event-pairs
    maxdct (float) --- Max. interevent dist. for ct event-pairs
    minwght (float) --- Min. acceptable data weight
    wt_ccp (float) --- Inversion weight for P cc data
    wt_ccs (float) --- Inversion weight for S cc data
    wt_ctp (float) --- Inversion weight for P ct data
    wt_cts (float) --- Inversion weight for S ct data
    dt_c1 (int array) --- Array of event 1 IDs
    dt_c2 (int array) --- Array of event 2 IDs
    dt_idx (int array) --- Array of data type indexes
    dt_qual (float array) --- Array of a prioi data weight
    dt_res (float array) --- Array of double-pair residuals
    dt_offse (float array) --- Array of interevent distances
    dt_offss (float array) --- Array of interstation distances
    reweight_switch (int) --- On/off switch for reweighting
    :::
    RETURNS:
    ineg (int) --- Switch indicating if data needs to be skipped
                   0 if no data to skip
                   1 if there is data with weight less than min. weight
                        this data will be removed in the skip data function
    dt_wt (float array) --- Array containing updated/reweighted data weights
    :::
    """

    """
    Get a priori data weights
    ---
    Pull and weight data based on inversion data type weighting and a priori weights
    All the quality transfer is done in getdata
    """
    ineg = 0 # flag if neg weights exist
    dt_wt = np.zeros(ndt,dtype='float')
    for i in range(ndt): 
        if dt_idx[i]==1:
            dt_wt[i] = wt_ccp*dt_qual[i]
        elif dt_idx[i]==2:
            dt_wt[i] = wt_ccs*dt_qual[i]
        elif dt_idx[i]==3:
            dt_wt[i] = wt_ctp*dt_qual[i]
        elif dt_idx[i]==4:
            dt_wt[i] = wt_cts*dt_qual[i]

        #for j in range(0,mbad):
        if reloctype==2:
            if dt_c1[i] in amcusp:
                dt_wt[i]=0.
                ineg=1
        else:
            if (dt_c1[i] in amcusp) or (dt_c2[i] in amcusp):
                dt_wt[i]=0.
                ineg=1
    """
    Reweighting
    ---
    Biweight function (Mosteller and Tukey 1977)
    """
    datet = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    dt_tmp = np.zeros(ndt)

    # This is triggered if data is:
    #       (a) Cross correlation data and a max res. cc or max. dist cc is specified in inversion instructions
    #       (b) Catalog data and a max res. ct or max. dist. ct is specified in inversion instructions
    if reloctype!=2:
        if ((idata==1 or idata==3) and (maxres_cross!=-9 or maxdcc!=-9)) or ((idata==2 or idata==3) and (maxres_net!=-9 or maxdct!=-9)):
            log.write('\nReweighting ... %s \n' % datet)

            """
            Get median and MAD of data residuals
            """
            if idata==3:    # If both cc and ct are used
                if maxres_cross>=1: # And there is a residual cutoff
                    """
                    Cross-correlation data
                    """
                    k=0
                    for i in range(ndt):
                        if dt_idx[i]<=2:
                            dt_tmp[k] = dt_res[i]
                            k = k+1
                    med_cc = np.median(dt_tmp[:k])
                    dt_tmp = np.abs(dt_tmp - med_cc)
                    mad_cc = np.median(dt_tmp[:k])
                    mad_cc = mad_cc/0.67449 # MAD for gaussian
                if maxres_net >= 1: # And there is a residual cutoff
                    """
                    Catalog data
                    """
                    k=0
                    for i in range(ndt):
                        if dt_idx[i]>=3:
                            dt_tmp[k] = dt_res[i]
                            k = k+1
                    med_ct = np.median(dt_tmp[:k])
                    dt_tmp = np.abs(dt_tmp - med_ct)
                    mad_ct = np.median(dt_tmp[:k])
                    mad_ct = mad_ct/0.67449 # MAD for gaussian
            elif (idata==1 and maxres_cross>=1) or (idata==2 and maxres_net>=1):    # One data type only
                dt_tmp = np.copy(dt_res)
                med_cc = np.median(dt_tmp[:ndt])
                dt_tmp = np.abs(dt_tmp-med_cc)
                mad_cc = np.median(dt_tmp[:ndt])
                mad_cc = mad_cc/0.67449
                if idata==2:    # If data type is catalog just transfer value to new name
                    mad_ct = mad_cc

            """
            Define residual cutoffs
            """
            maxres_cc = maxres_cross # absolute cutoff value
            maxres_ct = maxres_net # absolute cutoff value
            if maxres_cross >= 1:
                maxres_cc = mad_cc*maxres_cross
            if maxres_net >= 1:
                maxres_ct = mad_ct*maxres_net

            """
            Apply residual/offset dependent weights to a priori weights
            ---
            Factor of 2 difference in cc and ct weighting due to expected noise if
            there is a distance cutoff (cc differences need to be small for best coherency)
            """
            nncc = 0
            nnct = 0
            ncc = 0
            nct = 0
            for i in range(ndt):  # Loop over data
                if dt_idx[i]<=2:    # For cross-correlation data
                    ncc += 1        # Index cc data counter

                    """
                    bi ^5 offset weighting for cross data:
                    ---
                    exp needs to be uneven so weights become
                    negative for offsets greater than x km
                    """
                    if maxdcc!=-9:  # If there is a distance cutoff
                        dt_wt[i] = dt_wt[i]*(1-(dt_offse[i]/(maxdcc*1000))**5)**5
                        #elif reloctype==2:
                        #    dt_wt[i] = dt_wt[i]*(1-(dt_offss[i]/(maxdcc*1000))**5)**5
                        #elif reloctype==3:
                        #    dt_wt[i] = dt_wt[i]*(1-(dt_offse[i]/(maxdcc*1000))**5)**5
                    """
                    bi-cube residual weighting
                    ---
                    needs to be a cube so that res > cutoff
                    became negative
                    """
                    if maxres_cross>0 and dt_wt[i]>0.000001:
                        dt_wt[i] = dt_wt[i]*(1-(np.abs(dt_res[i])/maxres_cc)**3)**3
                    if dt_wt[i]<minwght:
                        nncc = nncc+1
                else:               # Catalog data
                    nct += 1        # Index ct data counter

                    """
                    bicube offset weighting for catalog data:
                    """
                    if maxdct!=-9:
                        dt_wt[i] = dt_wt[i]*(1-(dt_offse[i]/(maxdct*1000))**3)**3
                        #elif reloctype==2:
                        #    dt_wt[i] = dt_wt[i]*(1-(dt_offss[i]/(maxdct*1000))**3)**3
                        #elif reloctype==3:
                        #    dt_wt[i] = dt_wt[i]*(1-(dt_offse[i]/(maxdct*1000))**3)**3
                    """
                    bi-cube residual weighting
                    """
                    if dt_wt[i]>0.000001 and maxres_net>0:
                        dt_wt[i] = dt_wt[i]*(1-(np.abs(dt_res[i])/maxres_ct)**3)**3
                    if dt_wt[i]<minwght:
                        nnct = nnct+1

            """
            Check if there is data to be skipped
            ---
            If data weight < minwght, the skip data function will be called
            And data with small weights will be removed
            """
            ineg = 0
            for j in range(ndt):
                if dt_wt[j]<minwght:
                    ineg=1
                    break

            """
            Write summary to log
            """
            if idata==1 or idata==3:
                log.write('cc res/dist cutoff: %7.3f s/ %5.2f km (%5.1f%%)\n' %
                          (maxres_cc,maxdcc,(nncc*100./ncc)))
            if idata==2 or idata==3:
                log.write('ct res/dist cutoff: %7.3f s/ %5.2f km (%5.1f%%)\n' %
                          (maxres_ct,maxdct,(nnct*100./nct)))

    elif reloctype==2:
        if ((idata==1 or idata==3) and (maxres_cross!=-9)) or ((idata==2 or idata==3) and (maxres_net!=-9)):
            log.write('\nReweighting ... %s \n' % datet)

            """
            Get median and MAD of data residuals
            """
            if idata==3:    # If both cc and ct are used
                if maxres_cross>=1: # And there is a residual cutoff
                    """
                    Cross-correlation data
                    """
                    k=0
                    for i in range(ndt):
                        if dt_idx[i]<=2:
                            dt_tmp[k] = dt_res[i]
                            k = k+1
                    med_cc = np.median(dt_tmp[:k])
                    dt_tmp = np.abs(dt_tmp - med_cc)
                    mad_cc = np.median(dt_tmp[:k])
                    mad_cc = mad_cc/0.67449 # MAD for gaussian
                if maxres_net >= 1: # And there is a residual cutoff
                    """
                    Catalog data
                    """
                    k=0
                    for i in range(ndt):
                        if dt_idx[i]>=3:
                            dt_tmp[k] = dt_res[i]
                            k = k+1
                    med_ct = np.median(dt_tmp[:k])
                    dt_tmp = np.abs(dt_tmp - med_ct)
                    mad_ct = np.median(dt_tmp[:k])
                    mad_ct = mad_ct/0.67449 # MAD for gaussian
            elif (idata==1 and maxres_cross>=1) or (idata==2 and maxres_net>=1):    # One data type only
                dt_tmp = np.copy(dt_res)
                med_cc = np.median(dt_tmp[:ndt])
                dt_tmp = np.abs(dt_tmp-med_cc)
                mad_cc = np.median(dt_tmp[:ndt])
                mad_cc = mad_cc/0.67449
                if idata==2:    # If data type is catalog just transfer value to new name
                    mad_ct = mad_cc

            """
            Define residual cutoffs
            """
            maxres_cc = maxres_cross # absolute cutoff value
            maxres_ct = maxres_net # absolute cutoff value
            if maxres_cross >= 1:
                maxres_cc = mad_cc*maxres_cross
            if maxres_net >= 1:
                maxres_ct = mad_ct*maxres_net

            """
            Apply residual/offset dependent weights to a priori weights
            ---
            Factor of 2 difference in cc and ct weighting due to expected noise if
            there is a distance cutoff (cc differences need to be small for best coherency)
            """
            nncc = 0
            nnct = 0
            ncc = 0
            nct = 0
            for i in range(ndt):  # Loop over data
                if dt_idx[i]<=2:    # For cross-correlation data
                    ncc += 1        # Index cc data counter
                    """
                    bi-cube residual weighting
                    ---
                    needs to be a cube so that res > cutoff
                    became negative
                    """
                    if maxres_cross>0 and dt_wt[i]>0.000001:
                        dt_wt[i] = dt_wt[i]*(1-(np.abs(dt_res[i])/maxres_cc)**3)**3
                    if dt_wt[i]<minwght:
                        nncc = nncc+1
                else:               # Catalog data
                    nct += 1        # Index ct data counter
                    """
                    bi-cube residual weighting
                    """
                    if dt_wt[i]>0.000001 and maxres_net>0:
                        dt_wt[i] = dt_wt[i]*(1-(np.abs(dt_res[i])/maxres_ct)**3)**3
                    if dt_wt[i]<minwght:
                        nnct = nnct+1

            """
            Check if there is data to be skipped
            ---
            If data weight < minwght, the skip data function will be called
            And data with small weights will be removed
            """
            ineg = 0
            for j in range(ndt):
                if dt_wt[j]<minwght:
                    ineg=1
                    break

            """
            Write summary to log
            """
            if idata==1 or idata==3:
                log.write('cc res/dist cutoff: %7.3f s/ (%5.1f%%)\n' %
                          (maxres_cc,(nncc*100./ncc)))
            if idata==2 or idata==3:
                log.write('ct res/dist cutoff: %7.3f s/ (%5.1f%%)\n' %
                          (maxres_ct,(nnct*100./nct)))

    # If there's data to skip, mark as iteration with skipped data
    if ineg==0:
        kiter+=1
    # Return the ineg switch and the weight array
    return [ineg,dt_wt,kiter]


def skip(log,reloctype,kiter,minwght,ndt,nev,nsrc,nsta,ev_cusp,ev_date,ev_time,ev_mag, 
         ev_lat,ev_lon,ev_dep,ev_x,ev_y,ev_z,ev_herr,ev_zerr,ev_res,src_cusp,src_lat,
         src_lon,src_dep,src_lat0,src_lon0,src_x,src_y,src_z,src_t,src_x0,src_y0,src_z0,
         src_t0,src_xi,src_yi,src_zi,src_ti,sta_lab,sta_lat,sta_lon,sta_rmsc,sta_rmsn,
         sta_np,sta_ns,sta_nnp,sta_nns,dt_sta1,dt_sta2,dt_c1,dt_c2,dt_idx,dt_dt,dt_qual,
         dt_cal,dt_ista1,dt_ista2,dt_ic1,dt_ic2,dt_res,dt_wt,dt_offse,dt_offss,tmp_ttp,
         tmp_tts,tmp_xp,tmp_yp,tmp_zp,nct,ncc,amcusp,mbad):
    """
    Skip outliers and airquakes
    :::
    PARAMETERS:
    log (file obj) ---- Log file
    kiter (int) --- Iteration number
    minwght (float) --- Minimum data weight accepted
    [Pass in all other data arrays and counts.]
    :::
    RETURNS:
    retlist (list) --- List of all updated arrays
    :::
    """
    log.write('\n\nSkipping data...\n')

    """
    Skip data with small weights
    """
    # Record all counts
    ndtold = ndt
    nccold = int(0)
    nctold = int(0)
    nevold = nev
    nstaold = nsta
    ncc = int(0)
    nct = int(0)
    j = int(0)
    nccold = (dt_idx<=2).sum()
    nctold = (dt_idx>2).sum()

    """
    Exclude tmp_zp differences near
    zero (station-pairs only)
    """
    if reloctype==2:
        for i in range(nev):
            #try:
            mask = (dt_ic1==i)
            ic1 = dt_ic1[mask]
            ista1 = dt_ista1[mask]
            ista2 = dt_ista2[mask]
            zp1 = tmp_zp[ista1,ic1]
            zp2 = tmp_zp[ista2,ic1]
            zp = zp1-zp2
            zpsum = np.sum(zp**2)
            #zpsum = zpsum/ndt

            if zpsum==0.:
                import pdb; pdb.set_trace()
                log.write('Bad event (all angles equal) - Event index: %i \n' % i)
                print('Bad event (all angles equal) - Event index: %i' % i)
                #badindx = np.argwhere(amcusp==0)[0][0]
                amcusp[mbad]=ev_cusp[i] 
                mbad+=1
            #except:
            #    import pdb; pdb.set_trace()

    if reloctype==3:
        for i in range(nev):
            #try:
            #import pdb; pdb.set_trace()

            mask = (dt_ic1==i)
            ic1 = dt_ic1[mask]
            ista1 = dt_ista1[mask]
            ista2 = dt_ista2[mask]
            zp1 = tmp_zp[ista1,ic1]
            zp2 = tmp_zp[ista2,ic1]
            zp = zp1-zp2
            zpsum = np.sum(zp**2)
            count1 = np.count_nonzero(zp)

            mask2 = (dt_ic2==i)
            ic2 = dt_ic2[mask2]
            ista3 = dt_ista1[mask2]
            ista4 = dt_ista2[mask2]
            zp3 = tmp_zp[ista3,ic2]
            zp4 = tmp_zp[ista4,ic2]
            zpb = zp3-zp4
            
            zpsum2 = np.sum(zpb**2)
            count2 = np.count_nonzero(zpb)
            #zpsum = zpsum/ndt  

            zpsumtot = (zpsum+zpsum2)/ndt      

            # if i==554:
            #     log.write('Bad event (all angles equal) - Event index: %i \n' % i)
            #     print('Bad event (all angles equal) - Event index: %i' % i)
            #     #badindx = np.argwhere(amcusp==0)[0][0]
            #     amcusp[mbad]=ev_cusp[i] 
            #     mbad+=1

            #if zpsumtot<np.finfo(np.float).eps:
            if (zpsum==0. and zpsum2==0.) or ((count1+count2)<=1):
                log.write('Bad event (all angles equal) - Event index: %i \n' % i)
                print('Bad event (all angles equal) - Event index: %i' % i)
                #badindx = np.argwhere(amcusp==0)[0][0]
                amcusp[mbad]=ev_cusp[i] 
                mbad+=1
            #except:
            #    import pdb; pdb.set_trace()

    # if reloctype==3:
    #     for i in range(nev):
    #         if i==554:
    #             import pdb; pdb.set_trace()
    #         try:
    #             mask1 = (dt_ic1==i)
    #             #mask1 = (dt_ic1==i)
    #             #mask2 = (dt_ic2==i)
    #             ic1 = dt_ic1[mask]
    #             ic2 = dt_ic2[mask]
    #             ista3 = dt_ista1[mask]
    #             ista4 = dt_ista2[mask]
    #             zp1a = tmp_zp[ista3,ic1]
    #             zp2a = tmp_zp[ista4,ic1]
    #             zp1b = tmp_zp[ista3,ic2]
    #             zp2b = tmp_zp[ista4,ic2]
    #             zpb = (zp1a-zp2a) - (zp1b-zp2b)
    #             zpsum = np.sum(zpb**2)
    #             if zpsum==0.:
    #                 log.write('Bad event (all angles equal) - Event index: %i \n' % i)
    #                 #badindx = np.argwhere(amcusp==0)[0][0]
    #                 amcusp[mbad]=ev_cusp[i] 
    #                 mbad+=1
    #         except:
    #             import pdb; pdb.set_trace()
    #             log.write('Bad event (all angles equal) - Event index: %i \n' % i)
    #             #badindx = np.argwhere(amcusp==0)[0][0]
    #             amcusp[mbad]=ev_cusp[i] 
    #             mbad+=1

    """
    Mask of data above minwght
    ---
    Trim all other arrays
    """
    mask = (dt_wt>=minwght)
    dt_sta1 = dt_sta1[mask]
    if reloctype!=1:
        dt_sta2 = dt_sta2[mask]
        dt_offss = dt_offss[mask]
    dt_c1 = dt_c1[mask]
    if reloctype!=2:
        dt_c2 = dt_c2[mask]
        dt_offse = dt_offse[mask]
    dt_idx = dt_idx[mask]
    dt_qual = dt_qual[mask]
    dt_dt = dt_dt[mask]
    dt_cal = dt_cal[mask]
    dt_res = dt_res[mask]
    dt_wt = dt_wt[mask]
    ndt = len(dt_dt)

    """
    Exclude airquake data
    """
    k=0
    for i in range(ndt):
        if reloctype!=2:
            if dt_c1[i] in amcusp or dt_c2[i] in amcusp:
                continue   
            else: 
                dt_sta1[k] = dt_sta1[i]
                if reloctype==3:
                    dt_sta2[k] = dt_sta2[i]
                    dt_offss[k] = dt_offss[i]
                dt_offse[k] = dt_offse[i]
                dt_c1[k] = dt_c1[i]
                dt_c2[k] = dt_c2[i]
                dt_idx[k] = dt_idx[i]
                dt_qual[k] = dt_qual[i]
                dt_dt[k] = dt_dt[i]
                dt_cal[k] = dt_cal[i]
                dt_res[k] = dt_res[i]
                dt_wt[k] = dt_wt[i]
                k+=1
        else:
            if dt_c1[i] in amcusp:
                continue
            else:
                dt_sta1[k] = dt_sta1[i]
                dt_sta2[k] = dt_sta2[i]
                dt_offss[k] = dt_offss[i]
                dt_c1[k] = dt_c1[i]
                dt_idx[k] = dt_idx[i]
                dt_qual[k] = dt_qual[i]
                dt_dt[k] = dt_dt[i]
                dt_cal[k] = dt_cal[i]
                dt_res[k] = dt_res[i]
                dt_wt[k] = dt_wt[i]
                k+=1
    ndttmp = k 
    
    """
    Trim arrays
    """
    dt_sta1 = dt_sta1[:ndttmp].copy()
    dt_c1 = dt_c1[:ndttmp].copy()
    dt_idx = dt_idx[:ndttmp].copy()
    dt_qual = dt_qual[:ndttmp].copy()
    dt_dt = dt_dt[:ndttmp].copy()
    dt_cal = dt_cal[:ndttmp].copy()
    dt_res = dt_res[:ndttmp].copy()
    dt_wt = dt_wt[:ndttmp].copy()
    if reloctype!=1:
        dt_sta2 = dt_sta2[:ndttmp].copy()
        dt_offss = dt_offss[:ndttmp].copy()
    if reloctype!=2:
        dt_c2 = dt_c2[:ndttmp].copy()
        dt_offse = dt_offse[:ndttmp].copy()
    # New data counts
    ncc = (dt_idx<=2).sum()
    nct = (dt_idx>2).sum()
    ndt = ncc+nct
    if ndt!=ndttmp:
        raise Exception('(skip function): Airquake filtering')

    """
    Record no. of data removed/kept
    """
    log.write('# obs = %9i (%5.1f%%)\n' % (ndt,(ndt*100./ndtold)))
    if nccold>0. and nctold>0:
        log.write('# obs cc = %9i (%5.1f%%)\n' % (ncc,(ncc*100./nccold)))
        log.write('# obs ct = %9i (%5.1f%%)\n' % (nct,(nct*100./nctold)))

    """
    Clean station arrays not in cleanarrays
    """
    k=0
    for i in range(nsta):
        if reloctype==1:
            if sta_lab[i] in dt_sta1:
                sta_np[k] = sta_np[i]
                sta_ns[k] = sta_ns[i]
                sta_nnp[k] = sta_nnp[i]
                sta_nns[k] = sta_nns[i]
                sta_rmsc[k] = sta_rmsc[i]
                sta_rmsn[k] = sta_rmsn[i]
                # Update raytracing 2d arrays
                tmp_ttp[k,:] = tmp_ttp[i,:]
                tmp_tts[k,:] = tmp_tts[i,:]
                tmp_xp[k,:] = tmp_xp[i,:]
                tmp_yp[k,:] = tmp_yp[i,:]
                tmp_zp[k,:] = tmp_zp[i,:]
                k=k+1
            #else:
            #    print('StaRes Skip index: ',i)
        else:
            if sta_lab[i] in dt_sta1 or sta_lab[i] in dt_sta2:
                sta_np[k] = sta_np[i]
                sta_ns[k] = sta_ns[i]
                sta_nnp[k] = sta_nnp[i]
                sta_nns[k] = sta_nns[i]
                sta_rmsc[k] = sta_rmsc[i]
                sta_rmsn[k] = sta_rmsn[i]
                # Update raytracing 2d arrays
                tmp_ttp[k,:] = tmp_ttp[i,:]
                tmp_tts[k,:] = tmp_tts[i,:]
                tmp_xp[k,:] = tmp_xp[i,:]
                tmp_yp[k,:] = tmp_yp[i,:]
                tmp_zp[k,:] = tmp_zp[i,:]
                k=k+1
            #else:
            #    print('StaRes Skip index: ',i)
    nstatmp = k        # Update station number
    """
    Trim arrays
    """
    sta_np = sta_np[:nstatmp]
    sta_ns = sta_ns[:nstatmp]
    sta_nnp = sta_nnp[:nstatmp]
    sta_nns = sta_nns[:nstatmp]
    sta_rmsc = sta_rmsc[:nstatmp]
    sta_rmsn = sta_rmsn[:nstatmp]
    tmp_ttp = tmp_ttp[:nstatmp,:]
    tmp_tts = tmp_tts[:nstatmp,:]
    tmp_xp = tmp_xp[:nstatmp,:]
    tmp_yp = tmp_yp[:nstatmp,:]
    tmp_zp = tmp_zp[:nstatmp,:]

    """
    Update src arrays (starting location arrays)s
    """
    if nsrc!=1:
        k=0
        for i in range(0,nsrc):
            if reloctype!=2:
                # Skip if Airquake
                if src_cusp[i] in amcusp:
                    #print('Src Airquake index: ',i)
                    continue
                # The event ID has to be in one of the data arrays
                if src_cusp[i] in dt_c1 or src_cusp[i] in dt_c2:
                    # Update source arrays
                    src_cusp[k] = src_cusp[i]
                    src_lat[k] = src_lat[i]
                    src_lon[k] = src_lon[i]
                    src_lat0[k] = src_lat0[i]
                    src_lon0[k] = src_lon0[i]
                    src_dep[k] = src_dep[i]
                    src_x[k] = src_x[i]
                    src_y[k] = src_y[i]
                    src_z[k] = src_z[i]
                    src_t[k] = src_t[i]
                    src_x0[k] = src_x0[i]
                    src_y0[k] = src_y0[i]
                    src_z0[k] = src_z0[i]
                    src_t0[k] = src_t0[i]
                    src_xi[k] = src_xi[i]
                    src_yi[k] = src_yi[i]
                    src_zi[k] = src_zi[i]
                    src_ti[k] = src_ti[i]
                    # Update raytracing arrays
                    tmp_ttp[:,k] = tmp_ttp[:,i]
                    tmp_tts[:,k] = tmp_tts[:,i]
                    tmp_xp[:,k] = tmp_xp[:,i]
                    tmp_yp[:,k] = tmp_yp[:,i]
                    tmp_zp[:,k] = tmp_zp[:,i]
                    k += 1
                #else:
                #    print('Src skip index: ',i,' id: ',src_cusp[i])
            else:
                # Check if airquake
                if src_cusp[i] in amcusp:
                    #print('Src Airquake index: ',i)
                    continue
                if src_cusp[i] in dt_c1:
                    # Update event arrays
                    src_cusp[k] = src_cusp[i]
                    src_lat[k] = src_lat[i]
                    src_lon[k] = src_lon[i]
                    src_lat0[k] = src_lat0[i]
                    src_lon0[k] = src_lon0[i]
                    src_dep[k] = src_dep[i]
                    src_x[k] = src_x[i]
                    src_y[k] = src_y[i]
                    src_z[k] = src_z[i]
                    src_t[k] = src_t[i]
                    src_x0[k] = src_x0[i]
                    src_y0[k] = src_y0[i]
                    src_z0[k] = src_z0[i]
                    src_t0[k] = src_t0[i]
                    src_xi[k] = src_xi[i]
                    src_yi[k] = src_yi[i]
                    src_zi[k] = src_zi[i]
                    src_ti[k] = src_ti[i]
                    # Update raytracing arrays
                    tmp_ttp[:,k] = tmp_ttp[:,i]
                    tmp_tts[:,k] = tmp_tts[:,i]
                    tmp_xp[:,k] = tmp_xp[:,i]
                    tmp_yp[:,k] = tmp_yp[:,i]
                    tmp_zp[:,k] = tmp_zp[:,i]
                    k += 1
                #else:
                #    print('Src skip index: ',i)
        nsrc = k

    """
    Trim arrays
    """
    src_cusp = src_cusp[:nsrc]
    src_lat = src_lat[:nsrc]
    src_lon = src_lon[:nsrc]
    src_lat0 = src_lat0[:nsrc]
    src_lon0 = src_lon0[:nsrc]
    src_dep = src_dep[:nsrc]
    src_x = src_x[:nsrc]
    src_y = src_y[:nsrc]
    src_z = src_z[:nsrc]
    src_t = src_t[:nsrc]
    src_x0 = src_x0[:nsrc]
    src_y0 = src_y0[:nsrc]
    src_z0 = src_z0[:nsrc]
    src_t0 = src_t0[:nsrc]
    src_xi = src_xi[:nsrc]
    src_yi = src_yi[:nsrc]
    src_zi = src_zi[:nsrc]
    src_ti = src_ti[:nsrc]
    # Update raytracing arrays
    tmp_ttp = tmp_ttp[:,:nsrc]
    tmp_tts = tmp_tts[:,:nsrc]
    tmp_xp = tmp_xp[:,:nsrc]
    tmp_yp = tmp_yp[:,:nsrc]
    tmp_zp = tmp_zp[:,:nsrc]

    [nev,ev_cusp,ev_date,ev_time,ev_lat,ev_lon,ev_dep,ev_mag,ev_herr,ev_zerr,
     ev_res,nsta,sta_lab,sta_lat,sta_lon,dt_ic1,dt_ic2,dt_ista1,
     dt_ista2] = cleanarrays(log,reloctype,nev,ev_cusp,ev_date,ev_time,ev_lat,ev_lon,
                             ev_dep,ev_mag,ev_herr,ev_zerr,ev_res,nsta,sta_lab,sta_lat,
                             sta_lon,ndt,dt_c1,dt_c2,dt_sta1,dt_sta2,dt_ic1,dt_ic2,
                             dt_ista1,dt_ista2,amcusp)

    """
    Break if data cleared incorrectly
    """
    if nsta!=nstatmp:
        raise Exception('Inconsistent station data cleared in skip.')
    if nsrc!=nev:
        raise Exception('Inconsistent event data cleared in skip.')

    """
    Return all updated counts and data arrays
    """
    log.write('Skip function. No data deleted = %i \n' % (ndtold-ndt))
    log.write('Skip function. No events skipped = %i \n' % (nevold-nev))
    log.write('Skip function. No stations removed = %i \n' % (nstaold-nsta))
    retlist = [ndt,nev,nsrc,nsta,ev_cusp,ev_date,ev_time,ev_mag, ev_lat,ev_lon,ev_dep,ev_x,
               ev_y,ev_z,ev_herr,ev_zerr,ev_res,src_cusp,src_lat,src_lon,src_dep,src_lat0,
               src_lon0,src_x,src_y,src_z,src_t,src_x0,src_y0,src_z0,src_t0,src_xi,src_yi,
               src_zi,src_ti,sta_lab,sta_lat,sta_lon,sta_rmsc,sta_rmsn,sta_np,sta_ns,sta_nnp,
               sta_nns,dt_sta1,dt_sta2,dt_c1,dt_c2,dt_idx,dt_dt,dt_qual,dt_cal,dt_ista1,dt_ista2,
               dt_ic1,dt_ic2,dt_res,dt_wt,dt_offse,dt_offss,tmp_ttp,tmp_tts,tmp_xp,tmp_yp,
               tmp_zp,nct,ncc]
    return retlist


def resstat(log,reloctype,idata,ndt,nev,d,w,idx,rms_cc,rms_ct,rms_cc0,rms_ct0,rms_ccold,
            rms_ctold,rms_cc0old,rms_ct0old,dum=-999):
    """
    Calculate residual statistics
    :::
    Parameters:
    log (file obj) ---- Log file
    idata (int) --- Data type switch
    ndt (int) --- No. of data
    nev (int) --- No. of events
    d (float array) --- Residual Data
    w (float array) --- Data weights
    idx (int array) --- Data type indexes
    rms_cc (float) --- Previous iteration RMS cc residual
    rms_ct (float) --- Previous iteration RMS ct residual
    rms_cc0 (float) --- Previous iteration weighted RMS cc residual
    rms_ct0 (float) --- Previous iteration weighted RMS ct residual
    rms_ccold (float) --- Previous old RMS cc residual
    rms_ctold (float) --- Previous old RMS ct residual
    rms_cc0old (float) --- Previous old RMS cc residual
    rms_ct0old (float) --- Previous old RMS ct residual
    dum (float) --- Previous iteration average residual
    :::
    Returns:
    rms_cc (float) --- Current iteration RMS CC residual
    rms_ct (float) --- Current iteration RMS CT residual
    rms_cc0 (float) --- Current iteration weighted RMS CC residual
    rms_ct0 (float) --- Current iteration weighted RMS CT residual
    rms_ccold (float) --- Previous iteration RMS CC residual
    rms_ctold (float) --- Previous iteration RMS CT residual
    rms_cc0old (float) --- Previous iteration weighted RMS CC residual
    rms_ct0old (float) --- Previous iteration weighted RMS CT residual
    dum (float) --- This iterations avg. residual
    :::
    """

    """
    Save previous iter. rms values to old variables
    """
    rms_cc0old = rms_cc0
    rms_ct0old = rms_ct0
    rms_ccold = rms_cc
    rms_ctold = rms_ct 

    """
    Sum and scale weights over no. of data
    """
    j=0
    sw_cc=0.
    sw_ct=0.
    for i in range(0,ndt):
        if idx[i]<=2:
            sw_cc += w[i]
            j += 1
        else:
            sw_ct += w[i]
    if idata!=2:
        f_cc=j/sw_cc          # Factor to scale weights for rms value
    if idata!=1:
        f_ct=(ndt-j)/sw_ct    # Factor to scale weights for rms value

    """
    Reset current RMS variables
    """
    rms_cc0=0.
    rms_ct0=0.
    av_cc0=0.
    av_ct0=0.
    rms_cc=0.
    rms_ct=0.
    av_cc=0.
    av_ct=0.
    j=0
    """
    Calculate mean and RMS of data
    """
    for i in range(0,ndt):
        if idx[i]<=2:   # For CC data
            rms_cc0 += d[i]**2              # Unweighted
            av_cc0 += d[i]                  # Unweighted
            rms_cc += (f_cc*w[i]*d[i])**2   # Weighted and scaled
            av_cc += f_cc*w[i]*d[i]         # Weighted and scaled
            j = j+1
        else:   # For CT data
            rms_ct0 += d[i]**2              # Unweighted
            av_ct0 += d[i]                  # Unweighted
            rms_ct += (f_ct*w[i]*d[i])**2   # Weighted and scaled
            av_ct += f_ct*w[i]*d[i]         # Weighted and scaled
    if idata!=2:    # For cc data
        av_cc0 = av_cc0/j
        rms_cc0 = np.sqrt((rms_cc0-av_cc0**2/j)/(j-1))
        av_cc = av_cc/j
        rms_cc = np.sqrt((rms_cc-av_cc**2/j)/(j-1))
    if idata!=1:    # For ct data
        av_ct0 = av_ct0/(ndt-j)
        rms_ct0 = np.sqrt((rms_ct0-av_ct0**2/(ndt-j))/(ndt-j-1))
        av_ct = av_ct/(ndt-j)
        rms_ct = np.sqrt((rms_ct-av_ct**2/(ndt-j))/(ndt-j-1))

    """
    Calculate residual variance
    """
    if (abs(dum+999)<0.0001):
        dav = 0.
        dvar = 0.
        dav1 = 0.
        dvar1 = 0.
    # Update previous iter. variables if defined
    try:
        davold = dav
        dvarold = dvar
        dav1old = dav1 
        dvar1old = dvar1 
    except:
        davold = 0.
        dvarold = 0.
        dav1old = 0. 
        dvar1old = 0. 
    # Reset variables
    dav=0.   # unweighted
    dvar=0.  # unweighted
    dav1=0.  # unweighted
    dvar1=0. # unweighted
    sw=0.
    dav0=0.  # unweighted

    sw = np.sum(w)  # Sum weights
    f = ndt/sw      # Factor to scale weights for rms value
    # Sum data residuals
    dav = np.sum(d)         # unweighted sum
    dav0 = np.sum(w*d)      # weighted sum
    dav1 = np.sum(f*w*d)    # weighted and scaled sum
    # Average over no. of data
    dav = dav/ndt
    dav0 = dav0/ndt
    dav1 = dav1/ndt

    s = np.zeros(ndt)
    s1 = np.zeros(ndt)
    # Subtract out res. average
    s = d*1000 - dav*1000 # in msec
    s1 = w*d*1000 - dav0*1000 # weighted in msec
    # Sum and avg. the residual values
    ss = np.sum(s)
    ss1 = np.sum(s1)
    dvar = np.sum(s*s)
    dvar1 = np.sum(s1*s1)
    
    if reloctype!=3:
        if ndt>4*nev:
            dvar = (dvar - ss**2/ndt)/(ndt - 4*nev) # divide by number of degrees of freedom
            dvar1 = (dvar1 - ss1**2/ndt)/(ndt - 4*nev)
        else:
            dvar = (dvar/1) # divide by the number of degrees of freedom
            dvar1 = (dvar1/1)
            print('>>> Warning: ndt < 4*nev')
            log.write('>>> Warning: ndt < 4*nev \n')
    else:
        if ndt>3*nev:
            dvar = (dvar - ss**2/ndt)/(ndt - 3*nev) # divide by number of degrees of freedom
            dvar1 = (dvar1 - ss1**2/ndt)/(ndt - 3*nev)
        else:
            dvar = (dvar/1) # divide by the number of degrees of freedom
            dvar1 = (dvar1/1)
            print('>>> Warning: ndt < 4*nev')
            log.write('>>> Warning: ndt < 4*nev \n')

    """
    Record residual summary to log file
    """
    if (abs(dum+999)<0.0001): # original data
        log.write('\nResidual summary of initial data: \n')
        log.write(' absolute mean [s] = %7.4f \n' % dav)
        log.write(' weighted mean [s] = %7.4f \n' % dav1)
        log.write(' absolute variance [s] = %7.4f \n' % (dvar/1000))
        log.write(' weighted variance [s] = %7.4f \n' % (dvar1/1000))
        if idata==1 or idata==3:
            log.write(' absolute cc rms [s] = %7.4f \n' % rms_cc0)
            log.write(' weighted cc rms [s] (RMSCC) = %7.4f \n' % rms_cc)
        if idata==2 or idata==3:
            log.write(' absolute ct rms [s] = %7.4f \n' % rms_ct0)
            log.write(' weighted ct rms [s] (RMSCT) = %7.4f \n' % rms_ct)
    else:
        log.write('\nResidual summary: \n')
        if idata==1 or idata==3:
            log.write(' absolute cc rms [s] = %7.4f (%7.2f %%) \n' % (rms_cc0,100.*(rms_cc0-rms_cc0old)/rms_cc0old))
            log.write(' weighted cc rms [s] = %7.4f (%7.2f %%) \n' % (rms_cc,100.*(rms_cc-rms_ccold)/rms_ccold))
        if idata==2 or idata==3:
            log.write(' absolute ct rms [s] = %7.4f (%7.2f %%) \n' % (rms_ct0,100.*(rms_ct0-rms_ct0old)/rms_ct0old))
            log.write(' weighted ct rms [s] = %7.4f (%7.2f %%) \n' % (rms_ct,100.*(rms_ct-rms_ctold)/rms_ctold))

    dum = dvar1  # Update residual average

    """
    Return residual statistics
    """
    return [rms_cc, rms_ct, rms_cc0, rms_ct0, rms_ccold, 
            rms_ctold, rms_cc0old, rms_ct0old, dum]            


def dataprep(log,reloctype,fileout,fn_cc,fn_ct,fn_sta,fn_eve,idata,iphase,ncusp,icusp,
             maxdist,maxsep_ct,maxsep_cc):
    """
    This function reads in data files for hypoDD and organizes
    the data arrays.
    :::
    Parameters:
    log (file) - Log file
    fn_cc (str) - CC catalog file location
    fn_ct (str) - DT catalog file location
    fn_sta (str) - Station file location
    fn_eve (str) - Event file location
    idata (int) - Data type switch from input file
    iphase (int) - Phase type switch from input file
    ncusp (int) - Number of selected events
    icusp (int array) - IDs of selected events
    maxdist (float) - Max. event station separation
    maxsep_ct
    maxsep_cc
    :::
    Returns:
    retlist (list): list object containing all data parameters
    :::
    """

    """
    Begin write to log for data reading.
    """
    datet = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    log.write('\n\nReading data . . .  %s \n' % datet)

    """
    Read Events
    """
    [nev,ev_date,ev_time,ev_lat,ev_lon,ev_dep,ev_mag,ev_herr,ev_zerr,ev_res,
     ev_cusp] = readevent(log,fn_eve,ncusp,icusp)
    # Get center of event cluster for all events
    clat = np.sum(ev_lat)
    clon = np.sum(ev_lon)
    clat = clat/nev
    clon = clon/nev

    """
    Read Stations
    """
    [nsta,sta_lab,sta_lat,sta_lon] = readstat(log,fn_sta,maxdist,clat,clon,fileout)

    """
    Read in data files
    """
    [ndt,nccp,nccs,nctp,ncts,dt_sta1,dt_sta2,dt_dt,dt_qual,dt_offse,dt_offss,dt_c1,
     dt_c2,dt_idx] = readdata(log,reloctype,idata,iphase,fn_cc,fn_ct,nev,ev_cusp,
                              ev_lat,ev_lon,ev_dep,nsta,sta_lab,sta_lat,sta_lon,
                              maxsep_ct,maxsep_cc)
    ncc = nccp+nccs
    nct = nctp+ncts
    log.write('\n\nStarting number of data measurements:\n')
    log.write('NCCP: %i\n' % nccp)
    log.write('NCCS: %i\n' % nccs)
    log.write('Total Number of CC measurements: %i\n' % ncc)
    log.write('NCTP: %i\n' % nctp)
    log.write('NCTS: %i\n' % ncts)
    log.write('Total Number of CATALOG (CT) Measurements: %i\n' % nct)

    """
    Data Double Check
    ---
    Break if no data
    """
    log.write('Total Number of DTs total: %8i \n' % ndt)
    if ndt==0:
        raise Exception('Dtimes == 0. Stop triggered.')

    """
    Clean Data Arrays
    """
    dt_ic1 = np.zeros(ndt,dtype='int')
    dt_ic2 = np.zeros(ndt,dtype='int')
    dt_ista1 = np.zeros(ndt,dtype='int')
    dt_ista2 = np.zeros(ndt,dtype='int')
    [nev,ev_cusp,ev_date,ev_time,ev_lat,ev_lon,ev_dep,ev_mag,ev_herr,
     ev_zerr,ev_res,nsta,sta_lab,sta_lat,sta_lon,dt_ic1,dt_ic2,dt_ista1,
     dt_ista2] = cleanarrays(log,reloctype,nev,ev_cusp,ev_date,ev_time,ev_lat,ev_lon,
                             ev_dep,ev_mag,ev_herr,ev_zerr,ev_res,nsta,sta_lab,sta_lat,
                             sta_lon,ndt,dt_c1,dt_c2,dt_sta1,dt_sta2,dt_ic1,dt_ic2,
                             dt_ista1,dt_ista2)

    """
    List of Return Values
    """
    retlist = [ev_date,ev_time,ev_cusp,ev_lat,ev_lon,ev_dep,ev_mag,ev_herr,ev_zerr,
               ev_res,sta_lab,sta_lat,sta_lon,dt_sta1,dt_sta2,dt_dt,dt_qual,dt_c1,dt_c2,
               dt_idx,dt_ista1,dt_ista2,dt_ic1,dt_ic2,dt_offse,dt_offss,nev,nsta,ndt,
               ncc,nct,nccp,nccs,nctp,ncts]
    return retlist


def sigcoherency(log,reloctype,nct,ncc,ndt,idata,src_x,src_y,src_z,dt_ic1,dt_ic2,
                 dt_offse,dt_qual,dt_idx):
    """
    This function replaces the signal coherency loop in hypoDD.  Only
    needed if reloctype!=2.
    :::
    Parameters:
    log (file object) --- Log file
    reloctype (int) ---- Double-difference pairing type
    nct (int) --- No. of ct data
    ncc (int) --- No. of cc data
    ndt (int) --- No. of data
    idata (int) --- Data type switch
    src_x[nev] (float array) --- Src x locations
    src_y[nev] (float array) --- Src y locations
    src_z[nev] (float array) --- Src z locaitons
    dt_ic1[ndt] (int array) --- Event 1 IDs
    dt_ic2[ndt] (int array) --- Event 2 IDs
    dt_offse[ndt] (float array) --- Interevent offsets
    dt_qual[ndt] (float array) --- Data weights
    dt_idx[ndt] (int array) --- Data types
    :::
    Returns:
    ncc (int) --- No. of cc data
    nct (int) --- No. of ct data
    cohav (float) --- Correlation coherency
    picav (float) --- Pick coherency
    dt_offse[ndt] (float array) --- Updated interevent offsets
    :::
    """
    cohav = 0
    picav = 0
    j = nct
    k = ncc
    ncc = 0
    nct = 0

    for i in range(0,ndt):
        if reloctype!=2:
            dt_offse[i] = np.sqrt((src_x[int(dt_ic1[i])]-src_x[int(dt_ic2[i])])**2+
                                  (src_y[int(dt_ic1[i])]-src_y[int(dt_ic2[i])])**2+
                                  (src_z[int(dt_ic1[i])]-src_z[int(dt_ic2[i])])**2)
        if dt_idx[i]<=2:
            cohav = cohav + np.sqrt(dt_qual[i])
            ncc = ncc+1
        else:
            picav = picav + dt_qual[i]
            nct = nct+1
    
    log.write('\nMore: \n')
    if idata!=2:
        cohav = cohav/ncc
        log.write(' mean phase coherency = %5.3f \n' % cohav)
    if idata!=1:
        picav = picav/nct
        log.write(' mean pick quality = %5.3f \n' % picav)
    return ncc,nct,cohav,picav,dt_offse


def statres(log,nsta,ndt,idata,reloctype,sta_lab,dt_ista1,dt_ista2,dt_idx,dt_res):
    """
    Calculates residual information by station.
    :::
    Parameters:
    log (file object) --- Log file
    nsta (int) --- No. of stations
    ndt (int) --- No. of data
    idata (int) --- Data type switch
    reloctype (int) --- Double-difference pairing type
    sta_lab[nsta] (object array) --- Station codes
    dt_ista1[ndt] (int array) --- Station 1 indexes
    dt_ista2[ndt] (int array) --- Station 2 indexes
    dt_idx[ndt] (int array) --- Data type indexes
    dt_res[ndt] (float array) --- Data diff. time residuals
    :::
    sta_np[nsta] (int array) --- No. of P cc times per station
    sta_ns[nsta] (int array) --- No. of S cc times per station
    sta_nnp[nsta] (int array) --- No. of P ct times per station
    sta_nns[nsta] (int array) --- No. of S ct times per station
    sta_rmsc[nsta] (float array) --- Res. rms for cc times per station
    sta_rmsn[nsta] (float array) --- Res. rms for ct times per station
    tmpr1 (float) --- Largest cc station residual
    tmpr2 (float) --- Largest ct station residual
    :::
    """

    tmpr1 = float(0.)
    tmpr2 = float(0.)
    # Initalise arrays as zeros 
    sta_np = np.zeros(nsta,dtype='int')
    sta_ns = np.zeros(nsta,dtype='int')
    sta_nnp = np.zeros(nsta,dtype='int')
    sta_nns = np.zeros(nsta,dtype='int')
    sta_rmsc = np.zeros(nsta,dtype='float')
    sta_rmsn = np.zeros(nsta,dtype='float')

    for i in range(0,nsta):
        for j in range(0,ndt):
            if dt_ista1[j]==i:
                if dt_idx[j]<=2:
                    sta_rmsc[i] += dt_res[j]*dt_res[j] # Sum CC over each station
                elif dt_idx[j]>2:
                    sta_rmsn[i] += dt_res[j]*dt_res[j] # Sum CT over each station
                if dt_idx[j]==1:
                    sta_np[i] += 1  # No. of CC P Measurements
                elif dt_idx[j]==2:  
                    sta_ns[i] += 1   # No. of CC S Measurements
                elif dt_idx[j]==3:
                    sta_nnp[i] += 1   # No. of Cat P Measurements
                elif dt_idx[j]==4:
                    sta_nns[i] += 1   # No. of Cat S Measurements
            if reloctype==2 or reloctype==3:
                if dt_ista2[j]==i:
                    if dt_idx[j]<=2:
                        sta_rmsc[i] += dt_res[j]*dt_res[j] # Sum CC over each station
                    elif dt_idx[j]>2:
                        sta_rmsn[i] += dt_res[j]*dt_res[j] # Sum CT over each station
                    if dt_idx[j]==1:
                        sta_np[i] += 1  # No. of CC P Measurements
                    elif dt_idx[j]==2:  
                        sta_ns[i] += 1   # No. of CC S Measurements
                    elif dt_idx[j]==3:
                        sta_nnp[i] += 1   # No. of Cat P Measurements
                    elif dt_idx[j]==4:
                        sta_nns[i] += 1   # No. of Cat S Measurements
                    continue
        if (sta_np[i]+sta_ns[i])>0:
            sta_rmsc[i] = np.sqrt(sta_rmsc[i]/(sta_np[i]+sta_ns[i]))
        if (sta_nnp[i]+sta_nns[i])>0:
            sta_rmsn[i] = np.sqrt(sta_rmsn[i]/(sta_nnp[i]+sta_nns[i]))
        if sta_rmsc[i]>tmpr1:
            tmpr1 = sta_rmsc[i]
            k = i
        if sta_rmsn[i]>tmpr2:
            tmpr2 = sta_rmsn[i]
            l = i

    tmpr1 = tmpr1*1000.
    tmpr2 = tmpr2*1000.
    if idata==1 or idata==3:
        log.write(' station with largest cc rms: %s = %7.0f ms (RMSST) \n' % (sta_lab[k],tmpr1))
    if idata==2 or idata==3:
        log.write(' station with largest ct rms: %s = %7.0f ms (RMMST) \n' % (sta_lab[l],tmpr2))
    return sta_np,sta_ns,sta_nnp,sta_nns,sta_rmsc,sta_rmsn,tmpr1,tmpr2


def eventstats(nev,ndt,dt_ic1,dt_ic2,dt_res,dt_idx):
    """
    Res statistics for each event.
    :::
    Parameters:
    nev (int) --- No. of events
    ndt (int) --- No. of data
    dt_ic1[ndt] (int array) --- Data event 1 indexes
    dt_ic2[ndt] (int array) --- Data event 2 indexes
    dt_res[ndt] (float array) --- Data residuals
    dt_idx[ndt] (int array) --- Data type indexes
    :::
    Returns:
    src_np[nev] (int array) --- No. of CC P data per event
    src_ns[nev] (int array) --- No. of CC S data per event
    src_nnp[nev] (int array) --- No. of CT P data per event
    src_nns[nev] (int array) --- No. of CT S data per event
    src_rmsc[nev] (float array) --- RMS res. for cc data per event
    src_rmsn[nev] (float array) --- RMS res. for ct data per event
    :::
    """

    src_np = np.zeros(nev)
    src_ns = np.zeros(nev)
    src_nnp = np.zeros(nev)
    src_nns = np.zeros(nev)
    src_rmsc = np.zeros(nev)
    src_rmsn = np.zeros(nev)
    for i in range(0,ndt):
        if dt_idx[i]==1:
            src_np[int(dt_ic1[i])] = src_np[int(dt_ic1[i])]+1
            src_np[int(dt_ic2[i])] = src_np[int(dt_ic2[i])]+1
        if dt_idx[i]==2:
            src_ns[int(dt_ic1[i])] = src_ns[int(dt_ic1[i])]+1
            src_ns[int(dt_ic2[i])] = src_ns[int(dt_ic2[i])]+1
        if dt_idx[i]<=2:
            src_rmsc[int(dt_ic1[i])] = src_rmsc[int(dt_ic1[i])] + dt_res[i]**2
            src_rmsc[int(dt_ic2[i])] = src_rmsc[int(dt_ic2[i])] + dt_res[i]**2
        if dt_idx[i]==3:
            src_nnp[int(dt_ic1[i])] = src_nnp[int(dt_ic1[i])]+1
            src_nnp[int(dt_ic2[i])] = src_nnp[int(dt_ic2[i])]+1
        if dt_idx[i]==4:
            src_nns[int(dt_ic1[i])] = src_nns[int(dt_ic1[i])]+1
            src_nns[int(dt_ic2[i])] = src_nns[int(dt_ic2[i])]+1
        if dt_idx[i]>=3:
            src_rmsn[int(dt_ic1[i])] = src_rmsn[int(dt_ic1[i])] + dt_res[i]**2
            src_rmsn[int(dt_ic2[i])] = src_rmsn[int(dt_ic2[i])] + dt_res[i]**2
    for i in range(0,nev):
        if (src_np[i]+src_ns[i])>0:
            src_rmsc[i] = np.sqrt(src_rmsc[i]/(src_np[i]+src_ns[i]))
        else:
            src_rmsc[i] = -9
        if (src_nnp[i]+src_nns[i])>0:
            src_rmsn[i] = np.sqrt(src_rmsn[i]/(src_nnp[i]+src_nns[i]))
        else:
            src_rmsn[i] = -9

    return src_np,src_ns,src_nnp,src_nns,src_rmsc,src_rmsn










