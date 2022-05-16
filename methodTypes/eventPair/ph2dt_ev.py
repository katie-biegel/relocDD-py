# Builtin python packages
import os
import numpy as np
from datetime import datetime

# Other hypoDD functions
from ph2dt.ph2dt_subfunc import evpair_offsets
from utility.universal.misc import delaz

"""
Contains fucntions to build event-pair catalog
---
Includes functions:
    catalog_tofile
    catalog_toarray
    evpair_buildcatalog
---
"""

def catalog_tofile(log,datfol,outfol,makedata,maxdist,maxoffsete,mnb,limobs_pair,minobs_pair,
                   maxobs_pair,nsta,s_lab,s_lat,s_lon,nev,lat,lon,depth,cuspid,dates,times,
                   mag,herr,zerr,res,npha,nobs_ct,p_pha,p_sta,p_time,p_wghtr):

    """
    Declare counter variables (for log file)
    """
    n1 = int(0)      # Total No. pairs
    n2 = int(0)      # No. weakly linked event-pairs
    n3 = int(0)      # Total number of P pairs observed
    n4 = int(0)      # No. stations not listed in station.dat
    n5 = int(0)      # No. stations that are too far from event-pair centroids to be considered
    n6 = int(0)      # Total number of S pairs observed
    n7 = int(0)      # No. of retained P pairs
    n8 = int(0)      # No. of retained S pairs
    nerr = int(0)    # No. outliers (too slow given velocity assumption)

    """
    Open output files
    """
    dfile = os.path.join(datfol,'dte.ct')
    dts = open(dfile,'w')
    epfile = os.path.join(outfol,'evpairs.dat')
    evps = open(epfile,'w')
    if makedata==1:  # KB - Added for synthetic cc file generation
        ccfile = os.path.join(datfol,'dte.cc')
        ccs = open(ccfile,'w')  

    """
    Declare data arrays
    ---
    Naming conventions:
        - a -> unsorted data pairs
        - b -> data pairs sorted by centroid-station separation
        - 1 -> indicates event 1 in pair
        - 2 -> indicates event 2 in pair
    ---
    """
    MAXOBS = int(2*npha*nev)                # Rest maxobs since we have total # of phases
    a_lab = np.empty(MAXOBS,dtype='U7')     # Station label codes
    a_pha = np.empty(MAXOBS,dtype='U1')     # 'P' or 'S' phase codes
    b_lab = np.empty(MAXOBS,dtype='U7')
    b_pha = np.empty(MAXOBS,dtype='U1')
    a_time1 = np.zeros(MAXOBS,dtype='float')    # Arrivaltime for event1 in pair
    a_time2 = np.zeros(MAXOBS,dtype='float')    # Arrivaltime for event2 in pair
    b_time1 = np.zeros(MAXOBS,dtype='float')
    b_time2 = np.zeros(MAXOBS,dtype='float')
    a_dist = np.zeros(MAXOBS,dtype='float')     # Event-pair centroid to station separation in km
    b_dist = np.zeros(MAXOBS,dtype='float')
    a_wtr = np.zeros(MAXOBS,dtype='float')      # Phase arrival weight from phase.dat
    b_wtr = np.zeros(MAXOBS,dtype='float')

    """
    Set up event-pairing
    ---
    Take indicates how strongly event pair is linked.
        - 0 -> strong neighbour
        - 9 -> weak link
        - 1 -> Default (no linking)
    This is a lower half-matrix used to get rid of duplicate pairs.
    ---
    """
    take = np.ones((nev,nev),dtype='int')
    ipair = int(0)
    ipair_str = int(0)
    avoff = float(0)
    avoff_str = float(0)
    maxoff_str = float(0)
    tmp = int(0)
    aoffs = np.zeros(nev,dtype='float')
    stats_used = np.zeros(nsta,dtype='uint8')

    """
    Loop over the first event in the pair
    ---
    Loop over events (i indicates event1 in the event pair):
    ---
    """
    for i in range(0,nev):
        # Find inter-event separations for all other events
        indx,aoffs = evpair_offsets(aoffs,lat[i],lon[i],depth[i],nev,lat,lon,depth)
        # Reset neighbor indexes and obs. counter
        inb = int(0)
        nobs = int(0)

        """
        Now loop over the second event in the pair
        ---
        Looping over neighbours according to distance
        m does not indicate event but the event index in indx array
        """
        for m in range(nev-1):   # same event is last so excluded
            if inb>=mnb:         # If event i already has max # neighbors then break to next i
                break            # next event i

            """
            Declare event k (2nd event in pair)
            ---
            Pick event 2 by increasing distance in indx array
            """
            k = indx[m] # nearest event first

            # If this pair was already picked then skip to next k.
            # (If take is already set as weak or strong)
            if take[k,i] == 0: # already selected as a strong neighbor
                inb += 1
                continue
            elif take[k,i] == 9: # weak neighbor
                continue

            n1 += 1
            # Check max interevent offset
            # Since indx array sorted by distance, break to next i if not passed.
            if aoffs[indx[m]]>maxoffsete:
                break
            # Search for common stations/phases:
            iobs = int(0)
            iobsP = int(0)
            iimp = int(0) # obs needs to be included regardless of weight or dist

            """
            Loop over observations at event i
            ---
            j indicates observation index for event 1 (i)
            """
            for j in range(nobs_ct[i]):
                """
                Loop over observations at event k
                ---
                l indicates observation index for event 2 (k)
                """
                for l in range(nobs_ct[k]):

                    # Both observations must have same station and same phase (P or S)
                    if p_sta[i,j]==p_sta[k,l] and p_pha[i,j]==p_pha[k,l]: 
                        if p_pha[i,j]=='P':
                            n3+=1
                        if p_pha[i,j]=='S':
                            n6+=1

                        # Check that station file contains station label:
                        # p_* is data from phase.dat and s_* is from station.dat
                        try:
                            ii = np.argwhere(s_lab==p_sta[i,j])[0][0]
                        except:
                            n4+=1
                            break
                        alat = s_lat[ii]
                        alon = s_lon[ii]
                        blat = (lat[i] + lat[k])/2.
                        blon = (lon[i] + lon[k])/2.
                        delt,dist,az = delaz(alat,alon,blat,blon)

                        # If station - centroid separation too high indicate break
                        if dist>maxdist:
                            n5 += 1
                            break

                        # From ii loop above: We break loop when station is a match
                        # Now we found station index:
                        ista = ii
                        # Get avg. weight of 2 phase values
                        wtr = (np.abs(p_wghtr[i,j]) + np.abs(p_wghtr[k,l]))/2.

                        """
                        Remove outliers above the separation-delaytime line:
                        ---
                        Estimate unreasonably slow velocity, if the arrival separation is
                        much greater than this value the data is considered an outlier.
                        ---
                        * We do record outliers in the log file for review but they
                        * are not included in the ouptut files.
                        """
                        if p_pha[i,j]=='P':
                            vel = 4.
                        if p_pha[i,j]=='S':
                            vel = 2.3
                        if np.abs(p_time[i,j]-p_time[k,l])>(aoffs[indx[m]]/vel+0.5):
                            log.write('Outlier: %s %08i %08i %3.5f %3.5f %3.5f %3.5f \n' %
                                      (p_sta[i,j],cuspid[i],cuspid[k],aoffs[indx[m]],p_time[i,j],p_time[k,l],p_time[i,j]-p_time[k,l]))
                            nerr += 1
                            break # Go to next j

                        # Save data to unsorted arrays
                        a_lab[iobs] = s_lab[ista]
                        if stats_used[ista]==0:
                            stats_used[ista]+=1
                        a_time1[iobs] = p_time[i,j]
                        a_time2[iobs] = p_time[k,l]
                        a_wtr[iobs] = wtr
                        a_dist[iobs] = dist
                        a_pha[iobs] = p_pha[i,j]
                        # If weights are set negative in input, must include these picks.
                        # Convention set in hypoDD v1.3
                        if p_wghtr[i,j]<0. or p_wghtr[k,l]<0:
                            a_dist[iobs] = 0. # Set dist to 0 so it's selected first
                            iimp += 1

                        iobs += 1
                        if p_pha[i,j]=='P':
                            iobsP += 1
                        nobs += 1

                        break # Go to next j

            # Done with j loop (i.e., observations for event i)
            itmp = iobs
            if iobs>maxobs_pair: # itmp is the max obs limit (set to either maxobs + neg. weighted data
                                 # or all data for an event-pair - whichever is smaller)
                itmp = min(maxobs_pair+iimp,iobs)
            if iobs>=minobs_pair:  # If the pair has more than the min obs per pair (otherwise skip)
                if iobs>1:
                    # Sort by distance betwen station and event-pair centriod
                    iindx = np.argsort(a_dist[:iobs])
                    for kk in range(iobs):
                        # Save data to sorted arrays by centroid - station separation
                        # Closest stations first
                        b_lab[kk] = a_lab[iindx[kk]]
                        b_time1[kk] = a_time1[iindx[kk]]
                        b_time2[kk] = a_time2[iindx[kk]]
                        b_wtr[kk] = a_wtr[iindx[kk]]
                        b_pha[kk] = a_pha[iindx[kk]]
                        b_dist[kk] = a_dist[iindx[kk]]
                else:
                    b_lab[0] = a_lab[0]
                    b_time1[0] = a_time1[0]
                    b_time2[0] = a_time2[0]
                    b_wtr[0] = a_wtr[0]
                    b_pha[0] = a_pha[0]
                    b_dist[0] = a_dist[0]

                """
                Write to file
                """
                dts.write('# %8i %8i \n' % (cuspid[i],cuspid[k]))
                evps.write('%8i %8i \n' % (cuspid[i],cuspid[k]))
                if makedata==1:
                    """
                    Also write cc file if there is a synthetic model
                    """
                    ccs.write('# %8i %8i 0.0 \n' % (cuspid[i],cuspid[k])) # KB - ADDED FOR SIMULATIONS, to generate synthetic xcorr file
                    for kk in range(itmp):
                        stindx = np.where(s_lab==b_lab[kk])[0][0]
                        if b_pha[kk] == 'P':
                            ccs.write('%s %7.6f 1.0 %s \n' % (b_lab[kk],(b_time1[kk]-b_time2[kk]),b_pha[kk]))
                        if b_pha[kk] == 'S':
                            ccs.write('%s %7.6f 1.0 %s \n' % (b_lab[kk],(b_time1[kk]-b_time2[kk]),b_pha[kk]))
                for kk in range(itmp):
                    dts.write('%s %7.6f %7.6f %6.4f %s \n' % (b_lab[kk],b_time1[kk],b_time2[kk],b_wtr[kk],b_pha[kk]))
                    if b_pha[kk] == 'P':
                        n7 += 1
                    if b_pha[kk] == 'S':
                        n8 += 1
                avoff += aoffs[k]
                ipair += 1

            """
            Here is where event pairs are classified as weak or strong pairs
            ----
            * Notice we switch i,k indexes to k,i (symmetric matrix)
            * Do this to avoid including duplicate event-pairs as we countinue to 
            * loop over subsequent events.
            """
            if iobs>=limobs_pair: # select as strong neighbor, meets min. share obs.
                take[i,k] = 0
                inb += 1
                ipair_str += 1
                avoff_str += aoffs[k]
                if aoffs[k] > maxoff_str:
                    maxoff_str = aoffs[k]
            else: # weak neighbor, still linked but does not meet min. number of shared obs.
                take[i,k] = 9
        # Index over number of neighbors
        if inb<mnb:
            n2+=1

    npair = ipair #-1
    avoff = avoff/npair
    avoff_str = avoff_str/(ipair_str) #-1)

    #np.savetxt('event_pairing.txt',take.astype('int'))

    print('Stations used: ',np.count_nonzero(stats_used))
    """
    Print final index counters to terminal and log
    """
    print('> P-phase pairs total = ',n3)
    print('> S-phase pairs total = ',n6)
    print('> outliers = ',nerr,' (',nerr*100./(n3+n6),'%)')
    print('> phases at stations not in station list = ',n4)
    print('> phases at distances larger than MAXDIST = ',n5)
    if n3>0:
        print('> P-phase pairs selected = ',n7,'(',n7*100./n3,'%)')
    if n6>0:
        print('> S-phase pairs selected = ',n8,'(',n8*100./n6,'%)')
    print('> weakly linked events = ',n2,'(',n2*100./nev,'%)')
    print('> linked event pairs = ',ipair)
    print('> average links per pair = ',int((n7+n8)/ipair))
    print('> average offset (km) betw. linked events = ',avoff)
    print('> average offset (km) betw. strongly linked events = ',avoff_str)
    print('> maximum offset (km) betw. strongly linked events = ',maxoff_str)
    log.write('> P-phase pairs total = %i \n' % n3)
    log.write('> S-phase pairs total = %i \n' % n6)
    log.write('> outliers = %i (%3.2f percent) \n' % (nerr, nerr*100./(n3+n6)))
    log.write('> phases at stations not in station list = %i \n' % n4)
    log.write('> phases at distances larger than MAXDIST = %i \n' % n5)
    if n3>0:
        log.write('> P-phase pairs selected = %i (%3.2f percent) \n' % (n7,n7*100./n3))
    if n6>0:
        log.write('> S-phase pairs selected = %i (%3.2f percent) \n' % (n8,n8*100./n6))
    log.write('> weakly linked events = %i (%3.2f percent) \n' % (n2,n2*100./nev))
    log.write('> linked event pairs %i \n' % ipair)
    log.write('> average offset (km) betw. linked events = %4.5f \n' % avoff)
    log.write('> average offset (km) betw. strongly linked events = %4.5f \n' % avoff_str)
    log.write('> maximum offset (km) betw. strongly linked events = %4.5f \n' % maxoff_str)
    print('Output files: dt.ct; event.dat; event.sel; ph2dt.log')
    print('ph2dt parameters were: ')
    print('(maxdist,maxsep,maxngh,minlnk,minobs,maxobs)')
    print(maxdist,maxoffsete,mnb,limobs_pair,minobs_pair,maxobs_pair)

    """
    Close all files
    """
    dts.close()
    evps.close()
    if makedata==1:
        ccs.close()     # KB - ADDED FOR SYNTHETIC TESTING, xcorr file

    """
    Return empty list
    """
    return []


def catalog_toarrays(log,makedata,idata,
                     maxdist,maxoffsete,mnb,limobs_pair,minobs_pair,maxobs_pair,
                     nsta,s_lab,s_lat,s_lon,
                     nev,lat,lon,depth,cuspid,dates,times,mag,herr,zerr,res,
                     npha,nobs_ct,p_pha,p_sta,p_time,p_wghtr):

    """
    Declare counter variables (for log file)
    """
    n1 = int(0)      # Total No. pairs
    n2 = int(0)      # No. weakly linked event-pairs
    n3 = int(0)      # Total number of P pairs observed
    n4 = int(0)      # No. stations not listed in station.dat
    n5 = int(0)      # No. stations that are too far from event-pair centroids to be considered
    n6 = int(0)      # Total number of S pairs observed
    n7 = int(0)      # No. of retained P pairs
    n8 = int(0)      # No. of retained S pairs
    nerr = int(0)    # No. outliers (too slow given velocity assumption)

    """
    Declare data arrays
    ---
    Naming conventions:
        - a -> unsorted data pairs
        - b -> data pairs sorted by centroid-station separation
        - 1 -> indicates event 1 in pair
        - 2 -> indicates event 2 in pair
    """
    MAXOBS = int(npha*nev)                    # Max. no. of observations
    a_lab = np.empty(MAXOBS,dtype='U7')     # Station label codes
    a_pha = np.empty(MAXOBS,dtype='U1')     # 'P' or 'S' phase codes
    b_lab = np.empty(MAXOBS,dtype='U7')
    b_pha = np.empty(MAXOBS,dtype='U1')
    a_time1 = np.zeros(MAXOBS,dtype='float')    # Arrivaltime for event1 in pair
    a_time2 = np.zeros(MAXOBS,dtype='float')    # Arrivaltime for event2 in pair
    b_time1 = np.zeros(MAXOBS,dtype='float')
    b_time2 = np.zeros(MAXOBS,dtype='float')
    a_dist = np.zeros(MAXOBS,dtype='float')     # Event-pair centroid to station separation in km
    b_dist = np.zeros(MAXOBS,dtype='float')
    a_wtr = np.zeros(MAXOBS,dtype='float')      # Phase arrival weight from phase.dat
    b_wtr = np.zeros(MAXOBS,dtype='float')

    """
    Declare all data arrays for data passing
    ---
    ONLY IS fileout==1
    If we're skipping the file i/o between ph2dt and hypoDD
    then we need to initialise all the dt arrays that are returned from
    the hypoDD getdata function
    """

    # Dt arrays
    dt_sta1 = np.empty(2*MAXOBS,dtype='U7')
    dt_dt = np.zeros(2*MAXOBS,dtype='float')
    dt_qual = np.zeros(2*MAXOBS,dtype='float')
    dt_c1 = np.zeros(2*MAXOBS,dtype='int')
    dt_c2 = np.zeros(2*MAXOBS,dtype='int')
    dt_idx = np.zeros(2*MAXOBS,dtype='int')
    dt_ista1 = np.zeros(2*MAXOBS,dtype='int')
    dt_ic1 = np.zeros(2*MAXOBS,dtype='int')
    dt_ic2 = np.zeros(2*MAXOBS,dtype='int')
    dt_offse = np.zeros(2*MAXOBS,dtype='float')
    dt_idx_tmp = np.zeros(MAXOBS,dtype='int')
    ndt = int(0)
    nccp = int(0)
    nccs = int(0)
    nctp = int(0)
    ncts = int(0)

    """
    Set up event-pairing
    ---
    Take indicates how strongly event pair is linked.
        - 0 -> strong neighbour
        - 9 -> weak link
        - 1 -> Default (no linking)
    This is a lower half-matrix used to get rid of duplicate pairs.
    """
    take = np.ones((nev,nev),dtype='int')
    ipair = int(1)
    ipair_str = int(1)
    avoff = int(0)
    avoff_str = int(0)
    maxoff_str = int(0)
    tmp = int(0)
    aoffs = np.zeros(nev,dtype='float')

    """
    Loop over the first event in the pair
    ---
    Loop over events (i indicates event1 in the event pair):
    """
    for i in range(0,nev):
        # Find inter-event separations for all other events
        indx,aoffs = evpair_offsets(aoffs,lat[i],lon[i],depth[i],nev,lat,lon,depth)

        # Reset neighbor indexes and obs. counter
        inb = int(0)
        nobs = int(0)

        """
        Now loop over the second event in the pair
        ---
        Looping over neighbours according to distance
        m does not indicate event but the event index in indx array
        """
        for m in range(0,nev-1):    # same event is last so excluded
            if inb>=mnb:     # If event i already has max # neighbors then break to next i
                break # next event i

            """
            Declare event k (2nd event in pair)
            ---
            Pick event 2 by increasing distance in indx array
            """
            k = indx[m] # nearest event first

            # If this pair was already picked then skip to next k.
            # (If take is already set as weak or strong)
            if take[k,i] == 0: # already selected as a strong neighbor
                inb += 1
                continue
            elif take[k,i] == 9: # weak neighbor
                continue

            n1 += 1

            # Check max interevent offset
            # Since indx array sorted by distance, break to next i if not passed.
            if aoffs[indx[m]]>maxoffsete:
                break

            # Search for common stations/phases:
            iobs = int(0)
            iobsP = int(0)
            iimp = int(0) # obs needs to be included regardless of weight or dist

            """
            Loop over observations at event i
            ---
            j indicates observation index for event 1 (i)
            """
            for j in range(0,nobs_ct[i]):
                n5break=int(0)

                """
                Loop over observations at event k
                ---
                l indicates observation index for event 2 (k)
                """
                for l in range(0,nobs_ct[k]):

                    # Both observations must have same station and same phase (P or S)
                    if p_sta[i,j]==p_sta[k,l] and p_pha[i,j]==p_pha[k,l]: 
                        if p_pha[i,j]=='P':
                            n3 += 1
                        if p_pha[i,j]=='S':
                            n6 += 1

                        # Check that station file contains station label:
                        # p_* is data from phase.dat and s_* is from station.dat
                        okbreak = int(0)
                        for ii in range(nsta):
                            ok = int(0)
                            if p_sta[i,j]==s_lab[ii]:
                                ok = 1
                                stindx = ii
                            # For selected observation pair: Compute distance between station and 
                            # event-pair centroid in km
                            if ok==1:
                                alat = s_lat[ii]
                                alon = s_lon[ii]
                                blat = (lat[i] + lat[k])/2.
                                blon = (lon[i] + lon[k])/2.
                                delt,dist,az = delaz(alat,alon,blat,blon)

                                # If station - centroid separation too high indicate break
                                if dist>maxdist:
                                    n5 += 1
                                    n5break = 1
                                    break
                                okbreak = 1
                                break

                        if n5break==1:  # If station - centroid separation too high break to next j
                            break
                        if okbreak==0:  # If station not in file, break to next j
                            n4 += 1
                            break

                        # From ii loop above: We break loop when station is a match
                        # Now we found station index:
                        ista = ii
                        # Get avg. weight of 2 phase values
                        wtr = (np.abs(p_wghtr[i,j]) + np.abs(p_wghtr[k,l]))/2.

                        """
                        Remove outliers above the separation-delaytime line:
                        ---
                        Estimate unreasonably slow velocity, if the arrival separation is
                        much greater than this value the data is considered an outlier.
                        ---
                        * We do record outliers in the log file for review but they
                        * are not included in the ouptut files.
                        """
                        if p_pha[i,j]=='P':
                            vel = 4.
                        if p_pha[i,j]=='S':
                            vel = 2.3
                        if np.abs(p_time[i,j]-p_time[k,l])>(aoffs[indx[m]]/vel+0.5):
                            log.write('Outlier: %s %08i %08i %3.5f %3.5f %3.5f %3.5f \n' %
                                      (p_sta[i,j],cuspid[i],cuspid[k],aoffs[indx[m]],p_time[i,j],p_time[k,l],p_time[i,j]-p_time[k,l]))
                            nerr += 1
                            break # Go to next j

                        if p_pha[i,j]=='P':
                            iobsP += 1
                        nobs += 1

                        # Save data to unsorted arrays
                        a_lab[iobs] = s_lab[ista]
                        a_time1[iobs] = p_time[i,j]
                        a_time2[iobs] = p_time[k,l]
                        a_wtr[iobs] = wtr
                        a_dist[iobs] = dist
                        a_pha[iobs] = p_pha[i,j]
                        # If weights are set negative in input, must include these picks.
                        # Convention set in hypoDD v1.3
                        if p_wghtr[i,j]<0. or p_wghtr[k,l]<0:
                            a_dist[iobs] = 0. # Set dist to 0 so it's selected first
                            iimp += 1
                        iobs += 1
                        break # Go to next j

            # Done with j loop (i.e., observations for event i)
            itmp = iobs
            if iobs>maxobs_pair:   # itmp is the max obs limit (set to either maxobs + neg. weighted data
                                        # or all data for an event-pair - whichever is smaller)
                itmp = min(maxobs_pair+iimp,iobs)
            if iobs>=minobs_pair:  # If the pair has more than the min obs per pair (otherwise skip)
                if iobs>1:
                    # Sort by distance betwen station and event-pair centriod
                    iindx = np.argsort(a_dist[0:iobs])
                    for kk in range(0,iobs):
                        # Save data to sorted arrays by centroid - station separation
                        # Closest stations first
                        b_lab[kk] = a_lab[iindx[kk]]
                        b_time1[kk] = a_time1[iindx[kk]]
                        b_time2[kk] = a_time2[iindx[kk]]
                        b_wtr[kk] = a_wtr[iindx[kk]]
                        b_pha[kk] = a_pha[iindx[kk]]
                        b_dist[kk] = a_dist[iindx[kk]]
                else:
                    b_lab[0] = a_lab[0]
                    b_time1[0] = a_time1[0]
                    b_time2[0] = a_time2[0]
                    b_wtr[0] = a_wtr[0]
                    b_pha[0] = a_pha[0]
                    b_dist[0] = a_dist[0]
                """
                Save to data arrays
                ---
                Here is where we need to save to the same variables and order and return
                Those variables correctly to hypoDD for combined running
                """
                for kk in range(0,itmp):    ################# KB TODO: Remove loop
                    # Initialise Catalog Array
                    dt_sta1[ndt] = b_lab[kk]
                    dt_dt[ndt] = (b_time1[kk] - b_time2[kk])
                    dt_qual[ndt] = b_wtr[kk]
                    dt_c1[ndt] = cuspid[i]
                    dt_c2[ndt] = cuspid[k]
                    dt_offse[ndt] = b_dist[kk]
                    if b_pha[kk]=='P':
                        dt_idx[ndt] = 1
                        dt_idx_tmp[ndt] = 3
                        nccp += 1
                        nctp += 1
                    elif b_pha[kk]=='S':
                        dt_idx[ndt] = 2
                        dt_idx_tmp[ndt] = 4
                        ncts += 1
                        nccs += 1
                    # Index arrays
                    dt_ic1[ndt] = i
                    dt_ic2[ndt] = k
                    try:
                        dt_ista1[ndt] = np.where(s_lab==b_lab[kk])[0][0]
                    except: 
                        raise Exception('FATAL ERROR STATION INDEXING. PH2DT.')
                    ndt += 1

                avoff = avoff + aoffs[k]
                ipair += 1

            """
            Here is where event pairs are classified as weak or strong pairs
            ----
            * Notice we switch i,k indexes to k,i (symmetric matrix)
            * Do this to avoid including duplicate event-pairs as we countinue to 
            * loop over subsequent events.
            """
            if iobs>=limobs_pair: # select as strong neighbor, meets min. share obs.
                take[i,k] = 0
                inb += 1
                ipair_str += 1
                avoff_str += aoffs[k]
                if aoffs[k] > maxoff_str:
                    maxoff_str = aoffs[k]
            else: # weak neighbor, still linked but does not meet min. number of shared obs.
                take[i,k] = 9
        
        # Index over number of neighbors
        if inb<mnb:
            n2 += 1

    """
    Duplicate data arrays if synthetic
    """
    # Since the ccs are first, now duplicte eveything for cts
    if idata == 3:
        # CC Arrays
        dt_sta1[ndt:2*ndt] = dt_sta1[0:ndt]
        dt_qual[ndt:2*ndt] = dt_qual[0:ndt]
        dt_c1[ndt:2*ndt] = dt_c1[0:ndt]
        dt_c2[ndt:2*ndt] = dt_c2[0:ndt]
        dt_offse[ndt:2*ndt] = dt_offse[0:ndt]
        dt_dt[ndt:2*ndt] = dt_dt[0:ndt]
        dt_idx[ndt:2*ndt] = dt_idx_tmp[0:ndt]
        # Index arrays
        dt_ic1[ndt:2*ndt] = dt_ic1[0:ndt]
        dt_ic2[ndt:2*ndt] = dt_ic2[0:ndt]
        dt_ista1[ndt:2*ndt] = dt_ista1[0:ndt]
        # Counters
        ndt = 2*ndt
        nccp = nctp
        nccs = ncts 
    if idata == 2:
        # CC Arrays
        dt_idx = dt_idx_tmp[0:ndt]
        # Counters
        ndt = nctp+ncts
        nccp = 0
        nccs = 0
    if idata ==1:
        # Counters
        ndt = nccp+nccs
        nctp = 0
        ncts = 0

    """
    Trim to size
    """
    dt_sta1 = dt_sta1[0:ndt]
    dt_qual = dt_qual[0:ndt]
    dt_c1 = dt_c1[0:ndt]
    dt_c2 = dt_c2[0:ndt]
    dt_offse = dt_offse[0:ndt]
    dt_dt = dt_dt[0:ndt]
    dt_idx = dt_idx[0:ndt]
    dt_ic1 = dt_ic1[0:ndt]
    dt_ic2 = dt_ic2[0:ndt]
    dt_ista1 = dt_ista1[0:ndt]

    npair = ipair-1
    avoff = avoff/npair
    avoff_str = avoff_str/(ipair_str-1)

    # Other pairing arrays that are empty
    dt_sta2 = np.empty(ndt,dtype='U7')
    dt_offss = np.zeros(ndt,dtype='float')
    dt_ista2 = np.zeros(ndt,dtype='int')

    """
    Write to log and return data arrays
    """
    log.write('> P-phase pairs total = %i \n' % n3)
    log.write('> S-phase pairs total = %i \n' % n6)
    log.write('> outliers = %i (%3.2f percent) \n' % (nerr, nerr*100./(n3+n6)))
    log.write('> phases at stations not in station list = %i \n' % n4)
    log.write('> phases at distances larger than MAXDIST = %i \n' % n5)
    if n3>0:
        log.write('> P-phase pairs selected = %i (%3.2f percent) \n' % (nctp,nctp*100./n3))
    if n6>0:
        log.write('> S-phase pairs selected = %i (%3.2f percent) \n' % (ncts,ncts*100./n6))
    log.write('> weakly linked events = %i (%3.2f percent) \n' % (n2,n2*100./nev))
    log.write('> linked event pairs %i \n' % ipair)
    log.write('> average offset (km) betw. linked events = %4.5f \n' % avoff)
    log.write('> average offset (km) betw. strongly linked events = %4.5f \n' % avoff_str)
    log.write('> maximum offset (km) betw. strongly linked events = %4.5f \n' % maxoff_str)
    
    retlist = [dt_sta1,dt_sta2,dt_dt,dt_qual,dt_c1,dt_c2,
               dt_idx,dt_ista1,dt_ista2,dt_ic1,dt_ic2,dt_offse,dt_offss,
               ndt,nccp,nccs,nctp,ncts]
    return retlist


def evpair_buildcatalog(log,datfol,outfol,fileout,makedata,idata,maxdist,maxoffsete,mnb,
                        limobs_pair,minobs_pair,maxobs_pair,nsta,s_lab,s_lat,s_lon,nev,lat,
                        lon,depth,cuspid,dates,times,mag,herr,zerr,res,npha,nobs_ct,p_pha,
                        p_sta,p_time,p_wghtr):
    """
    This function calls the correct catalog building function for event-pairing
    and returns the data arrays (if fileout==1).
    ---
    Parameters:
    ---
    Returns:
    ---
    """
    log.write('\n\nBuilding event-pair catalog.\n\n')
    if fileout==0:
        print('\n\nBuilding event-pair catalog.')
    
    """
    If writing to file call catalog_tofile
    """
    retlist = []
    if fileout==0:
        retlist = catalog_tofile(log,datfol,outfol,makedata,maxdist,maxoffsete,mnb,limobs_pair,
                                 minobs_pair,maxobs_pair,nsta,s_lab,s_lat,s_lon,nev,lat,lon,depth,
                                 cuspid,dates,times,mag,herr,zerr,res,npha,nobs_ct,p_pha,p_sta,
                                 p_time,p_wghtr)
    # """
    # Else call catalog_toarrays
    # """
    # if fileout==1:
    #     tempstart = datetime.now()
    #     retlist = catalog_toarrays(log,makedata,idata,
    #                                maxdist,maxoffsete,mnb,limobs_pair,minobs_pair,maxobs_pair,
    #                                nsta,s_lab,s_lat,s_lon,
    #                                nev,lat,lon,depth,cuspid,dates,times,mag,herr,zerr,res,
    #                                npha,nobs_ct,p_pha,p_sta,p_time,p_wghtr)
    #     print('\n\nCatalog_toarrays completed.')
    #     print('Time to run catalog_toarrays (evpair_buildcatalog): ',datetime.now()-tempstart)
    #     log.write('\n\nCatalog_toarrays completed.')
    #     log.write('Time to run catalog_toarrays (evpair_buildcatalog): %i\n' % (datetime.now()-tempstart).seconds)

    """
    Close Function
    """
    datet = str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    log.write('ph2dt done %s \n\n\n' % datet)
    print('ph2dt done %s' % datet)
    return retlist
