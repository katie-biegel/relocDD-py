import numpy as np
import os
from datetime import datetime
from utility.universal.misc import delaz


def catalog_tofile(log,datfol,outfol,makedata,maxdist,minoffsets,mnb,limobs_pair,minobs_pair,
                   maxobs_pair,nsta,s_lab,s_lat,s_lon,nev,lat,lon,depth,cuspid,dates,times,
                   mag,herr,zerr,res,npha,nobs_ct,p_pha,p_sta,p_time,p_wghtr):
    """
    Overall hypoDD variables
    """
    PI = 3.141593
    KMPERDEG = 111.1949266

    """
    Declare counter variables (for log file)
    """
    n1 = int(0)      # Total No. pairs
    n3 = int(0)      # Total number of P pairs observed
    n4 = int(0)      # No. stations not listed in station.dat
    n5 = int(0)      # No. events that are too far from station-pair centroids to be considered
    n6 = int(0)      # Total number of S pairs observed
    n7 = int(0)      # No. of retained P pairs
    n8 = int(0)      # No. of retained S pairs
    nerr = int(0)    # No. outliers (too slow given velocity assumption)
    obsbreak=0

    """
    Open output files
    """
    dfile = os.path.join(datfol,'dts.ct')
    dts = open(dfile,'w')
    if makedata==1:  # KB - Added for synthetic cc file generation
        ccfile = os.path.join(datfol,'dts.cc')
        ccs = open(ccfile,'w')  

    """
    Declare data arrays
    ---
    Naming conventions:
        - a -> unsorted data pairs
        - b -> data pairs sorted by centroid-station separation
        - 1 -> indicates event 1 in pair
        - 2 -> indicates event 2 in pair
    """
    MAXOBS = int(2*nsta*nsta)                   # Rest maxobs since we have total # of phases
    a_lab1 = np.empty(MAXOBS,dtype='U7')        # Station label codes
    a_lab2 = np.empty(MAXOBS,dtype='U7')        # Station label codes
    a_pha = np.empty(MAXOBS,dtype='U1')         # 'P' or 'S' phase codes
    b_lab1 = np.empty(MAXOBS,dtype='U7')
    b_lab2 = np.empty(MAXOBS,dtype='U7')
    b_pha = np.empty(MAXOBS,dtype='U1')
    a_time1 = np.zeros(MAXOBS,dtype='float')    # Arrivaltime for event1 in pair
    a_time2 = np.zeros(MAXOBS,dtype='float')    # Arrivaltime for event2 in pair
    b_time1 = np.zeros(MAXOBS,dtype='float')
    b_time2 = np.zeros(MAXOBS,dtype='float')
    a_dist = np.zeros(MAXOBS,dtype='float')     # Event-pair centroid to station separation in km
    b_dist = np.zeros(MAXOBS,dtype='float')
    a_wtr = np.zeros(MAXOBS,dtype='float')      # Phase arrival weight from phase.dat
    b_wtr = np.zeros(MAXOBS,dtype='float')
    a_offss = np.zeros(MAXOBS,dtype='float') 
    b_offss = np.zeros(MAXOBS,dtype='float') 
    a_idx1 = np.zeros(MAXOBS,dtype='uint16')
    a_idx2 = np.zeros(MAXOBS,dtype='uint16')
    b_idx1 = np.zeros(MAXOBS,dtype='uint16')
    b_idx2 = np.zeros(MAXOBS,dtype='uint16')

    """
    Set up station-pairing strength test
    ---
    Take counts the number of data links between station
    pairs.
    Will be used later to count the number of weak/strong pairs.
    """
    take = np.zeros((nsta,nsta),dtype='uint16') 
    stp = 0

    """
    Declare loop variables
    """
    ipair = int(0)
    ipair_str = int(0)
    maxoff_str = int(0)
    tmp = int(0)
    aoffs = np.zeros(nsta,dtype='float')
    avoff = float(0.)
    avoff_str = float(0)
    
    """
    Loop over events
    ---
    One event at a time recorded to file
    Then station observations / station-pairs are indexed and
    recorded.
    """
    for i in range(nev):
        itmp = int(0)
        nobs = int(0)
        iobs = int(0)
        iimp = int(0) # obs needs to be included regardless of weight or dist
        count = int(0)

        dts.write('# %i \n' % (cuspid[i]))
        if makedata==1:
            ccs.write('# %i 0.0 \n' % (cuspid[i]))

        """
        Loop over the number of observations twice
        """
        for j in range(int(nobs_ct[i])-1):
            """
            Station 1 Index
            """
            try:
                sta1 = np.where(s_lab==p_sta[i,j])[0][0]
            except:
                log.write('Station not in station file: %s' % p_sta[i,j])
                continue
            alat = s_lat[sta1]
            alon = s_lon[sta1]
            # Station Offsets
            #[indx,aoffs] = stpair_offsets(aoffs,lati,loni,nsta,lat,lon,minoffsets)

            for l in range(j+1,int(nobs_ct[i])):
                """
                Station 2 index
                """
                try:
                    sta2 = np.where(s_lab==p_sta[i,l])[0][0]
                except:
                    log.write('Station not in station file: %s' % p_sta[i,l])
                    continue
                """
                Check for valid phase
                """
                if p_pha[i,j]==p_pha[i,l]:
                    if p_pha[i,j]=='P':
                        n3 += 1
                    if p_pha[i,j]=='S':
                        n6 += 1
                else:
                    continue

                alat = s_lat[sta1]
                alon = s_lon[sta1]
                blat = s_lat[sta2]
                blon = s_lon[sta2]
                clat = lat[i]
                clon = lon[i]
                """
                Maxdist check
                ----
                Calculate station-pair centroid to event separation
                Check if < maxdist
                """
                dlat = (alat+blat)/2
                dlon = (alon+blon)/2
                delt,dist,az = delaz(clat,clon,dlat,dlon)
                if dist>maxdist:
                    n5 += 1
                    continue
                if dist<minoffsets:
                    continue
                delta,dista,aza = delaz(alat,alon,blat,blon)
                if dista>=minoffsets: # and dista<150.:
                    #print('Good dist: ',dista,' sta1: ',sta1,' and sta2 ',sta2)

                    """
                    Calculate avg. weight
                    """
                    wtr = (np.abs(p_wghtr[i,j])+np.abs(p_wghtr[i,l]))/2

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
                        vel=4.
                    if p_pha[i,j]=='S':
                        vel=2.3
                    dist2 = np.sqrt(((alat-blat)*KMPERDEG)**2 +
                                    ((alon-blon)*np.cos(alon*PI/180.)*KMPERDEG)**2)

                    if np.abs(p_time[i,j]-p_time[i,l])>(dist2/vel+0.5):
                        log.write('Outlier: %s %s %9i %9.3f %9.3f %9.3f %9.3f \n' % 
                                  (s_lab[sta1],s_lab[sta2],cuspid[i],dist2,p_time[i,j],
                                   p_time[i,l],(p_time[i,j]-p_time[i,l])))
                        nerr += 1 
                        continue
                    
                    """
                    Save to arrays
                    """
                    a_lab1[itmp] = s_lab[sta1]
                    a_lab2[itmp] = s_lab[sta2]
                    a_pha[itmp] = p_pha[i,j]
                    a_time1[itmp] = p_time[i,j]
                    a_time2[itmp] = p_time[i,l]
                    a_dist[itmp] = dist
                    a_wtr[itmp] = wtr
                    a_offss[itmp] = dista
                    a_idx1[itmp] = sta1
                    a_idx2[itmp] = sta2
                    if p_wghtr[i,j]<0. or p_wghtr[i,l]<0.:
                        a_dist[itmp]=0.
                        iimp+=1
                    itmp+=1

        # Done with j loop (i.e., observations for event i)
        iobs = itmp
        if iobs>maxobs_pair: 
            itmp = min(maxobs_pair+iimp,iobs)
        if iobs>=minobs_pair:  # If the pair has more than the min obs per pair (otherwise skip)
            # Sort by distance betwen station and event-pair centriod
            iindx = np.argsort(a_dist[0:iobs])
            for kk in range(0,iobs):
                # Save data to sorted arrays by centroid - station separation
                # Closest stations first
                b_lab1[kk] = a_lab1[iindx[kk]]
                b_lab2[kk] = a_lab2[iindx[kk]]
                b_time1[kk] = a_time1[iindx[kk]]
                b_time2[kk] = a_time2[iindx[kk]]
                b_wtr[kk] = a_wtr[iindx[kk]]
                b_pha[kk] = a_pha[iindx[kk]]
                b_dist[kk] = a_dist[iindx[kk]]
                b_idx1[kk] = a_idx1[iindx[kk]]
                b_idx2[kk] = a_idx2[iindx[kk]]
                b_offss[kk] = a_offss[iindx[kk]]

            for i in range(itmp):
                dts.write('%s %s %7.6f %7.6f %6.4f %s \n' %
                          (b_lab1[i],b_lab2[i],b_time1[i],b_time2[i],
                           b_wtr[i],b_pha[i]))
                if makedata==1:
                    ccs.write('%s %s %7.6f 1.0 %s \n' % (b_lab1[i],b_lab2[i],
                                                         b_time1[i]-b_time2[i],
                                                         b_pha[i]))
        
                if b_pha[i]=='P':
                    n7+=1
                if b_pha[i]=='S':
                    n8+=1
                if take[b_idx1[i],b_idx2[i]]==0:
                    stp+=1
                    avoff += b_offss[i]
                take[b_idx1[i],b_idx2[i]]+=1
                if take[b_idx1[i],b_idx2[i]]>=limobs_pair:
                    n1+=1
                    avoff_str += b_offss[i]
        #else:
        #    continue

    """
    Here is where station pairs are classified as weak or strong pairs
    ----
    * Notice we switch i,k indexes to k,i (symmetric matrix)
    * Do this to avoid including duplicate station-pairs as we continue to 
    * loop over subsequent stations
    """
    #taketmp = take.flatten()
    #ipair = (taketmp>0).sum()
    #n1 = (taketmp>limobs_pair).sum()
    
    #avoff = float(0.)
    #avoff_str = float(0.)
    #for ind,value in enumerate(take):
    #    if value>0:
    #        avoff += stpair_offs[ind]
    #    if value>limobs_pair:
    #        avoff_str += stpair_offs[ind]

    avoff = float(avoff/stp)
    if n1>0:
        avoff_str = float(avoff_str/n1)

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
    print('> linked station pairs = ',stp)
    if n1>0:
        print('> strongly linked station pairs = ',n1)
    print('> average links per pair = ',int((n7+n8)/stp))
    print('> average offset (km) betw. linked stations = ',avoff)
    print('> average offset (km) betw. strongly linked stations = ',avoff_str)
    log.write('> P-phase pairs total = %i \n' % n3)
    log.write('> S-phase pairs total = %i \n' % n6)
    log.write('> outliers = %i (%3.2f percent) \n' % (nerr, nerr*100./(n3+n6)))
    log.write('> phases at stations not in station list = %i \n' % n4)
    log.write('> phases at distances larger than MAXDIST = %i \n' % n5)
    if n3>0:
        log.write('> P-phase pairs selected = %i (%3.2f percent) \n' % (n7,n7*100./n3))
    if n6>0:
        log.write('> S-phase pairs selected = %i (%3.2f percent) \n' % (n8,n8*100./n6))
    log.write('> linked event pairs %i \n' % stp)
    log.write('> average links per pair = %i' % int((n7+n8)/stp))
    log.write('> average offset (km) betw. linked stations = %f \n' % avoff)
    log.write('> average offset (km) betw. strongly linked stations = %f \n' % avoff_str)
    """
    Close all files
    """
    dts.close()
    if makedata==1:
        ccs.close()     # KB - ADDED FOR SYNTHETIC TESTING, xcorr file

    """
    Return empty list
    """
    return []

# #@profile(stream=open('mem_logs/catalog_tofile.mem','w+'))
# def catalog_tofile(log,makedata,
#                    maxdist,maxoffsets,mnb,limobs_pair,
#                    minobs_pair,maxobs_pair,
#                    nsta,s_lab,s_lat,s_lon,
#                    nev,lat,lon,depth,cuspid,dates,times,
#                    mag,herr,zerr,res,
#                    npha,nobs_ct,p_pha,p_sta,p_time,p_wghtr):
#     """
#     Overall hypoDD variables
#     """
#     PI = 3.141593
#     KMPERDEG = 111.1949266

#     """
#     Declare counter variables (for log file)
#     """
#     n1 = int(0)      # Total No. pairs
#     n3 = int(0)      # Total number of P pairs observed
#     n4 = int(0)      # No. stations not listed in station.dat
#     n5 = int(0)      # No. events that are too far from station-pair centroids to be considered
#     n6 = int(0)      # Total number of S pairs observed
#     n7 = int(0)      # No. of retained P pairs
#     n8 = int(0)      # No. of retained S pairs
#     nerr = int(0)    # No. outliers (too slow given velocity assumption)
#     obsbreak=0

#     dts = open('dt.ct','w')
#     if makedata==1:
#         ccs = open('dt.cc','w')         # KB - Added for synthetic cc file generation

#     """
#     Build stationpairs array
#     """
#     print('Starting statpairing')
#     [statpairs,stpair_offs,stp] = statpairing(nsta,mnb,maxoffsets,s_lat,s_lon,s_lab)

#     """
#     Set up station-pairing strength test
#     ---
#     Take counts the number of data links between station
#     pairs.
#     Will be used later to count the number of weak/strong pairs.
#     """
#     take = np.zeros((stp),dtype='int') 

#     """
#     Declare loop variables
#     """
#     ipair = int(1)
#     ipair_str = int(1)
#     maxoff_str = int(0)
#     tmp = int(0)
#     aoffs = np.zeros(nsta,dtype='float')
    
#     """
#     Loop over events
#     ---
#     One event at a time recorded to file
#     Then station observations / station-pairs are indexed and
#     recorded.
#     """
#     for i in range(nev):
#         nobs = int(0)
#         iobs = int(0)
#         iimp = int(0) # obs needs to be included regardless of weight or dist
#         count = int(0)

#         dts.write('# %i \n' % (cuspid[i]))
#         if makedata==1:
#             ccs.write('# %i 0.0 \n' % (cuspid[i]))

#         """
#         Loop over the number of observations twice
#         """
#         for j in range(int(nobs_ct[i])-1):
#             for l in range(j+1,int(nobs_ct[i])):
#                 stbreak=0
#                 """
#                 Set station indexes to impossible value
#                 """
#                 sta1 = -999
#                 sta2 = -999
#                 """
#                 Loop over to find our station pair at this event
#                 and confirm these stations are in station.dat
#                 """
#                 for ii in range(0,nsta):
#                     if p_sta[i,j]==s_lab[ii]:
#                         sta1=ii
#                     if p_sta[i,l]==s_lab[ii]:
#                         sta2=ii
#                 """
#                 Record bad stations
#                 """        
#                 if sta1==-999:
#                     log.write('Station not in station file: %s' % p_sta[i,j])
#                     n4 += 1
#                     stbreak=1
#                     break
#                 if sta2==-999:
#                     log.write('Station not in station file: %s' % p_sta[i,l])
#                     n4 += 1
#                     continue
#                 """
#                 Check for station pair in station-pair array
#                 """
#                 spcheck=-999
#                 for spind,pair in enumerate(statpairs):
#                     if sta1==pair[0] and sta2==pair[1]:
#                         spcheck=spind
#                         break
#                     if sta1==pair[1] and sta2==pair[0]:
#                         spcheck=spind
#                         break
#                 if spcheck==-999:
#                     continue
                
#                 """
#                 Check for valid phase
#                 """
#                 if p_pha[i,j]==p_pha[i,l]:
#                     if p_pha[i,j]=='P':
#                         n3 += 1
#                     if p_pha[i,j]=='S':
#                         n6 += 1
#                 else:
#                     continue

#                 """
#                 Save to arrays
#                 """
#                 alat = s_lat[sta1]
#                 alon = s_lon[sta1]
#                 blat = s_lat[sta2]
#                 blon = s_lon[sta2]
#                 clat = lat[i]
#                 clon = lon[i]
#                 """
#                 Maxdist check
#                 ----
#                 Calculate station-pair centroid to event separation
#                 Check if < maxdist
#                 """
#                 dlat = (alat+blat)/2
#                 dlon = (alon+blon)/2
#                 delt,dist,az = delaz(clat,clon,dlat,dlon)
#                 if dist>maxdist:
#                     n5 += 1
#                     continue

#                 """
#                 Calculate avg. weight
#                 """
#                 wtr = (np.abs(p_wghtr[i,j])+np.abs(p_wghtr[i,l]))/2

#                 """
#                 Remove outliers above the separation-delaytime line:
#                 ---
#                 Estimate unreasonably slow velocity, if the arrival separation is
#                 much greater than this value the data is considered an outlier.
#                 ---
#                 * We do record outliers in the log file for review but they
#                 * are not included in the ouptut files.
#                 """
#                 if p_pha[i,j]=='P':
#                     vel=4.
#                 if p_pha[i,j]=='S':
#                     vel=2.3
#                 dist2 = np.sqrt(((alat-blat)*KMPERDEG)**2 +
#                                 ((alon-blon)*np.cos(alon*PI/180.)*KMPERDEG)**2)

#                 if np.abs(p_time[i,j]-p_time[i,l])>(dist2/vel+0.5):
#                     log.write('Outlier: %s %s %9i %9.3f %9.3f %9.3f %9.3f \n' % 
#                               (s_lab[sta1],s_lab[sta2],cuspid[i],dist2,p_time[i,j],
#                                p_time[i,l],(p_time[i,j]-p_time[i,l])))
#                     nerr += 1
#                     continue
#                 """
#                 Write to file
#                 """
#                 if count<maxobs_pair:
#                     dts.write('%s %s %7.6f %7.6f %6.4f %s \n' %
#                               (s_lab[sta1],s_lab[sta2],p_time[i,j],p_time[i,l],
#                                wtr,p_pha[i,j]))
#                     if makedata==1:
#                         ccs.write('%s %s %7.6f 1.0 %s \n' % (s_lab[sta1],s_lab[sta2],
#                                                              p_time[i,j]-p_time[i,l], 
#                                                              p_pha[i,j]))
#                     if p_pha[i,j]=='P':
#                         n7+=1
#                     if p_pha[i,l]=='S':
#                         n8+=1
#                     take[spind]+=1
#                     if wtr>0:
#                         count+=1
#                 else:
#                     obsbreak==1
#             # Next outerloop if bad station
#             if stbreak==1:
#                 continue
#             if obsbreak==1:
#                 break
#         if obsbreak==1:
#             obsbreak=0
#             continue

#     """
#     Here is where station pairs are classified as weak or strong pairs
#     ----
#     * Notice we switch i,k indexes to k,i (symmetric matrix)
#     * Do this to avoid including duplicate station-pairs as we continue to 
#     * loop over subsequent stations
#     """
#     ipair = (take>0).sum()
#     n1 = (take>limobs_pair).sum()
#     avoff = float(0.)
#     avoff_str = float(0.)
#     for ind,value in enumerate(take):
#         if value>0:
#             avoff += stpair_offs[ind]
#         if value>limobs_pair:
#             avoff_str += stpair_offs[ind]
#     avoff = float(avoff/ipair)
#     avoff_str = float(avoff_str/n1)

#     """
#     Print final index counters to terminal and log
#     """
#     print('> P-phase pairs total = ',n3)
#     print('> S-phase pairs total = ',n6)
#     print('> outliers = ',nerr,' (',nerr*100./(n3+n6),'%)')
#     print('> phases at stations not in station list = ',n4)
#     print('> phases at distances larger than MAXDIST = ',n5)
#     if n3>0:
#         print('> P-phase pairs selected = ',n7,'(',n7*100./n3,'%)')
#     if n6>0:
#         print('> S-phase pairs selected = ',n8,'(',n8*100./n6,'%)')
#     print('> linked station pairs = ',ipair)
#     print('> strongly linked station pairs = ',n1)
#     print('> average links per pair = ',int((n7+n8)/ipair))
#     print('> average offset (km) betw. linked stations = ',avoff)
#     print('> average offset (km) betw. strongly linked stations = ',avoff_str)
#     log.write('> P-phase pairs total = %i \n' % n3)
#     log.write('> S-phase pairs total = %i \n' % n6)
#     log.write('> outliers = %i (%3.2f percent) \n' % (nerr, nerr*100./(n3+n6)))
#     log.write('> phases at stations not in station list = %i \n' % n4)
#     log.write('> phases at distances larger than MAXDIST = %i \n' % n5)
#     if n3>0:
#         log.write('> P-phase pairs selected = %i (%3.2f percent) \n' % (n7,n7*100./n3))
#     if n6>0:
#         log.write('> S-phase pairs selected = %i (%3.2f percent) \n' % (n8,n8*100./n6))
#     log.write('> linked event pairs %i \n' % ipair)
#     log.write('> average links per pair = %i' % int((n7+n8)/ipair))
#     log.write('> average offset (km) betw. linked stations = %f \n' % avoff)
#     log.write('> average offset (km) betw. strongly linked stations = %f \n' % avoff_str)
#     """
#     Close all files
#     """
#     dts.close()
#     if makedata==1:
#         ccs.close()     # KB - ADDED FOR SYNTHETIC TESTING, xcorr file

#     """
#     Return empty list
#     """
#     return []


#@profile(stream=open('mem_logs/catalog_toarrays.mem','w+'))
def catalog_toarrays(log,makedata,idata,
                     maxdist,maxoffsets,mnb,limobs_pair,minobs_pair,maxobs_pair,
                     nsta,s_lab,s_lat,s_lon,
                     nev,lat,lon,depth,cuspid,dates,times,mag,herr,zerr,res,
                     npha,nobs_ct,p_pha,p_sta,p_time,p_wghtr):

    """
    Overall hypoDD variables
    """
    PI = 3.141593
    KMPERDEG = 111.1949266

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
    Declare all data arrays for data passing
    ---
    ONLY IS fileout==1
    If we're skipping the file i/o between ph2dt and hypoDD
    then we need to initialise all the dt arrays that are returned from
    the hypoDD getdata function
    """
    
    # Dt arrays
    MAXOBS = int(npha*nsta) 
    dt_sta1 = np.empty(2*MAXOBS,dtype='U7')
    dt_sta2 = np.empty(2*MAXOBS,dtype='U7')
    dt_dt = np.zeros(2*MAXOBS,dtype='float')
    dt_qual = np.zeros(2*MAXOBS,dtype='float')
    dt_c1 = np.zeros(2*MAXOBS,dtype='int')
    dt_idx = np.zeros(2*MAXOBS,dtype='int')
    dt_ista1 = np.zeros(2*MAXOBS,dtype='int')
    dt_ista2 = np.zeros(2*MAXOBS,dtype='int')
    dt_ic1 = np.zeros(2*MAXOBS,dtype='int')
    dt_offss = np.zeros(2*MAXOBS,dtype='float')
    dt_idx_tmp = np.zeros(MAXOBS,dtype='int')
    ndt = int(0)
    nccp = int(0)
    nccs = int(0)
    nctp = int(0)
    ncts = int(0)

    """
    Build stationpairs array
    """
    [statpairs,stpair_offs,stp] = statpairing(nsta,mnb,maxoffsets,s_lat,s_lon,s_lab)

    """
    Set up station-pairing strength test
    ---
    Take counts the number of data links between station
    pairs.
    Will be used later to count the number of weak/strong pairs.
    """
    take = np.zeros((stp),dtype='int') 

    """
    Declare loop variables
    """
    ipair = int(1)
    ipair_str = int(1)
    maxoff_str = int(0)
    tmp = int(0)
    aoffs = np.zeros(nsta,dtype='float')
    
    """
    Loop over events
    ---
    One event at a time recorded to file
    Then station observations / station-pairs are indexed and
    recorded.
    """
    for i in range(nev):
        nobs = int(0)
        iobs = int(0)
        iimp = int(0) # obs needs to be included regardless of weight or dist

        """
        Loop over the number of observations twice
        """
        for j in range(int(nobs_ct[i])):
            for l in range(int(nobs_ct[i])):
                if j==l:        # Move to next observation if j==l
                    continue
            
                stbreak=0
            
                """
                Set station indexes to impossible value
                """
                sta1 = -999
                sta2 = -999

                """
                Loop over to find our station pair at this event
                and confirm these stations are in station.dat
                """
                for ii in range(0,nsta):
                    if p_sta[i,j]==s_lab[ii]:
                        sta1=ii
                    if p_sta[i,l]==s_lab[ii]:
                        sta2=ii

                """
                Record bad stations
                """        
                if sta1==-999:
                    log.write('Station not in station file: %s' % p_sta[i,j])
                    n4 += 1
                    stbreak=1
                    break
                if sta2==-999:
                    log.write('Station not in station file: %s' % p_sta[i,l])
                    n4 += 1
                    continue

                """
                Check for station pair in station-pair array
                """
                spcheck=-999
                for spind,pair in enumerate(statpairs):
                    if sta1==pair[0] and sta2==pair[1]:
                        spcheck=spind
                        break
                    if sta1==pair[1] and sta2==pair[0]:
                        spcheck=spind
                        break
                if spcheck==-999:
                    continue
                
                """
                Check for valid phase
                """
                if p_pha[i,j]==p_pha[i,l]:
                    if p_pha[i,j]=='P':
                        n3 += 1
                    if p_pha[i,j]=='S':
                        n6 += 1
                else:
                    continue

                """
                Save to arrays
                """
                alat = s_lat[sta1]
                alon = s_lon[sta1]
                blat = s_lat[sta2]
                blon = s_lon[sta2]
                clat = lat[i]
                clon = lon[i]

                """
                Maxdist check
                ----
                Calculate station-pair centroid to event separation
                Check if < maxdist
                """
                dlat = (alat+blat)/2
                dlon = (alon+blon)/2
                delt,dist,az = delaz(clat,clon,dlat,dlon)
                if dist>maxdist:
                    n5 += 1
                    continue

                """
                Calculate avg. weight
                """
                wtr = (np.abs(p_wghtr[i,j])+np.abs(p_wghtr[i,l]))/2

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
                    vel=4.
                if p_pha[i,j]=='S':
                    vel=2.3
                dist2 = np.sqrt(((alat-blat)*KMPERDEG)**2 +
                                ((alon-blon)*np.cos(alon*PI/180.)*KMPERDEG)**2)

                if np.abs(p_time[i,j]-p_time[i,l])>(dist2/vel+0.5):
                    log.write('Outlier: %s %s %9i %9.3f %9.3f %9.3f %9.3f \n' % 
                              (s_lab[sta1],s_lab[sta2],cuspid[i],dist2,p_time[i,j],
                               p_time[i,l],(p_time[i,j]-p_time[i,l])))
                    nerr += 1
                    continue

                """
                Save to data arrays
                ---
                Here is where we need to save to the same variables and order and return
                Those variables correctly to hypoDD for combined running
                """
                # Initialise Catalog Array
                dt_sta1[ndt] = s_lab[sta1]
                dt_sta2[ndt] = s_lab[sta2]
                dt_c1[ndt] = cuspid[i]

                dt_offss[ndt] = dist2
                dt_dt[ndt] = p_time[i,j]-p_time[i,l]
                dt_qual[ndt] = (p_wghtr[i,j]+p_wghtr[i,l])/2.

                if p_pha[i,j]=='P':
                    dt_idx[ndt] = 1
                    dt_idx_tmp[ndt] = 3
                    nccp += 1
                    nctp += 1
                elif p_pha[i,j]=='S':
                    dt_idx[ndt] = 2
                    dt_idx_tmp[ndt] = 4
                    ncts += 1
                    nccs += 1

                # Index arrays
                dt_ic1[ndt] = i
                dt_ista1[ndt] = sta1
                dt_ista2[ndt] = sta2                      
                ndt += 1

                if p_pha[i,j]=='P':
                    n7+=1
                if p_pha[i,l]=='S':
                    n8+=1
                take[spind]+=1

            # Next outerloop if bad station
            if stbreak==1:
                continue

    """
    Here is where station pairs are classified as weak or strong pairs
    ----
    * Notice we switch i,k indexes to k,i (symmetric matrix)
    * Do this to avoid including duplicate station-pairs as we continue to 
    * loop over subsequent stations
    """
    ipair = (take>0).sum()
    ipair_str = int(0)
    n1 = (take>limobs_pair).sum()

    avoff = float(0.)
    avoff_str = float(0.)
    for ind,value in enumerate(take):
        if value>0:
            avoff += stpair_offs[ind]
        if value>limobs_pair:
            avoff_str += stpair_offs[ind]
            ipair_str +=1

    avoff = float(avoff/ipair)
    avoff_str = float(avoff_str/n1)

    """
    Duplicate data arrays if synthetic
    """
    # Since the ccs are first, now duplicte eveything for cts
    if idata == 3:
        # CC Arrays
        dt_sta1[ndt:2*ndt] = dt_sta1[0:ndt]
        dt_sta2[ndt:2*ndt] = dt_sta2[0:ndt]
        dt_qual[ndt:2*ndt] = dt_qual[0:ndt]
        dt_c1[ndt:2*ndt] = dt_c1[0:ndt]
        dt_offss[ndt:2*ndt] = dt_offss[0:ndt]
        dt_dt[ndt:2*ndt] = dt_dt[0:ndt]
        dt_idx[ndt:2*ndt] = dt_idx_tmp[0:ndt]
        # Index arrays
        dt_ic1[ndt:2*ndt] = dt_ic1[0:ndt]
        dt_ista1[ndt:2*ndt] = dt_ista1[0:ndt]
        dt_ista2[ndt:2*ndt] = dt_ista2[0:ndt]
        # Counters
        ndt = 2*ndt
        nccp = nctp
        nccs = ncts 
    if idata == 2:
        # CC Arrays
        dt_idx = dt_idx_tmp[0:ndt]
        # Counters
        ndt = ndt
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
    dt_sta2 = dt_sta2[0:ndt]
    dt_qual = dt_qual[0:ndt]
    dt_c1 = dt_c1[0:ndt]
    dt_offss = dt_offss[0:ndt]
    dt_dt = dt_dt[0:ndt]
    dt_idx = dt_idx[0:ndt]
    dt_ic1 = dt_ic1[0:ndt]
    dt_ista1 = dt_ista1[0:ndt]
    dt_ista2 = dt_ista2[0:ndt]

    npair = ipair-1
    avoff = avoff/npair
    avoff_str = avoff_str/(ipair_str-1)

    # Other pairing arrays that are empty
    dt_c2 = np.zeros(ndt,dtype='int')
    dt_offse = np.zeros(ndt,dtype='float')
    dt_ic2 = np.zeros(ndt,dtype='int')

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
    print('> linked station pairs = ',ipair)
    print('> strongly linked station pairs = ',n1)
    print('> average links per pair = ',int((n7+n8)/ipair))
    print('> average offset (km) betw. linked stations = ',avoff)
    print('> average offset (km) betw. strongly linked stations = ',avoff_str)

    log.write('> P-phase pairs total = %i \n' % n3)
    log.write('> S-phase pairs total = %i \n' % n6)
    log.write('> outliers = %i (%3.2f percent) \n' % (nerr, nerr*100./(n3+n6)))
    log.write('> phases at stations not in station list = %i \n' % n4)
    log.write('> phases at distances larger than MAXDIST = %i \n' % n5)
    if n3>0:
        log.write('> P-phase pairs selected = %i (%3.2f percent) \n' % (n7,n7*100./n3))
    if n6>0:
        log.write('> S-phase pairs selected = %i (%3.2f percent) \n' % (n8,n8*100./n6))
    log.write('> linked event pairs %i \n' % ipair)
    log.write('> average links per pair = %i' % int((n7+n8)/ipair))
    log.write('> average offset (km) betw. linked stations = %f \n' % avoff)
    log.write('> average offset (km) betw. strongly linked stations = %f \n' % avoff_str)

    retlist = [dt_sta1,dt_sta2,dt_dt,dt_qual,dt_c1,dt_c2,
               dt_idx,dt_ista1,dt_ista2,dt_ic1,dt_ic2,dt_offse,dt_offss,
               ndt,nccp,nccs,nctp,ncts]
    return retlist


def statpair_buildcatalog(log,datfol,outfol,fileout,makedata,idata,maxdist,minoffsets,
                          mnb,limobs_pair,minobs_pair,maxobs_pair,nsta,s_lab,s_lat,
                          s_lon,nev,lat,lon,depth,cuspid,dates,times,mag,herr,zerr,res,
                          npha,nobs_ct,p_pha,p_sta,p_time,p_wghtr):

    log.write('\n\nBuilding station-pair catalog.\n\n')
    
    if fileout==0:
        print('\n\nBuilding station-pair catalog.')
        retlist = catalog_tofile(log,datfol,outfol,makedata,maxdist,minoffsets,mnb,
                                 limobs_pair,minobs_pair,maxobs_pair,nsta,s_lab,s_lat,
                                 s_lon,nev,lat,lon,depth,cuspid,dates,times,mag,herr,
                                 zerr,res,npha,nobs_ct,p_pha,p_sta,p_time,p_wghtr)
    """
    Else call catalog_toarrays
    """
    # if fileout==1:
    #     retlist = catalog_toarrays(log,makedata,idata,
    #                                maxdist,maxoffsets,mnb,limobs_pair,
    #                                minobs_pair,maxobs_pair,
    #                                nsta,s_lab,s_lat,s_lon,
    #                                nev,lat,lon,depth,cuspid,dates,times,
    #                                mag,herr,zerr,res,
    #                                npha,nobs_ct,p_pha,p_sta,p_time,p_wghtr)

    """
    Close Function
    """
    datet = str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    log.write('ph2dt done %s \n\n\n' % datet)
    print('ph2dt done %s' % datet)
    return retlist

