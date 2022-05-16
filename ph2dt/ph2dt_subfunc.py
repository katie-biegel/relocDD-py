import numpy as np

from ph2dt.ph2dt_files import ph2dt_input,readphase
from utility.universal.misc import readstat, delaz
from datetime import datetime

"""
This script contains the subfunctions need for the ph2dt command.
---
This script contains the functions:
    evpair_offset:  Returns a sorted array of all possible evpair hypocentral
                    distances given one event lat and lon
    statpairing:    Pairs up stationspairs and returns array of station pair codes,
                    station pair offsets, and no. of station pairs
    ph2dt_prep:     Reads in the ph2dt input files and data needed for ph2dt
                    or
                    Sets up all the data for a synthetic model
---
"""

def evpair_offsets(aoffs,lati,loni,depi,nev,lat,lon,dep,KMPERDEG=111.1949266,PI=3.141593):
    """
    Calculates all possible eventpair hypocentral distances
    and returns a sorted array of indexes of increasing offsets
    :::
    Parameters:
    lati (float) --- Reference event latitude
    loni (float) --- Reference event longitude
    depi (float) --- Reference event depth
    nev (int) --- Number of events
    lat[nev] (float array) --- Latitudes for all events
    lon[nev] (float array) --- Longitudes for all events
    dep[nev] (float array) --- Depths for all events
    :::
    Returns:
    indx[nev] (int array) --- Ev indexes of sorted interevent offsets
    aoffs[nev] (float array) --- Interevent distances unsorted
    :::
    """

    """
    Loop over events and calculate hypocentral dist. from ref. event
    """
    for j in range(0,nev):
        dlat = lati - lat[j]
        dlon = loni - lon[j]
        x = dlat*KMPERDEG
        y = dlon*(np.cos(lati*PI/180.)*KMPERDEG)
        z = depi-dep[j]
        aoffs[j] = np.sqrt(x*x + y*y + z*z)
        if aoffs[j]<= 0.001:
            aoffs[j] = 99999    # Set same event to large dist (end of array)
    indx = np.argsort(aoffs)    # Sort events by distance to event i

    # Return indexes and offs
    return indx,aoffs


# #@profile(stream=open('mem_logs/evpair_offsets.mem','w+'))
# def stpair_offsets(aoffs,lati,loni,nsta,lat,lon,minoffsets,KMPERDEG=111.1949266,PI=3.141593):
#     """
#     Calculates all possible eventpair hypocentral distances
#     and returns a sorted array of indexes of increasing offsets
#     :::
#     Parameters:
#     lati (float) --- Reference event latitude
#     loni (float) --- Reference event longitude
#     depi (float) --- Reference event depth
#     nev (int) --- Number of events
#     lat[nev] (float array) --- Latitudes for all events
#     lon[nev] (float array) --- Longitudes for all events
#     dep[nev] (float array) --- Depths for all events
#     :::
#     Returns:
#     indx[nev] (int array) --- Ev indexes of sorted interevent offsets
#     aoffs[nev] (float array) --- Interevent distances unsorted
#     :::
#     """

#     """
#     Loop over events and calculate hypocentral dist. from ref. event
#     """
#     for j in range(nsta):
#         dlat = lati - lat[j]
#         dlon = loni - lon[j]
#         x = dlat*KMPERDEG
#         y = dlon*(np.cos(lati*PI/180.)*KMPERDEG)
#         aoffs[j] = np.sqrt(x*x + y*y)
#         if aoffs[j]<=minoffsets:
#             aoffs[j] = 99999    # Set same event to large dist (end of array)
#     indx = np.argsort(aoffs)    # Sort events by distance to event i

#     # Return indexes and offs
#     return indx,aoffs


# @profile(stream=open('mem_logs/statpairing.mem','w+'))
# def statpairing(nsta,mnb,minoffsets,s_lat,s_lon,s_lab,KMPERDEG=111.1949266,PI=3.141593):
#     print('Starting statpairing')
#     print('Nsta: ',nsta,' nsta**2: ',nsta*nsta)
#     statpairs = np.zeros((nsta*nsta,2),dtype='uint16')
#     statpairdist = np.zeros((nsta*nsta),dtype='float64')
#     stp = 0
#     """
#     First all possible station pairs and offsets
#     """
#     for sta1 in range(nsta-1):
#         lat1 = s_lat[sta1]
#         lon1 = s_lon[sta1]
#         for sta2 in range(sta1+1,nsta):
#             lat2 = s_lat[sta2]
#             lon2 = s_lon[sta2]
#             dels,dists,azs = delaz(lat1,lon1,lat2,lon2)
#             if dists>minoffsets:
#                 statpairs[stp,0] = sta1
#                 statpairs[stp,1] = sta2
#                 statpairdist[stp] = dists
#                 stp += 1
#     statpairs = statpairs[:stp,:].copy()
#     statpairdist = statpairdist[:stp].copy()
#     print('Stp 1: ',stp)

#     """
#     Now sort by offsets
#     """
#     offsort = np.argsort(statpairdist)
#     statpairs = statpairs[offsort,:]
#     statpairdist = statpairdist[offsort]
#     #mask = (statpairdist>minoffsets)
#     #statpairs = statpairs[mask,:]
#     #statpairdist = statpairdist[mask]
#     stpold = len(statpairs)

#     """
#     Loop over station pairs and only keep
#     5 neighbors per station
#     """
#     nongbr = np.zeros(nsta,dtype='int')
#     sneigh = np.zeros((nsta,nsta),dtype='int')
#     staprs = open('statpairs.txt','w')
#     stp = 0
#     for ipair in range(stpold):
#         i = statpairs[ipair,0]
#         j = statpairs[ipair,1]

#         # Check if already paired (this would be an error)
#         if sneigh[i,j]>0 or sneigh[j,i]>0:
#             print('Station Pairing Issue: Duplicate Pairs')
#             continue    # Already a station pair so don't save but also don't break
#         # Check if maxnbs already for both events
#         if nongbr[i]>=mnb:
#             continue
#         if nongbr[j]>=mnb:
#             continue

#         # If we pass all these checks we save the station-pair
#         # KB NOTE::::::: This preferences station-pairs with smaller separations
#         # Update all the neighbor counts
#         sneigh[i,j] += 1
#         sneigh[j,i] += 1
#         nongbr[i] += 1
#         nongbr[j] += 1
#         staprs.write('%s %s \n' % (s_lab[i],s_lab[j]))
#         statpairs[stp,0] = i
#         statpairs[stp,1] = j
#         statpairdist[stp] = statpairdist[ipair]
#         stp += 1

#     statpairs = statpairs[:stp,:]
#     statpairdist = statpairdist[:stp]
#     staprs.close()
#     return [statpairs,statpairdist,stp] # Return pair IDs, pair offsets, and no. of pairs

def ph2dt_prep(log,pinputs,datfol='datafiles',fileout=0,reloctype=1):
    """
    Read in data files for ph2dt_prep
    :::
    Parameters:
    log (file object)   --- Log file
    pinputs (list)      --- ph2dt input file location
    datfol (str)        --- Data file folder
    fileout (int)       --- File output switch
    reloctype (int)     --- Pairing type
    :::
    Returns:
    retlist (list)      --- Contains data required for ph2dt
    :::
    """

    """
    Read in ph2dt input file
    """
    [statfile,phasefile,minwght,maxdist,maxoffsete,maxoffsets,
     mnb,limobs_pair,minobs_pair,maxobs_pair] = pinputs

    """
    Read station file
    """
    [nsta,s_lab,s_lat,s_lon] = readstat(log,statfile,fileout=fileout) 

    """
    Phase/event data section
    """
    [nev,lat,lon,depth,cuspid,dates,times,mag,herr,zerr,res,npha,nobs_ct,p_pha,p_sta,
     p_time,p_wghtr] = readphase(log,phasefile,datfol,minwght,minobs_pair,nsta,fileout)

    """
    Return list of variables
    """
    retlist = [nsta,s_lab,s_lat,s_lon,nev,lat,lon,depth,cuspid,dates,times,
               mag,herr,zerr,res,npha,nobs_ct,p_pha,p_sta,p_time,p_wghtr]
    return retlist

