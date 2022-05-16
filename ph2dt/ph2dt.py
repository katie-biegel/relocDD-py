#!/usr/bin/env python
import numpy as np
from datetime import datetime

# Import needed hypoDD functions
from utility.universal.misc import delaz
from ph2dt.ph2dt_subfunc import ph2dt_prep
from methodTypes.eventPair.ph2dt_ev import evpair_buildcatalog
from methodTypes.doublePair.ph2dt_doub import doubpair_buildcatalog
from methodTypes.stationPair.ph2dt_stat import statpair_buildcatalog

"""
The script contains the main ph2dt function.  
---
This includes the catalog building part of the function.  All
precursory data is read in the ph2dt_subfunc.ph2dt_prep function.
"""

def ph2dt(log,pinputs,makedata=0,datfol='datfiles',outfol='outputs',reloctype=1,fileout=0,
          idata=3,MAXEV=10000,datavals=[]):
    """
    Catalog Building Function
    :::
    Purpose:
        - Reads input file (defaults to ph2dt.inp) to get the input station and 
          phase file names, and the pairing parameters.
        - Reads and filters absolute travel-time data from network catalogs 
          to form travel-time data for: 
            - pairs of earthquakes (reloctype=1), pairs of stations (reloctype=2), or 
              double pairs of earthquakes and stations (reloctype=3). 
    :::   
    If fileout=0:
        Writes the files (input files to hypoDD):
            dt.ct       the travel times for the common stations for close 
                        earthquake pairs
            event.dat   earthquake list in dd format
            event.sel   selected earthquake list in dd format (should be the same 
                        as event.dat)
    Else:
        Returns data arrays for running hypoDD  
    :::
    Parameters:
        log (file obj)      --- Log file
        pinputs (list)      --- Inputs from ph2dt.inp file
        datfol (str)        --- Folder path location for data files
        outfol (str)        --- Folder path location for output files
        reloctype (int)     --- Data pairing type [default=1]
        fileout (int)       --- Switch write to file [default=0]
        idata (int)         --- Type of data in relocation
        MAXEV (int)         --- Max. no. of events [default=1000]
        ph2dtinput (list)   --- List of input variables (only used if fileout==1)
    :::
    Returns:
        catalog (list)      --- Either empty or contains data arrays
    :::
    """

    datet = str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    log.write('\n\n(ph2dt) starting ph2dt %s \n' % datet)
    print('\n\n(ph2dt) starting ph2dt %s' % datet)

    """
    Read in data
    """
    [statfile,phasefile,minwght,maxdist,maxoffsete,minoffsets,
     mnb,limobs_pair,minobs_pair,maxobs_pair] = pinputs

    if fileout==0:
        # Read in data files
        [nsta,s_lab,s_lat,s_lon,nev,lat,lon,depth,cuspid,dates,times,mag,herr,zerr,res,npha,
         nobs_ct,p_pha,p_sta,p_time,p_wghtr] = ph2dt_prep(log,pinputs,datfol,fileout,reloctype)
    elif fileout==1:
        # Set data parameters if read-in already
        [nsta,s_lab,s_lat,s_lon,nev,lat,lon,depth,cuspid,dates,times,mag,herr,zerr,res,
         npha,nobs_ct,p_pha,p_sta,p_time,p_wghtr] = datavals
    else:
       raise Exception('(ph2dt) ph2dt: Invalid fileout value.')

    """
    Build catalogs or data arrays
    ---
    This begins where data is sorted, paired, quality-checked, and then writen to file.
    Some additions for synthetic testing (indicated by comments.)
    """
    # Form dtimes
    if fileout==0:
        print('(ph2dt) Forming dtimes ...')
    log.write('(ph2dt) Forming dtimes ... \n')

    if reloctype == 1:
        """
        Event-pair catalog
        ---
        Returned catalog is either empty list or data arrays
        depending on fileout value
        """
        tempstart = datetime.now()
        catalog = evpair_buildcatalog(log,datfol,outfol,fileout,makedata,idata,maxdist,
                                      maxoffsete,mnb,limobs_pair,minobs_pair,maxobs_pair,
                                      nsta,s_lab,s_lat,s_lon,nev,lat,lon,depth,cuspid,
                                      dates,times,mag,herr,zerr,res,npha,nobs_ct,p_pha,
                                      p_sta,p_time,p_wghtr)
    elif reloctype == 2:
       """
       Station-pair catalog
       ---
       Returned catalog is either empty list or data arrays
       depending on fileout value
       """
       tempstart = datetime.now()
       catalog = statpair_buildcatalog(log,datfol,outfol,fileout,makedata,idata,maxdist,
                                       minoffsets,mnb,limobs_pair,minobs_pair,maxobs_pair,
                                       nsta,s_lab,s_lat,s_lon,nev,lat,lon,depth,cuspid,
                                       dates,times,mag,herr,zerr,res,npha,nobs_ct,p_pha,
                                       p_sta,p_time,p_wghtr)
    elif reloctype == 3:
       """
       Double-pair catalog 
       ---
       Returned catalog is either empty list or data arrays
       depending on fileout value
       """
       tempstart = datetime.now()
       catalog = doubpair_buildcatalog(log,datfol,outfol,fileout,makedata,idata,maxdist,
                                       maxoffsete,minoffsets,mnb,limobs_pair,minobs_pair,
                                       maxobs_pair,nsta,s_lab,s_lat,s_lon,nev,lat,lon,
                                       depth,cuspid,dates,times,mag,herr,zerr,res,npha,
                                       nobs_ct,p_pha,p_sta,p_time,p_wghtr)

    log.write('(ph2dt) Done with ph2dt.\n')
    if fileout==0:
        print('(ph2dt) Done with ph2dt. \n\n\n')
        return []
    elif fileout==1:
        [dt_sta1,dt_sta2,dt_dt,dt_qual,dt_c1,dt_c2,dt_idx,dt_ista1,dt_ista2,dt_ic1,dt_ic2,
         dt_offse,dt_offss,ndt,nccp,nccs,nctp,ncts] = catalog
        return [dates,times,cuspid,lat,lon,depth,mag,herr,zerr,res,s_lab,s_lat,s_lon,dt_sta1,
                dt_sta2,dt_dt,dt_qual,dt_c1,dt_c2,dt_idx,dt_ista1,dt_ista2,dt_ic1,dt_ic2,
                dt_offse,dt_offss,nev,nsta,ndt,nccp,nccs,nctp,ncts]

