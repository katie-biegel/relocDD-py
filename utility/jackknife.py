import numpy as np
import datetime
import sys
import matplotlib.pyplot as plt
import os


def removestation(log,reloctype,jks,
                  sta_lab,sta_lat,sta_lon,dt_sta1,
                  dt_sta2,dt_dt,dt_qual,dt_c1,dt_c2,dt_idx,
                  dt_ista1,dt_ista2,dt_ic1,dt_ic2,dt_offse,dt_offss,
                  nsta,ndt,ncc,nct,nccp,nccs,nctp,ncts):

    log.write('(Jackknife Solution): Removing station %s' % sta_lab[jks])
    print('(Jackknife Solution): Removing station %s' % sta_lab[jks])

    log.write('NDT before: %i \n' % ndt)
    log.write('NCT Before: %i and NCC Before: %i \n' % (nct,ncc))
    print('NDT before: %i \n' % ndt)
    print('NCT Before: %i and NCC Before: %i \n' % (nct,ncc))

    # Remove station from data arrays
    rmind = int(jks)
    if reloctype==1:
        dt_rmind = np.where(dt_ista1[:ndt]==jks)[0]
    else:
        dt_rmind = np.where(np.logical_or(dt_ista1[:ndt]==jks,dt_ista2[:ndt]==jks))[0]

    # Remove from all dt arrays
    dt_sta1 = dt_sta1[~np.isin(np.arange(ndt), dt_rmind)]
    if reloctype!=1:
        dt_sta2 = dt_sta2[~np.isin(np.arange(ndt), dt_rmind)]
    dt_dt = dt_dt[~np.isin(np.arange(ndt), dt_rmind)]
    dt_qual = dt_qual[~np.isin(np.arange(ndt), dt_rmind)]
    dt_c1 = dt_c1[~np.isin(np.arange(ndt), dt_rmind)]
    if reloctype!=2:
        dt_c2 = dt_c2[~np.isin(np.arange(ndt), dt_rmind)]
    dt_idx = dt_idx[~np.isin(np.arange(ndt), dt_rmind)]
    dt_ista1 = dt_ista1[~np.isin(np.arange(ndt), dt_rmind)]
    if reloctype!=1:
        dt_ista2 = dt_ista2[~np.isin(np.arange(ndt), dt_rmind)]
    dt_ic1 = dt_ic1[~np.isin(np.arange(ndt), dt_rmind)]
    if reloctype!=2:
        dt_ic2 = dt_ic2[~np.isin(np.arange(ndt), dt_rmind)]
    if reloctype!=2:
        dt_offse = dt_offse[~np.isin(np.arange(ndt), dt_rmind)]
    if reloctype!=1:
        dt_offss = dt_offss[~np.isin(np.arange(ndt), dt_rmind)]

    # Reset indices for stations
    for i in range(len(dt_sta1)):
        if reloctype==1:
            if dt_ista1[i]>jks:
                dt_ista1[i]=dt_ista1[i]-1
        else:
            if dt_ista1[i]>jks:
                dt_ista1[i]=dt_ista1[i]-1
            if dt_ista2[i]>jks:
                dt_ista2[i]=dt_ista2[i]-1

    # Remove from station arrays
    sta_lab = sta_lab[~np.isin(np.arange(nsta), rmind)]
    sta_lat = sta_lat[~np.isin(np.arange(nsta), rmind)]
    sta_lon = sta_lon[~np.isin(np.arange(nsta), rmind)]

    # Update counts
    nsta = len(sta_lab)
    ndt = len(dt_dt)
    nctp = len(dt_idx[dt_idx==3])    
    ncts = len(dt_idx[dt_idx==4])  
    nct = nctp+ncts
    nccp = len(dt_idx[dt_idx==1])    
    nccs = len(dt_idx[dt_idx==2])
    ncc = nccp+nccs  

    log.write('NDT after: %i \n' % ndt)
    log.write('NCT after: %i and NCC after: %i \n' % (nct,ncc))
    print('NDT after: %i \n' % ndt)
    print('NCT after: %i and NCC after: %i \n' % (nct,ncc))
    print('Station Removed.\n\n')
    log.write('Station Removed. \n\n')

    return [sta_lab[:nsta],sta_lat[:nsta],sta_lon[:nsta],
            dt_sta1[:ndt],dt_sta2[:ndt],dt_dt[:ndt],dt_qual[:ndt],
            dt_c1[:ndt],dt_c2[:ndt],dt_idx[:ndt],dt_ista1[:ndt],dt_ista2[:ndt],
            dt_ic1[:ndt],dt_ic2[:ndt],dt_offse[:ndt],dt_offss[:ndt],
            nsta,ndt,ncc,nct,nccp,nccs,nctp,ncts]


def removeevent(log,reloctype,jke,
                ev_date,ev_time,ev_cusp,ev_lat,ev_lon,ev_dep,ev_mag,
                ev_herr,ev_zerr,ev_res,dt_sta1,dt_sta2,
                dt_dt,dt_qual,dt_c1,dt_c2,dt_idx,
                dt_ista1,dt_ista2,dt_ic1,dt_ic2,dt_offse,dt_offss,
                nev,ndt,ncc,nct,nccp,nccs,nctp,ncts):

    log.write('(Jackknife Solution): Removing event %i \n' % ev_cusp[jke])
    print('(Jackknife Solution): Removing event %i' % ev_cusp[jke])

    log.write('NDT before: %i \n' % ndt)
    log.write('NCT Before: %i and NCC Before: %i \n' % (nct,ncc))
    print('NDT before: %i \n' % ndt)
    print('NCT Before: %i and NCC Before: %i \n' % (nct,ncc))

    # Remove station from data arrays
    rmind = jke*np.ones(1)

    dt_rmind = np.where(np.logical_or(dt_ic1==jke,dt_ic2==jke))[0]
    # Remove from all dt arrays
    dt_sta1 = dt_sta1[~np.isin(np.arange(ndt), dt_rmind)]
    if reloctype!=1:
        dt_sta2 = dt_sta2[~np.isin(np.arange(ndt), dt_rmind)]
    dt_dt = dt_dt[~np.isin(np.arange(ndt), dt_rmind)]
    dt_qual = dt_qual[~np.isin(np.arange(ndt), dt_rmind)]
    dt_c1 = dt_c1[~np.isin(np.arange(ndt), dt_rmind)]
    if reloctype!=2:
        dt_c2 = dt_c2[~np.isin(np.arange(ndt), dt_rmind)]
    dt_idx = dt_idx[~np.isin(np.arange(ndt), dt_rmind)]
    dt_ista1 = dt_ista1[~np.isin(np.arange(ndt), dt_rmind)]
    if reloctype!=1:
        dt_ista2 = dt_ista2[~np.isin(np.arange(ndt), dt_rmind)]
    dt_ic1 = dt_ic1[~np.isin(np.arange(ndt), dt_rmind)]
    if reloctype!=2:
        dt_ic2 = dt_ic2[~np.isin(np.arange(ndt), dt_rmind)]
    if reloctype!=2:
        dt_offse = dt_offse[~np.isin(np.arange(ndt), dt_rmind)]
    if reloctype!=1:    
        dt_offss = dt_offss[~np.isin(np.arange(ndt), dt_rmind)]

    # Reset indices for events
    for i in range(len(dt_sta1)):
        if reloctype==2:
            if dt_ic1[i]>jke:
                dt_ic1[i]=dt_ic1[i]-1
        else:
            if dt_ic1[i]>jke:
                dt_ic1[i]=dt_ic1[i]-1
            if dt_ic2[i]>jke:
                dt_ic2[i]=dt_ic2[i]-1

    # Remove from station arrays
    ev_date = ev_date[~np.isin(np.arange(nev), rmind)]
    ev_time = ev_time[~np.isin(np.arange(nev), rmind)]
    ev_cusp = ev_cusp[~np.isin(np.arange(nev), rmind)]
    ev_lat = ev_lat[~np.isin(np.arange(nev), rmind)]
    ev_lon = ev_lon[~np.isin(np.arange(nev), rmind)]
    ev_dep = ev_dep[~np.isin(np.arange(nev), rmind)]
    ev_mag = ev_mag[~np.isin(np.arange(nev), rmind)]
    ev_herr = ev_herr[~np.isin(np.arange(nev), rmind)]
    ev_zerr = ev_zerr[~np.isin(np.arange(nev), rmind)]
    ev_res = ev_res[~np.isin(np.arange(nev), rmind)]

    # Update counts
    ndt = len(dt_dt)
    nctp = len(dt_idx[dt_idx==3])    
    ncts = len(dt_idx[dt_idx==4])  
    nct = nctp+ncts
    nccp = len(dt_idx[dt_idx==1])    
    nccs = len(dt_idx[dt_idx==2])
    ncc = nccp+nccs  
    nev = len(ev_date)

    log.write('NDT after: %i \n' % ndt)
    log.write('NCT after: %i and NCC after: %i \n' % (nct,ncc))
    print('NDT after: %i \n' % ndt)
    print('NCT after: %i and NCC after: %i \n' % (nct,ncc))
    log.write('Event Removed. \n\n')

    return [ev_date[:nev],ev_time[:nev],ev_cusp[:nev],ev_lat[:nev],ev_lon[:nev],
            ev_dep[:nev],ev_mag[:nev],ev_herr[:nev],ev_zerr[:nev],ev_res[:nev],
            dt_sta1[:ndt],dt_sta2[:ndt],dt_dt[:ndt],dt_qual[:ndt],dt_c1[:ndt],
            dt_c2[:ndt],dt_idx[:ndt],dt_ista1[:ndt],dt_ista2[:ndt],
            dt_ic1[:ndt],dt_ic2[:ndt],dt_offse[:ndt],dt_offss[:ndt],
            nev,ndt,ncc,nct,nccp,nccs,nctp,ncts]
