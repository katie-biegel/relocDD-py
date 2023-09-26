import obspy
import glob
import os
import numpy as np
from obspy.clients.fdsn import Client
from obspy.clients.iris import Client as tt

# Catalog
catfile = '/Users/katie/yukon/catalogs/kluanelake/master_catalog_FEB18.xml'
wavefol = '/Users/katie/yukon/waveforms/bash_pull_kluanelake/dir_dmc/mseeds'
waveformfol = '/Users/katie/yukon/waveforms/bash_pull_kluanelake/dir_dmc/waveforms'
os.makedirs(waveformfol,exist_ok=True)
client = Client('IRIS')

catalog = obspy.read_events(catfile)

minlat = 50.
maxlat = 70.
minlon = -160.
maxlon = -120.

# Write Event File
for event in catalog:
    # Pull Event Origin
    origin = event.origins[0]
    time = origin.time

    HR = origin.time.hour
    MN = origin.time.minute
    SEC = origin.time.second
    MS = origin.time.microsecond

    lat = origin.latitude
    lon = origin.longitude
    dep = origin.depth/1000.

    ID = '%02i%02i%02i%02i' % (HR,MN,SEC,int(MS*1e-4))

    # Wfile
    wfile = os.path.join(wavefol,'%s.miniseed' % ID)
    waveforms = obspy.read(wfile)


    for trace in waveforms:
        if len(trace.data)==0:
            continue
        net = trace.stats.network
        stat = trace.stats.station
        channel = trace.stats.channel
        start = trace.stats.starttime
        end = trace.stats.endtime

        station = client.get_stations(network=net,station=stat)[0][0]
        slat = station.latitude
        slon = station.longitude

        if slat<minlat or slat>maxlat:
            continue
        if slon<minlon or slon>maxlon:
            continue

        wfol2 = os.path.join(waveformfol,ID)
        os.makedirs(wfol2,exist_ok=True)

        if channel[-1]=='Z':
            parr = 0.
            try:
                result = tt_iris.traveltime(model='iasp91',phases=['P'],evdepth=dep,evloc=(lat,lon),staloc=(slat,slon))
                result = list(filter(None,str(result).split(' '))) 
                parr = float(result[26])

                stp = client.get_waveforms(network=net,station=stat,channel=channel,starttime=time+parr-15,endtime=time+parr+15)
                stp.write('%s/%s.%s.%s.%s.%s.MSEED' % (wfol2,ID,net,stat,channel,'P'),format='MSEED')
            except:
                continue
        else:
            sarr = 0.
            try:
                result = tt_iris.traveltime(model='iasp91',phases=['S'],evdepth=dep,evloc=(lat,lon),staloc=(slat,slon))
                result = list(filter(None,str(result).split(' '))) 
                sarr = float(result[26])

                stp = client.get_waveforms(network=net,station=stat,channel=channel,starttime=time+sarr-15,endtime=time+sarr+15)
                stp.write('%s/%s.%s.%s.%s.%s.MSEED' % (wfol2,ID,net,stat,channel,'S'),format='MSEED')
            except:
                continue

