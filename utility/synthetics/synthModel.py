# Import General Python Packages
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import utm

# Import other hypoDD functions
from utility.universal.misc import readstat
from utility.universal.raytrace import partials
from utility.universal.geodetics import setorg,sdc2

"""
This script includes the functions necessary to make a synthetic
model from the synth.inp file.
---
The functions included are the subfunctions:
    synth_input:        Reads into memeory the synth.inp input file,
    write_orig:         Write original locations to file,
    define_events:      Define and return synthetic model events given synth.inp
                        specifications,
    write_phase:        Writes a phase.dat file for a synthetic model;
And the main synthetic model function:
    synthetic_generate:     Generate synthetic model data arrays from synth.inp
                            specifications.
---
"""

"""
Synthetic Model Subfunctions
"""
def synth_input(log,syninp='synth.inp'):
    """
    Read in synthetic model input file
    :::
    Parameters:
    log (file object) --- Log file
    syninp (str) --- Synthetic model input file location
    :::
    Returns:
    retlist (list) --- List of synthetic model input parameters
    :::
    """

    """
    Open and read file
    """
    sinputfile = open(syninp,'r')
    inputs = sinputfile.readlines()
    log.write('Reading in synthetic model input file.\n')

    l=int(0)
    for line in inputs:
        if line[0]=='*':                # Comment Lines
            continue
        if l==0:                        # Station.dat file location
            try:                            # For all variables try to read in
                line = line.split('\n')
                statfile = str(line[0])        # If read-in ok, set variable
                if len(statfile)==0:          # Else if value empty set default
                    statfile = str('station.dat')  
                    log.write('Default stationfile = station.dat.\n')
            except:                         # If read-in bad,
                statfile = str('station.dat')   # Also set to default
                log.write('Default stationfile = station.dat.\n')
        elif l==1:
            try:
                line = line.split('\n')
                nev = int(line[0])          # Number of events
            except:
                raise Exception('Syn.inp: Bad nev line.')
            else:
                if nev<=0:
                    raise Exception('Syn.inp: Nev must be positive and >0.')
        elif l==2: 
            try: 
                line = line.split(' ')
                line = list(filter(None,line))
                low = float(line[0])
                high = float(line[1])
                depth = np.random.uniform(low=low,high=high,size=nev)  # Event Depths
            except:
                raise Exception('Syn.inp: Bad event depth line.')
            else:
                if len(depth)!=nev:
                    raise Exception('Syn.inp: Incorrect number of event depths.')
        elif l==3:
            try:
                line = line.split('\n')
                stepsize = float(line[0])   # X/Y event step size
            except:
                raise Exception('Syn.inp: Bad stepsize line.')
            else:
                if stepsize<0:
                    raise Exception('Syn.inp: Stepsize must be positive.')  
        elif l==4:
            try:
                line = line.split('\n')
                latmid = float(line[0])     # Middle Latitude
            except:
                raise Exception('Syn.inp: Bad middle latitude line.')
        elif l==5:
            try:
                line = line.split('\n')
                lonmid = float(line[0])     # Middle Longitude
            except:
                raise Exception('Syn.inp: Bad middle longitude line.')
        elif l==6:
            try:
                line = line.split('\n')
                latstep = int(line[0])      # Latitude changes switch
                if latstep!=0 and latstep!=1:
                    latstep=0
                    log.write('Default latstep set to 0 (no latitude changes in synthetic model. \n')
            except:
                latstep=0
                log.write('Default latstep set to 0 (no latitude changes in synthetic model. \n')
        elif l==7:
            try:
                line = line.split('\n')
                lonstep = int(line[0])      # Longitude changes switch
                if lonstep!=0 and lonstep!=1:
                    lonstep=0
                    log.write('Default lonstep set to 0 (no longitude changes in synthetic model. \n')
            except:
                lonstep=0
                log.write('Default lonstep set to 0 (no longitude changes in synthetic model. \n')                
        elif l==8:
            try:
                line = line.split('\n')
                nl = int(line[0])       # Raytracing velocity model number of layers
            except:
                raise Exception('Syn.inp: Bad vel. model no. of layers line.')
            else:
                if nl<1:
                    raise Exception('Syn.inp: Invalid number of vel. model layers. Must be >=1.')
        elif l==9:
            try:
                line = line.split('\n')
                ratio = float(line[0])  # Raytracing vpvs ratio
            except:
                raise Exception('Syn.inp: Bad vel. model ratio line.')
            else:
                if ratio<1:
                    raise Exception('Syn.inp: Invalid vpvs ratio. Must be >1.')
        elif l==10:
            try:
                line = line.split(' ')
                line = list(filter(None,line))
                top = [float(ele) for ele in line[0:-1]]    # Raytracing layer depths
            except:
                raise Exception('Syn.inp: Bad vel. model top line.')
            else:
                if len(top)!=nl:
                    raise Exception('Syn.inp: Invalid no. of vel. layers. Top length must match nl.')
        elif l==11:
            try:
                line = line.split(' ')
                line = list(filter(None,line))
                v = [float(ele) for ele in line[0:-1]]      # Raytracing layer velocities
            except:
                raise Exception('Syn.inp: Bad vel. model velocity line.')
            else:
                if len(v)!=nl:
                    raise Exception('Syn.inp: Invalid no. of vel. layers. V length must match nl.')
        l+=1

    """
    Record variables to input file
    """
    log.write('\nSynth.inp file variables. \n\n')
    log.write('Statfile: %s \n' % statfile)
    log.write('Number of events: %i \n' % nev)
    if latstep==1:
        log.write('Latitude changes for synthetic event locations. \n')
    if lonstep==1:
        log.write('Longitude changes for synthetic event locations. \n')
    log.write('Raytracing velocity model: no. of layers = %i, VPVS ratio = %f \n' % (nl,ratio))

    """
    Return input parameters
    """
    retlist = [statfile,
               nev,depth,latmid,lonmid,
               stepsize,latstep,lonstep,
               nl,ratio,top,v]
    return retlist


def write_orig(nev,lat,lon,depth,dates,times,ids=None,filename=None):
    """
    This function writes the true event locations to file.
    :::
    PARAMETERS:
    nev (int) - number of events
    lat[nev] (float array) - numpy array containing event latitudes
    lon[nev] (float array) - numpy array containing event longitudes
    depth[nev] (float array) - numpy array containing event depths
    datetimes[nev] (object array) - numpy array containing event origin datetime.datetime objects
    :::
    Returns:
    None
    :::
    """

    """
    Write each true event lat/lon oriing to file
    ---
    The event format is:
    lat lon depth YYYYMODY HRMNSECS

    Where YYYYMODY is year, month, day in the format %4i%02i%02i
    And HRMNSECS is hour, minute, seconds and milliseconds in the format %02i%02i%04i
    The seconds portion of this string is assumed to have 2 decimal places if %2.2f
    with no decimal placeholder.
    """
    if filename:
        origLocfile = open(filename,'w') 
    else:
        origLocfile = open('origLatLon.loc','w') 
    for iev in range(nev):
        if len(ids)>0:
            origLocfile.write('%f %f %f %4i%02i%02i %02i%02i%04i \n' % (lat[iev],lon[iev],depth[iev],
                                                                        dates[iev].year,dates[iev].month,
                                                                        dates[iev].day,times[iev].hour,
                                                                        times[iev].minute,
                                                                        int((times[iev].second + times[iev].microsecond/1000000)*100)))
        else:
            origLocfile.write('%i %f %f %f %4i%02i%02i %02i%02i%04i \n' % (ids[iev],lat[iev],lon[iev],depth[iev],
                                                                           dates[iev].year,dates[iev].month,
                                                                           dates[iev].day,times[iev].hour,
                                                                           times[iev].minute,
                                                                           int((times[iev].second + times[iev].microsecond/1000000)*100)))
    origLocfile.close()

    """
    Temporary calculate cluster centroid coordinates for one
    cluster of all events
    ---
    HypoDD initiates the cluster centroid as the cartesian
    coordinate system origin.

    Here we do this to save the true event cartesian coordinates
    to file as these will be replaced by the hypoinverse catalog.
    """
    sdc0_lat = 0.
    sdc0_lon = 0.
    sdc0_dep = 0.
    sdc0_lat = np.sum(lat)
    sdc0_lon = np.sum(lon)
    sdc0_dep = np.sum(depth)
    sdc0_lat = sdc0_lat/nev
    sdc0_lon = sdc0_lon/nev
    sdc0_dep = sdc0_dep/nev
    setorg(sdc0_lat,sdc0_lon,0.0,0)

    """
    Write true cartesian locations to origXY.loc (needed for plotting)
    ---
    The format is:
    X Y Z (all in meters from cluster centroid)
    """
    if filename:
        origXYfile = open(filename+'.xy','w')
    else:
        origXYfile = open('origXY.loc','w')
    for i in range(0,nev):
        [evx,evy] = sdc2(lat[i],lon[i],-1)
        origXYfile.write('%f %f %f 0.0 \n' % (evx*1000,evy*1000,(depth[i]-sdc0_dep)*1000))
    origXYfile.close()

    return None


def define_events(log,nev,latmid,lonmid,stepsize,latstep,lonstep,depth,icusp):
    """
    Define and return synthetic model events given synth.inp file variables.
    :::
    Parameters:
    log (file object) --- Log file
    nev (int) --- No. of events
    latmid (float) --- Middle Latitude
    lonmid (float) --- Middle Longitude
    stepsize (float) --- Change in lat/lon step size betweene events
    latstep (float) --- Switch for changing latitudes
    lonstep (float) --- Switch for changing longitudes
    :::
    Returns:
    cuspid[nev] (int array) --- Event IDs
    lat[nev] (float array) --- Event latitudes
    lon[nev] (float array) --- Event longitudes
    dates[nev] (object array) --- Event dates (datetime.date objects)
    times[nev] (object array) --- Event times (datetime.time objects)
    mag[nev] (float array) --- Event magnitudes
    herr[nev] (float array) --- Event horizontal location errors
    verr[nev] (float array) --- Event vertical location errors
    res[nev] (float array) --- Event residuals
    :::
    """

    # readevent=bool(True)
    # if not readevent:
    # #     """
    #     First Define Event Info Arrays
    #     """
    #     # Event IDS - linspace array
    #     cuspid = np.arange(1,nev+1,dtype='int')

    #     # Define Lat/Lon based on the synth model parameters
    #     # Step sizes
    #     step = float(stepsize/nev)
    #     halfrange = step*(int(nev/2.))
    #     # Define range
    #     hr2 = halfrange + step
    #     latmin = latmid - halfrange
    #     latmax = latmid + hr2
    #     lonmin = lonmid - halfrange
    #     lonmax = lonmid + hr2
    #     decshift = 0.0005

    #     # Define Latitudes
    #     if latstep==1:
    #         #lat = np.arange(latmax,latmin,-1*step,dtype='float')
    #         latmin = latmid - 4*decshift
    #         latmax = latmid + 4*decshift
    #         lat = np.random.uniform(low=latmin,high=latmax,size=nev)
    #     else:
    #         #lat = np.ones(nev,dtype='float')*latmid
    #         lat = np.random.normal(loc=latmid,scale=decshift,size=nev)
    #     # Define Longitudes
    #     if lonstep==1:
    #         lonmin = lonmid - 4*decshift
    #         lonmax = lonmid + 4*decshift
    #         #lon = np.arange(lonmax,lonmin,-1*step,dtype='float')
    #         lon = np.random.uniform(low=lonmin,high=lonmax,size=nev)
    #     else:
    #         #lon = np.ones(nev,dtype='float')*lonmid
    #         lon = np.random.normal(loc=lonmid,scale=decshift,size=nev)

    #     # Set times, magnitudes, errors, residuals to defaults (all 0)
    #     dates = np.empty(nev,dtype='object')
    #     times = np.empty(nev,dtype='object')
    #     dt = datetime(2021,1,1)
    #     for i in range(0,nev):
    #         dates[i] = dt.date()
    #         times[i] = dt.time() 
    #     mag = np.ones(nev,dtype='float') 
    #     herr = np.zeros(nev,dtype='float')
    #     verr = np.zeros(nev,dtype='float') 
    #     res = np.zeros(nev,dtype='float') 

    #     log.write('\nSynthetic model events generated.\n')
    # else:
    """
    File
    """
    evf = open('event.sel','r')
    events = evf.readlines()
    evf.close()

    """
    Declare Arrays
    """
    cuspid = np.zeros(len(events),dtype='int')
    lat = np.zeros(len(events),dtype='float')
    lon = np.zeros(len(events),dtype='float')
    dates = np.empty(len(events),dtype='object')
    times = np.empty(len(events),dtype='object')
    mag = np.ones(len(events),dtype='float') 
    herr = np.zeros(len(events),dtype='float')
    verr = np.zeros(len(events),dtype='float') 
    res = np.zeros(len(events),dtype='float') 
    depth = np.zeros(len(events),dtype='float')

    # """
    # Populate Arrays
    # """
    # evcount = 0

    # minlat = 4.844E5
    # maxlat = 4.846E5
    # minlon = 6.0223E6
    # maxlon = 6.0225E6

    evcount=0
    for event in events:
        event = event.split(' ')
        event = list(filter(None,event))

        if int(event[9]) in icusp:

            # # Do DateTimes Later (these are bad in the hypoDD file)
            # u = utm.from_latlon(float(event[2]),float(event[3]))
            # if u[0]<minlat or u[0]>maxlat:
            #     continue
            # if u[1]<minlon or u[1]>maxlon:
            #     continue

            lat[evcount] = float(event[2])
            lon[evcount] = float(event[3])
            depth[evcount] = float(event[4])
            mag[evcount] = float(event[5])
            herr[evcount] = float(event[6])
            verr[evcount] = float(event[7])
            res[evcount] = float(event[8])
            cuspid[evcount] = int(event[9])

            evcount+=1

    # # Now save only a certain number of events
    # if evcount > nev:
    #     savemask = np.random.randint(0,evcount,size=nev,dtype=int)

    #     cuspid = cuspid[savemask]
    #     lat = lat[savemask]
    #     lon = lon[savemask]
    #     dates = dates[savemask]
    #     times = times[savemask]
    #     mag = mag[savemask]
    #     herr = herr[savemask]
    #     verr = verr[savemask]
    #     res = res[savemask]
    #     depth = depth[savemask]

    dt = datetime(2021,1,1)
    for i in range(0,nev):
        dates[i] = dt.date()
        times[i] = dt.time() 

    return [nev,cuspid[0:nev],lat[0:nev],lon[0:nev],dates,times,
            mag,herr,verr,res,depth]


def write_phase(log,nev,lat,lon,depth,cuspid,dates,times,mag,herr,verr,res,
               nsta,s_lab,tmp_ttp,tmp_tts,phasefile='phase.dat'):
    """
    """

    log.write('Writing phase.dat file.\n')
    """
    Write phase.dat file
    """
    phases = open(phasefile,'w')
    # Now right a phase file using this information
    for j in range(0,nev):
        phases.write('# %04i %02i %02i %02i %02i %04i %f %f %f %2.1f %2.2f %2.2f %2.2f %i \n' % 
                     (int(dates[j].year),int(dates[j].month),int(dates[j].day),
                      int(times[j].hour),int(times[j].minute),
                      int((times[j].second+times[j].microsecond/1000000)*100),
                      lat[j],lon[j],depth[j],mag[j],herr[j],verr[j],res[j],cuspid[j]))
        for m in range(0,nsta):
            phases.write('%s %f 1.0 P \n' % (s_lab[m],tmp_ttp[m,j]))
    # Copying same format to write as previous scripts
    for j in range(0,nev):
        phases.write('# %04i %02i %02i %02i %02i %04i %f %f %f %2.1f %2.2f %2.2f %2.2f %i \n' % 
                     (int(dates[j].year),int(dates[j].month),int(dates[j].day),
                      int(times[j].hour),int(times[j].minute),
                      int((times[j].second+times[j].microsecond/1000000)*100),
                      lat[j],lon[j],depth[j],mag[j],herr[j],verr[j],res[j],cuspid[j]))
        for m in range(0,nsta):
            phases.write('%s %f 1.0 S \n' % (s_lab[m],tmp_tts[m,j]))
    phases.close()

    log.write('Done writing phase.dat file. \n')

    return None


def update_phase(log,nev,lat,lon,depth,cuspid,dates,times,mag,herr,verr,res,
                 nsta,s_lab,tmp_ttp,tmp_tts,phasefile='phase.dat'):

    log.write('Writing phase.dat file.\n')

    oldpha = open(phasefile,'r')
    phas_old = oldpha.readlines()
    oldpha.close()

    pha_new = open('phase_synthetic.dat','w')

    """
    Write phase.dat file
    """
    for pha in phas_old:
        phatmp = pha
        pha = pha.strip()
        pha = pha.split(' ')
        pha = list(filter(None,pha))
        if pha[0]=='#':
            eid = int(pha[-1])
            eindx = np.argwhere(cuspid==eid)[0][0]
            pha_new.write(phatmp)

            #pha_new.write('# %s %s %s %s %s %s %s %s %f %s %s %s %s %s \n' % 
            #              (pha[1],pha[2],pha[3],pha[4],pha[5],pha[6],pha[7],pha[8],
            #               depth[eindx],pha[10],pha[11],pha[12],pha[13],pha[14]))  
        else:
            sid = str(pha[0])

            try:
                sindx = np.argwhere(s_lab==sid)[0][0]
            except:
                continue

            if pha[-1]=='P':
                pha_new.write('%s %f %s %s \n' % (pha[0],tmp_ttp[sindx,eindx],pha[2],pha[3]))
            elif pha[-1]=='S':
                pha_new.write('%s %f %s %s \n' % (pha[0],tmp_ttp[sindx,eindx],pha[2],pha[3]))
            else:
                continue

    pha_new.close()
    log.write('Done writing phase.dat file. \n')

    return None


def phasearrays(log,nev,lat,lon,depth,cuspid,dates,times,mag,herr,verr,res,
                nsta,s_lab,tmp_ttp,tmp_tts):
    """
    Sort data into phase arrays:
    ---
    For synth models all data values for all types are defined so counters needed since 
    there are no empty values
    """
    i = int(0)
    ii = int(0)
    npha = int(2*nsta*nev)

    # Declare all catalog arrays
    nobs_ct = 2*nsta*np.ones(nev,dtype='int')         # No. picks defined in synth model
    rtime = np.zeros(nev)
    p_pha = np.empty((nev,nsta*2),dtype='object')     # Data phase - character code 'P' or 'S'
    p_sta = np.empty((nev,nsta*2),dtype='object')     # Data station - alphanumeric station code ID
    p_time = np.zeros((nev,nsta*2),dtype='float')     # Arrival times
    # SET ALL DATA WEIGHTS TO 1 FOR SYNTH DATA
    p_wghtr = np.ones((nev,nsta*2),dtype='float')     # Data weights
    
    # Save traveltime formats to correct data arrays
    for iev in range(0,nev):
        for ista in range(0,nsta):
            p_pha[iev,ista] = 'P'
            p_sta[iev,ista] = s_lab[ista]
            p_time[iev,ista] = tmp_ttp[ista,iev]
            p_pha[iev,ista+nsta] = 'S'
            p_sta[iev,ista+nsta] = s_lab[ista]
            p_time[iev,ista+nsta] = tmp_tts[ista,iev]

    return [npha,nobs_ct,p_pha,p_sta,p_time,p_wghtr]


"""
Main Synthetic Model Function
"""
def synth_generate(log,syninp,datfol,fileout):
    """
    Define a synthetic model and write phase.dat file
    :::
    Parameters:
    log (file)      --- Log file
    syninp (str)    --- Synthetic model input file
    datfol (str)    --- Data folder path location
    fileout (int)   --- File output switch     
    :::
    Returns:
    retlist (list)  --- List of returned variables
    :::
    """
    
    log.write('\n\nGenerating synthetic model.\n\n')
    print('\n\nGenerating Sythetic Model')

    """
    Read in synthetic model input file
    """
    log.write('Read in synth.inp input file.\n')
    [statfile,nev,depth,latmid,lonmid,stepsize,latstep,lonstep,
     mod_nl,mod_ratio,mod_top,mod_v] = synth_input(log,syninp)

    """
    Read in station.dat file
    """
    statfile = os.path.join(datafol,'station.dat')
    readstat==bool(True)
    if readstat:
        log.write('Read in staton.dat. \n')
        [nsta,s_lab,s_lat,s_lon] = readstat(log,statfile)
    # else:
    #     nsta = 64
    #     decshift = 0.005
    #     chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    #     c1 = 0
    #     c2 = 0
    #     s_lab = np.empty(nsta,dtype='object')
    #     for i in range(64):
    #         s_lab[i] = chars[c1]+chars[c2]
    #         c2 += 1
    #         if c2>=len(chars):
    #             c1 += 1
    #             c2 = 0

    #     s_lat = np.random.uniform(low=latmid-15*decshift,high=latmid+15*decshift,size=nsta)
    #     s_lon = np.random.uniform(low=lonmid-15*decshift,high=lonmid+15*decshift,size=nsta)

    #     stats = open('station.dat','w')
    #     for i in range(nsta):
    #         stats.write('%s %f %f \n' % (s_lab[i],s_lat[i],s_lon[i]))
    #     stats.close()

    """
    Make the synthetic model starting here
    ---
    First Define Event Info Arrays
    """
    log.write('Define events. \n')
    [nev,cuspid,lat,lon,dates,times,
     mag,herr,verr,res,depth] = define_events(log,nev,latmid,lonmid,
                                              stepsize,latstep,lonstep,depth)
    # plt.plot(s_lon,s_lat,'rv')
    # plt.plot(lon,lat,'bx')
    # plt.savefig('synth_map.png')
    # plt.close('all')

    """
    Write true locations to file 
    """
    write_orig(nev,lat,lon,depth,dates,times)

    """
    Raytracing for the above problem
    ---
    Using the vmodel defined in the synth.inp file
    """
    log.write('Raytrace.')
    [tmp_ttp,tmp_tts,tmp_xp,tmp_yp,tmp_zp,
     dist,az,ang] = partials(nev,cuspid,lat,lon,depth,nsta,s_lab,s_lat,
                             s_lon,mod_nl,mod_ratio,mod_v,mod_top,return_all=True)

    """
    Write phase file for synthetic model
    ---
    Only if diskio==0
    """
    log.write('Writing phase.dat file.')
    if fileout==0:
        phasefile = os.path.join(datafol,'phase.dat')
        write_phase(log,nev,lat,lon,depth,cuspid,
                    dates,times,mag,herr,verr,res,nsta,s_lab,
                    tmp_ttp,tmp_tts,phasefile)

        print('Finished Generating Synthetic Model.\n\n')
        return []  # Nothing to return
    else:
        [npha,nobs_ct,
         p_pha,p_sta,p_time,p_wghtr] = phasearrays(log,nev,lat,lon,depth,cuspid,dates,times,
                                                   mag,herr,verr,res,nsta,s_lab,tmp_ttp,tmp_tts)

        return [nsta,s_lab,s_lat,s_lon,dist,az,ang,
                nev,lat,lon,depth,cuspid,dates,times,mag,herr,verr,res,
                npha,nobs_ct,p_pha,p_sta,p_time,p_wghtr]



