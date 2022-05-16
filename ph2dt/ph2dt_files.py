# Built in python packages
import numpy as np
import os
import sys
import datetime

# Other hypoDD functions
from utility.universal.misc import atoangle
#from memory_profiler import profile

"""
This script contains ph2dt subfunctions for reading in files or
for writing to file if needed.
----
This script includes the following functions:
    ph2dt_input:    Reading in to memory the ph2dt input file (ph2dt.inp)
    readeventsel:   Read in event.select file for event subsets in ph2dt
    readphase:      Read in the phase.dat file
----
"""
#@profile(stream=open('mem_logs/ph2dt_input.mem','w+'))
def ph2dt_input(log,inputfile='ph2dt.inp',dd_version=1):
    """
    Reads ph2dt input file and returns variables
    or 
    Prompt user to enter ph2dt input file variables
    :::
    Test open input parameter file and read in variables
    If the file specified above cannot be opened, then prompt user for all inputs.
    Variables in the input files also have default values that can be seen here (in case read-in breaks).
    :::
    PARAMETERS:
    log (file object) ---- Log file
    inputfile (str) ---- File string to input file location [default='ph2dt.inp']
    dd_version (int) ---- Switch for the type of double-difference relocation.
                          [default=1]
                          1 for event-pair
                          2 for station-pair
                          3 for double-pair
    :::
    RETURNS:
    retlist (list) ---- Contains the variables read in from the input file
    :::
    """

    """
    Read input file
    """
    inp = open(inputfile,'r')
    inputs = inp.readlines()
    inp.close()
    l = 1

    """
    Set maxoffset variables to default
    :::
    The maxoffset is the pairing offset. Depending on the DD
    relocation type one of these values may be unused.
    If unused the value is -9.
    """
    maxoffsete = -9
    maxoffsets = -9
    try:
        """
        Read in lines from input file
        """
        for line in inputs:
            if line[0:1] == '*' or line[1:2]=='*': # Comment lines
                continue
            else:
                if l==1:    # Station file declaration
                    line = line.split('\n')
                    statfile = str(line[0])
                if l==2:    # Phase file declaration
                    line = line.split('\n')
                    phasefile = str(line[0])
                if l==3:    # All specified pairing variables
                    line = line.split()
                    line = list(filter(None,line)) 
                    minwght = float(line[0])    # Min datum weight to be included
                    maxdist = float(line[1])    # Max dist b/w event-pair centroid and station for data pairing
                    maxoffsete = float(line[2])  # Max interevent distance for event-pair pairing
                    maxoffsets = float(line[3])  # Max interevent distance for station-pair pairing
                    mnb = int(line[4])           # Max no. of neighbors for clustering
                    limobs_pair = int(line[5])   # No. of shared data for a strong pair classifying
                    minobs_pair = int(line[6])   # Min no. of shared data to be a pair
                    maxobs_pair = int(line[7])   # Max no. of shared data that is written to file
                l+=1  
    except: # If there is a formatting error in the ph2dt input file we break.
        raise Exception('Error reading input file. Line %i \n %s' % (l,line))

    """
    Update log file
    """
    log.write('Ph2dt input file vairables: \n')
    log.write('Data file locations: statfile=%s, phasefile=%s \n' %
              (statfile,phasefile))
    log.write('Pairing variables: minwght=%f, maxdist=%f, maxoffsete=%f, maxoffsets=%f, \n' %
              (minwght,maxdist,maxoffsete,maxoffsets))
    log.write('mnb=%i, limobs_pair=%i, minobs_pair=%i, maxobs_pair=%i \n' %
              (mnb,limobs_pair,minobs_pair,maxobs_pair))
    log.write('Returning ph2dt input variables.\n\n')

    """
    Build return list from input variables
    """
    retlist = [statfile, phasefile, minwght, maxdist, maxoffsete, maxoffsets,
               mnb, limobs_pair, minobs_pair, maxobs_pair]
    return retlist

#@profile(stream=open('mem_logs/readeventsel.mem','w+'))
def readeventsel(log,fileout,selev='events.select'):
    """
    Determine if event subset is to be used from event.select file (fn9 file)
    (Not used everytime - only used if you want to look at a subset of a larger catalog without recreating input files)
    :::
    Parameters:
    log (file object) --- Log file
    fileout (int) --- Input/output switch
    selev (str) --- File location for selected event list (if defined)
    :::
    Returns:
    None
    or 
    ncusp (int) --- Number of selected events
    icusp[ncusp] (int array) --- Array of event IDs for selected events
    :::
    """
    ncusp = 0
    try:
        sel_ev = open(selev,'r')
        try:
            sel_ev = sel_ev.readlines()
        except:
            raise Exception('Error reading event.select file.')
        ncusp = len(sel_ev)
        icusp = np.zeros(nev)
        for i,ev in enumerate(sel_ev):
            try:
                icusp[i] = int(ev)
            except:
                raise Exception('Error reading ID number. Line %i' % i)
        return ncusp,icusp
    except:
        log.write('No events.select file. \n')
        if fileout==0:
            print('No event.select file.')
        return 0, []


def readphase(log,phasefile,datfol,minwght,minobs_pair,nsta,fileout=0,MAXEV=10000):
    """ 
    Opens and reads the phasefile specified
    Returns phase arrays.
    Writes event files from phasefile information
    :::
    PARAMETERS:
    log (file object) ---- Code log file
    phasefile (str) ---- File location for phase file
    fileout (int) ---- Integer switch for file inputs/outputs
                      If 0, run traditional hypoDD
                      If 1, limit input/output.  Used for bootstrapping.
    :::
    RETURNS:
    retlist (list) ---- List of event and phase information
    :::
    """
    
    """
    Check for selected event list
    """
    selev = readeventsel(log,fileout)
    if selev:
        ncusp = selev[0]
        icusp = selev[1]

    """
    Open phase file to read it into memory
    """
    phases = open(phasefile,'r')
    if fileout==0:   # Only save events to file if i/o turned on.
        efile = os.path.join(datfol,'event.dat')
        events = open(efile,'w')  # Event list
        eselfile = os.path.join(datfol,'event.sel')
        ev_sel = open(eselfile,'w')  # Subselection of event list (not always different)

    """
    Read absolute network arrivaltimes from phase.dat file and save into phase information 
    arrays and event information arrays
    ---
    arrays starting p_* indicate phase information from this loop
    """
    # Try reading file and break if unable
    try:
        if fileout==0:
            print('reading phase file...')
        phas = phases.readlines()
    except:
        raise Exception('Error reading phase file.')
    """
    Declare all event arrays
    """
    dates = np.empty(MAXEV,dtype='object')    
    times = np.empty(MAXEV,dtype='object') 
    lat = np.zeros(MAXEV,dtype='float')     # Event latitudes (for pre-determined event locations)
    lon = np.zeros(MAXEV,dtype='float')     # Event longitudes
    depth = np.zeros(MAXEV,dtype='float')   # Event depths
    cuspid = np.zeros(MAXEV,dtype='int')    # Event IDS - usually a network generated code (must be numeric)
    mag = np.zeros(MAXEV,dtype='float')     # Event Magnitudes
    herr = np.zeros(MAXEV,dtype='float')    # Horizontal Location Errors
    verr = np.zeros(MAXEV,dtype='float')    # Vertical Location Errors  
    res = np.zeros(MAXEV,dtype='float')     # Starting event traveltime residual (RMS)
    """
    Declare all phase arrays
    """
    nobs_ct = np.zeros(MAXEV,dtype='int')   # No. of arrivaltime data per event (allows for variable number of data per event)
    MAXOBS = int(nsta*2)              # Total no. of possible data
    p_pha = np.empty((MAXEV,MAXOBS),dtype='U1')     # Data phase - character code 'P' or 'S'
    p_sta = np.empty((MAXEV,MAXOBS),dtype='U7')     # Data station - alphanumeric station code ID
    p_time = np.zeros((MAXEV,MAXOBS),dtype='float')     # Arrival times
    p_wghtr = np.zeros((MAXEV,MAXOBS),dtype='float')    # Data weights

    """
    Declare all loop variables
    """
    i = int(0)
    ii = int(0)
    npha = int(0)
    count = int(0)
    ip = int(0)
    """
    Read and process phases
    """
    for pha in phas:
        pha = pha.strip()
        #print(pha)
        pha = pha.split()
        pha = list(filter(None,pha))
        #print(pha)
        #import pdb; pdb.set_trace()

        if pha[0]=='#':     # A '#' indicates an event line
            if count!=0:
                skip = int(0) # Indicator of duplicate event; changes to event index if event
                             # has already been read in.  A warning will appear for all duplicate
                             # events it is the user's responsibility to check to make sure there 
                             # are not duplicate arrivals.  This will mess up relocation.
                if i==0: 
                    """
                    If this is the first event:
                    """
                    # Store Previous Event
                    nobs_ct[i] = k
                    itake = 1
                    # Event Selection
                    if ncusp > 0:
                        itake = 0
                        if cuspid[i] in icusp:
                            itake = 1
                    if fileout==0:
                        # Right previous event to file
                        # Write header to total event list
                        events.write('%08i %08i %5.6f %5.6f %5.6f %1.2f %6.2f %6.2f %6.2f %8i \n' %
                                     (evdate,evtime,lat[i],lon[i],depth[i],mag[i],herr[i],verr[i],res[i],cuspid[i]))
                    ii = ii + 1
                    # If the event has the minimum no. of obs also right to the event.sel file
                    # Keep event only if min_obs met
                    if itake==1: #and nobs_ct[i]>=int(minobs_pair):
                        if fileout==0:
                            ev_sel.write('%08i %08i %5.6f %5.6f %5.6f %1.2f %6.2f %6.2f %6.2f %8i \n' %
                                         (evdate,evtime,lat[i],lon[i],depth[i],mag[i],herr[i],verr[i],res[i],cuspid[i]))
                        npha = npha+nobs_ct[i]
                        i = i+1
                else:
                    """"
                    Check for duplicate events
                    ---
                    If event is a duplicate the data is saved but a warning is printed.
                    User has to check for duplicate measurements themselves.
                    """
                    if cuspid[i] in cuspid[0:i]:  
                        """
                        Do not duplicate events but this still save the data
                        """
                        log.write('Duplicate event: %i. Check i/o files for duplicate arrivals.' % cuspid[i])
                        if fileout==0:
                            print('Duplicate event: %i. Check i/o files for duplicate arrivals.' % cuspid[i])
                        skip = np.where(cuspid[0:i]==cuspid[i])[0][0]
                        # Store Previous Event Phases in Correct Location
                        tmp = int(nobs_ct[skip])
                        p_sta[skip,tmp:tmp+k] = p_sta[i,0:k]
                        p_time[skip,tmp:tmp+k] = p_time[i,0:k]
                        p_wghtr[skip,tmp:tmp+k] = p_wghtr[i,0:k]
                        p_pha[skip,tmp:tmp+k] = p_pha[i,0:k]
                        nobs_ct[skip] = nobs_ct[skip]+k
                        npha = npha+k
                    else:
                        """
                        Store Previous Event
                        """
                        nobs_ct[i] = k
                        itake = 1
                        # Event Selection
                        if ncusp > 0:
                            itake = 0
                            for k in range(0,ncusp):
                                if cuspid[i] == icusp[k]:
                                    itake = 1
                        if fileout==0:
                            # Right previous event to file
                            # Write header to total event list
                            events.write('%08i %08i %5.6f %5.6f %5.6f %1.2f %6.2f %6.2f %6.2f %8i \n' %
                                         (evdate,evtime,lat[i],lon[i],depth[i],
                                          mag[i],herr[i],verr[i],res[i],cuspid[i]))
                        ii = ii + 1
                        # If the event has the minimum # of obs also right to the event.sel file
                        # Keep event only if min_obs met
                        if itake==1: # and nobs_ct[i]>=minobs_pair:
                            if fileout==0:
                                ev_sel.write('%08i %08i %5.6f %5.6f %5.6f %1.2f %6.2f %6.2f %6.2f %8i \n' %
                                             (evdate,evtime,lat[i],lon[i],depth[i],
                                              mag[i],herr[i],verr[i],res[i],cuspid[i]))
                            npha = npha+nobs_ct[i]
                            i = i+1
            
            """
            Read in new event variables
            """
            # Event Date
            yr = int(pha[1])
            mo = int(pha[2])
            dy = int(pha[3])
            evdate = int(yr*10000 + mo*100 + dy)
            dates[i] = datetime.date(year=yr,month=mo,day=dy)

            # Event Time
            hr = int(pha[4])
            minute = int(pha[5])
            #sec = float(pha[6])
            sec = int(pha[6][0:2])+float(pha[6][2:])/100.
            evtime = int(hr*1000000 + minute*10000 + sec*100)
            times[i] = datetime.time(hour=hr,minute=minute,second=int(sec),microsecond=int((sec%1)*1000000))

            # Event Location
            lat[i] = atoangle(pha[7])
            lon[i] = atoangle(pha[8])
            depth[i] = float(pha[9])

            # Other Event information
            mag[i] = float(pha[10])
            herr[i] = float(pha[11])
            verr[i] = float(pha[12])
            res[i] = float(pha[13])
            cuspid[i] = float(pha[14])

            k = 0  # Phase counter
            count += 1   
            #i += 1 
        else:
            """
            Read in Phase Line
            """
            p_sta[i,k] = str(pha[0])
            p_time[i,k] = float(pha[1])
            p_wghtr[i,k] = float(pha[2])
            p_pha[i,k] = str(pha[3])
            if p_pha[i,k]=='P' or p_pha[i,k]=='S':
                if p_wghtr[i,k]>=minwght or p_wghtr[i,k]<0:
                    if p_pha[i,k]=='P':
                        ip += 1
                    k += 1

    """
    Store Last Event 
    ---
    (Since events are written to file post-loop)
    """
    skip = int(0)
    if cuspid[i] in cuspid[0:i]:  
        """
        Duplicate events are not saved separately
        ---
        User must check for duplicate data.
        """
        log.write('Duplicate event: %i. Check i/o files for duplicate arrivals.' % cuspid[i])
        if fileout==0:
            print('Duplicate event: %i. Check i/o files for duplicate arrivals.' % cuspid[i])
        skip = np.where(cuspid[0:i]==cuspid[i])[0][0]
        # Store Previous Event Phases in Correct Location
        tmp = int(nobs_ct[skip])
        p_sta[skip,tmp:tmp+k] = p_sta[i,0:k]
        p_time[skip,tmp:tmp+k] = p_time[i,0:k]
        p_wghtr[skip,tmp:tmp+k] = p_wghtr[i,0:k]
        p_pha[skip,tmp:tmp+k] = p_pha[i,0:k]
        nobs_ct[skip] = nobs_ct[skip]+k
        npha = npha+k
    else:
        """
        If not duplicate event then just store it
        """
        nobs_ct[i] = k
        itake = 1
        # Event Selection
        if ncusp > 0:
            itake = 0
            if cuspid[i] in icusp:
                itake = 1
        if fileout==0:
            # Right previous event to file
            # Write header to total event list
            events.write('%08i %08i %5.6f %5.6f %5.6f %1.2f %6.2f %6.2f %6.2f %8i \n' %
                         (evdate,evtime,lat[i],lon[i],depth[i],
                          mag[i],herr[i],verr[i],res[i],cuspid[i]))
        ii = ii + 1
        # If the event has the minimum # of obs also right to the event.sel file
        # Keep event only if min_obs met
        if itake==1: # and nobs_ct[i]>=minobs_pair:
            if fileout==0:
                ev_sel.write('%08i %08i %5.6f %5.6f %5.6f %1.2f %6.2f %6.2f %6.2f %8i \n' %
                             (evdate,evtime,lat[i],lon[i],depth[i],
                              mag[i],herr[i],verr[i],res[i],cuspid[i]))
            npha = npha+nobs_ct[i]
            i = i+1

    """
    Log Phases
    --
    Print to terminal/log file data statistics (total number of events and phases):
    """
    nev = int(i)
    if fileout==0:
        print('> events total = %i' % (ii))
        print('> events selected = %i' % nev)
        print('> phases = %i' % npha)
        print('> P phases = %i' % ip)
        print('> S phases = %i' % (npha-ip))
    log.write('> events total = %i \n' % (ii))
    log.write('> events selected = %i \n' % nev)
    log.write('> phases = %i \n' % npha)
    log.write('> P phases = %i \n' % ip)
    log.write('> S phases = %i \n' % (npha-ip))

    """
    Trim Arrays To Clear Memory
    """
    dates = dates[0:nev]
    times = times[0:nev]
    lat = lat[0:nev]
    lon = lon[0:nev]
    depth = depth[0:nev]
    cuspid = cuspid[0:nev]
    mag = mag[0:nev]
    herr = herr[0:nev]
    verr = verr[0:nev]
    res = res[0:nev]
    nobs_ct = nobs_ct[0:nev]

    # The phase arrays are saved to the max no. of data per event
    # There can be zeros in this array.
    obs = np.amax(nobs_ct)
    npha = np.sum(nobs_ct)
    #print('obs: ',obs)
    #print('npha: ',npha)
    p_pha = p_pha[0:nev,0:obs]
    p_sta = p_sta[0:nev,0:obs]
    p_time = p_time[0:nev,0:obs]
    p_wghtr = p_wghtr[0:nev,0:obs]

    #npha=0
    #for i in range(nev):
    #    npha+=nobs_ct[i]
    #print('P phases: ',np.count_nonzero(p_pha=='P'))
    #print('S phases: ',np.count_nonzero(p_pha=='S'))

    """
    Close Files
    """
    phases.close()
    if fileout==0:
        events.close()
        ev_sel.close()

    """
    Return phase and event information
    """
    retlist = [nev,lat,lon,depth,cuspid,dates,times,mag,herr,verr,res,
               npha,nobs_ct,p_pha,p_sta,p_time,p_wghtr]
    return retlist


