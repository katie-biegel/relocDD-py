import numpy as np

""" 
This script contains universal functions called from 
multiple other functions.
---
This script contains the functions:
    atoangle:   Convert from string to decdeg format if needed.
    delaz:      Computes distance and azimuth on a sphere.
    readstat:   Read in the station.dat file
---
"""

def atoangle(locstr):
    """
    Convert string from form "degrees:minutes:seconds"
    to angle if needed. Otherwise just return float of angle.
    :::
    PARAMETERS:
    locstr (str) ---- Lat or lon string from file
    :::
    RETURNS:
    loc (float) ---- Lat or lon in angle form
    :::
    """
    if ':' in locstr:   # If in degree minutes form split and convert 
                        # to decimal degree
        loc = locstr.split(':')
        loc = float(loc[0]) + (float(loc[1]) + float(loc[2])/60.)/60.
    else:       # Otherwise just read in float value
        loc = float(locstr)

    return loc      # Return value
    

def delaz(alat, alon, blat, blon):
    """
    This function computes distance and azimuth on a sphere
    :::
    PARAMETERS:
    alat (float) ---- Latitude of first point
    alon (float) ---- Longitude of first point
    blat (float) ---- Latitude of second point
    blon (float) ---- Longitude of second point
    :::
    RETURNS:
    delt (float) ---- Central angle (degrees)
    dist (float) ---- Distance (km)
    az (float) ---- Azimuth from a to b (degrees)
    ::
    """
    # Variables declared in delaz2.f (original fortran)
    # Kept consistent to retain same values no
    pi2 = 1.570796
    rad = 1.745329e-2
    flat = .993231

    # Convert lat/lon from degrees to radians
    alatr = alat*rad
    alonr = alon*rad
    blatr = blat*rad
    blonr = blon*rad

    # Convert radian latitudes to geocentric colatitudes
    # North pole 0 to 2pi south pole (0 to 180 in degrees)
    tana = flat*np.tan(alatr)
    geoa = np.arctan(tana)
    acol = pi2 - geoa
    tanb = flat*np.tan(blatr)
    geob = np.arctan(tanb)
    bcol = pi2 - geob

    # Calculate delta (the central angle) in radians
    diflon = blonr-alonr
    cosdel = np.sin(acol)*np.sin(bcol)*np.cos(diflon) + np.cos(acol)*np.cos(bcol)
    delr = np.arccos(cosdel)

    # Calculate azimuth from a to b
    top = np.sin(diflon)
    den = (np.sin(acol)/np.tan(bcol)) - np.cos(diflon)*np.cos(acol)
    azr = np.arctan2(top,den)

    # Convert back into degrees
    delt = delr/rad
    az = azr/rad
    if az < 0.0:
        az = 360.+az

    # Compute distance in km
    colat = pi2 - (alatr+blatr)/2.
    # The equatorial radius of the Earth is 6378.137 km (IUGG value)
    # The mean equatorial radius from Bott 1982 is 6378.140 km
    # KB Note: The Bott conversion is used in the original hypoDDv1.3
    # but we did notice an integer division issue here.  Mathematically,
    # line 132 below is correct but may produce slightly different results 
    # due to the integer division issue in the original fortran function.
    radius = 6378.140*(1.0+3.37853e-3*(1./3.-((np.cos(colat))**2)))
    dist = delr*radius

    return delt,dist,az # Return the distance values


def readstat(log,statfile,maxdist=-9,clat=0.,clon=0.,fileout=0):
    """
    Open and read station file
    :::
    PARAMETERS:
    log (file object) ---- Code log file
    statfile (str) ---- File location for stat file
    maxdist (float) ---- Max. cluster centroid to sta. separation
    clat (float) --- Cluster centroid latitude
    clon (float) --- Cluster centroid longitude
    fileout (int) ---- Integer switch for file inputs/outputs
                        If 0, run traditional hypoDD
                        If 1, limit input/output.  Used for bootstrapping.
    :::
    RETURNS:
    nsta (int) --- Number of stations
    s_lab[nsta] (object array) --- List array of station codes
    s_lat[nsta] (float array) --- Array of station latitudes
    s_lon[nsta] (float array) --- Array of station longitudes
    :::
    """

    stations = open(statfile,'r')

    """
    Read station.dat file to create station information arrays
    """
    try:
        stats = stations.readlines()
    except:     # Break if can't read file.
        raise Exception('Error reading station file.')
    """
    Save stations to arrays
    """
    i = int(0)
    nsta = len(stats)                           # nsta  ---> No.of stations total in file (can include duplicates)
    s_lab = np.empty(nsta,dtype='U7')           # s_lab ---> Station label (usually character string or alphanumeric code)
    s_lat = np.zeros(nsta,dtype='float')        # s_lat ---> Station latitudes
    s_lon = np.zeros(nsta,dtype='float')        # s_lon ---> Station longitudes
    for ind,sta in enumerate(stats):
        sta = sta.split()
        sta = list(filter(None,sta))
        try:
            s_lab[i] = str(sta[0])
            s_lat[i] = atoangle(sta[1])         # Run through atoangle function (top of script) to ensure correct 
            s_lon[i] = atoangle(sta[2])         # decimal degree formatting
        except:
            raise Exception('Error reading station file. Line %i' % i)
        else:
            if maxdist==-9:
                i += 1
                continue
            else:
                # Skip at distances larger than maxdist
                # Only save stations with a good event-station separation
                delt,dist,azim = delaz(clat,clon,s_lat[i],s_lon[i])
                if dist <= maxdist:
                    i += 1
    """
    Correct nsta and lengths
    ---
    Only if maxdist defined
    """
    if maxdist!=-9:
        log.write('> Percent stations kept (<maxdist): %f' % (float(i/nsta)))
        nsta = i
        s_lab = s_lab[0:nsta]
        s_lat = s_lat[0:nsta]
        s_lon = s_lon[0:nsta]
    """   
    Log total number of stations
    """
    log.write('> stations = %i \n' % nsta)
    stations.close()
    
    """
    Check for double stations (if double entry exists, exit code and user
    is required to fix that in file):
    """
    for i in range(nsta):
        if s_lab[i] in s_lab[:i]:
            log.write('This station is listed twice: %s \n' % s_lab[i])
            raise Exception('Station listed twice: %s' % s_lab[i])

    """
    Return station arrays
    """
    return nsta,s_lab,s_lat,s_lon






    