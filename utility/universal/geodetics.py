import numpy as np

"""
This script contains the universal coordinate system scripts that build
the certesian coordinate system and also convert from latlon to/form the
cartesian coordinate system.
---
This script contains the functions:
    setorg:     Initialises a cartesian coordinate system around a given point.
    dist:       Convert lat/lon to cartesian.
    redist:     Convert cartesian to lat/lon.
    sdc2:       Passes given coordinates to the needed conversion function.
---
"""

"""
##########################################################################################
######## KB Note: These functions are a remnant of original fortran code.  I decided
########            not to change it for consistency purposes across packages but it might 
########            be worthwhile looking into supported python geodetic packages (such as 
########            utm) for future version of hypoDD-py.
##########################################################################################
"""


"""
Define global geodetic coordinate system variables
"""
rearth = float(0.)
ellip = float(0.)
rlatc = float(0.)
rad = float(0.)
olat = float(0.)
olon = float(0.)
aa = float(0.)
bb = float(0.)
bc = float(0.)
sint = float(0.)
cost = float(0.)
rotate = float(0.)
icoordsystem = int(0)


def setorg(orlat,orlon,rota=0.0,ifil=0):
    """
    This function defines a cartesian coordinate system based on the given lat,lon
    The given lat,lon are set to the origin point.
    :::
    Set up cartesian coordinate system by short distance conversion
    Unrotated coordinate system with pos. x-axis toward WEST
    and pos. y-axis toward NORTH
    pos. z-axis toward EARTH'S CENTER
    :::
    PARAMETERS:
    orlat (float) --- Origin latitude
    orlon (float) --- Origin longitude
    rota (float) --- Rotation of system defaults to 0.
    ifil (int) --- Switch indicating write to file.
                   If 0, do not write coordinate system to file.
                   If not 0, write coordinate system to file. Filename defined as ifil.
    :::
    """

    # Pull in global geodetic vairables
    global rearth, ellip, rlatc, rad
    global olat, olon, aa, bb, bc
    global sint, cost, rotate
    global icoordsystem

    """
    Define origin lat/lon and rotation
    ---
    If orlat or orlon are both set to zero, the Swiss Cartesian 
    coordinate system will be used (this system cannot be rotated)
    """
    rad = 0.017453292
    if orlat==0. and orlon==0.: # Use Swiss Cartesian
        olat = 46.95240 # BERN North
        olon = -7.439583 # BERN West
        rotate = 0.
    else:   # Else define the rotation (rotation defaults to 0.)
        olat = orlat
        olon = orlon
        rotate = rota*rad

    olat = olat*60. # Minutes N
    olon = olon*60. # Minutes 

    """
    Define Earth Radius and Flattening
    ---
    New Ellipsoid for Whol Earth: WGS72
    Also set rlatc according to origin
    """
    rearth = 6378.135
    ellip = 298.26 # Flattening

    """
    Calculate Rlatc
    ---
    Conversion from geographical lat to geocentrical lat
    """
    phi = olat*rad/60.                  # phi = geogr. lat
    beta = phi-np.sin(phi*2.)/ellip     # beta = geoc. lat
    rlatc = np.tan(beta)/np.tan(phi)

    """
    Write system to file if wanted
    """
    if ifil>0:
        ifile = open(ifil,'w')
        ifile.write('SHORT DISTANCE CONVERSION on ELLIPSOID of')
        ifile.write('WORLD GEODETIC SYSTEM 1972 (WGS72) \n')
        ifile.write('======================================')
        ifile.write('======================================\n \n')
        ifile.write('(Radius at equator (rearth) = %10.5f km) \n' % rearth)
        ifile.write('(1./(ellipticity) = %10.3f \n \n' % ellip)
        ifile.write('Origin of cartesian coordinates [degrees]: \n\n')
        if orlat==0. and orlon==0.:
            ifile.write('SWISS COORDINATE SYSTEM \n\n')
            ifile.write('(Origin = city of BERN, Switzerland)\n\n')
            ifile.write('no rotation of grid, pos. y-axis toward N \n')
            ifile.write('                     pos. x-axis toward E \n\n')
        else:
            ifile.write('( %12.7f N     %12.7f W )\n\n' % (olat/60.,olon/60.))
            ifile.write(' without rotation of grid, \n')
            ifile.write('               pos. x-axis toward W \n')
            ifile.write('               pos. y-axis toward N \n\n')
            ifile.write(' Rotation of y-axis from North anticlockwise \n')
            ifile.write(' with pos. angle given in degrees \n\n')
            if rota>=0:
                ifile.write(' Rotation of grid anticlockwise by \n')
                ifile.write(' %10.5f degrees \n\n' % rota)
            else:
                ifile.write(' Rotation of grid clockwise by \n')
                arota = -1.*rota
                ifile.write(' %10.5f degrees \n\n' % arota)

    """
    Calculate aa and bb
    ---
    Length of one minute of lat and lon in km at given Earth position
    """
    lat1 = np.arctan(rlatc*np.tan(olat*rad/60.))    # geoc. lat for OLAT
    lat2 = np.arctan(rlatc*np.tan((olat+1.)*rad/60.))   # geoc. lat for (OLAT + 1 min)
    dela = lat2-lat1
    r = rearth*(1. - (np.sin(lat1)**2)/ellip)       # spherical radius for lat=olat
    aa = r*dela     # aa = 1 min geogr. lat
    delb = np.arccos(np.sin(lat1)**2 + np.cos(rad/60.)*np.cos(lat1)**2)
    bc = r*delb     # bc = 1 min geogr. lon
    bb = r*delb/np.cos(lat1)

    """
    Update file
    """
    if ifil>0:
        ifile.write('( Radius of sphere of OLAT = %10.3f km )\n\n' % r)
        ifile.write('Conversion of GEOGRAPHICAL LATITUDE to GEOCENTRICAL LATITUDE \n')
        ifile.write('RLATC = TAN(GEOCENTR.LAT) / TAN(GEOGRAPH.LAT) \n')
        ifile.write('( RLATC = %12.8f ) \n\n' % rlatc)
        ifile.write('Short distance conversions: \n')
        ifile.write('one min lat ~ %7.4f km\n' % aa)
        ifile.write('one min lon ~ %7.4 km\n' % bc)
        ifile.close()

    """
    Rotate coordinates if a rotation was defined
    """
    sint = np.sin(rotate)
    cost = np.cos(rotate)
    
    return None     # Coordinate system is defined by global geodetic variables
                    # Don't need to return anything
    

def dist(xlat,xlon):
    """
    Convert latitude and longitude to kilometers relative to
    center of coordinates by short distance conversion.
    :::
    PARAMETERS:
    xlat (float) ---- Latitude of point
    ylon (float) ---- Longitude of point
    :::
    RETURNS:
    xkm (float) ---- X coordinate in cartesian using cluster centroid coordinate system
    ykm (float) ---- Y coordinate in cartesian using cluster centroid coordinate system
    :::
    """

    """
    Pull in geodetic global variables
    """
    global rearth, ellip, rlatc, rad
    global olat, olon, aa, bb, bc
    global sint, cost, rotate
    global icoordsystem
    
    """
    Set up short distance conversion by subtracting SETORG coordinates
    ---
    The coordinate system need to be predefined to run.
    """
    q = 60.*xlat - olat
    yp = q+olat
    lat1 = np.arctan(rlatc*np.tan(rad*yp/60.))
    lat2 = np.arctan(rlatc*np.tan(rad*olat/60.))
    lat3 = (lat2+lat1)/2.
    xx = 60.*xlon-olon
    q = q*aa
    xx = xx*bb*np.cos(lat3)
    if rotate!=0:
        yp = cost*q*sint*xx
        xx = cost*xx-sint*q
        q=yp

    xkm=xx
    ykm=q

    return xkm,ykm  # Return x and y in cartesian


def redist(xkm,ykm):
    """
    Convert from local Cartesian coordinates to lat and lon
    :::
    PARAMETERS:
    xkm (float) ---- X location using cluster centroid coordinate system
    ykm (float) ---- Y location using cluster centroid coordinate system
    :::
    RETURNS:
    xlat (float) ---- Latitude of point (x,y)
    xlon (float) ---- Longitude of point (x,y)
    :::
    """

    """
    Pull in global geodetic variables
    """
    global rearth, ellip, rlatc, rad
    global olat, olon, aa, bb, bc
    global sint, cost, rotate
    global icoordsystem

    xx = xkm
    yy = ykm
    """
    Rotate coordinates anticlockwise back into lat/lon
    """
    y = yy*cost-xx*sint
    x = yy*sint+xx*cost
    if abs(aa)>0.0000001:
        q = y/aa
        lat = (q+olat)/60.
        xlat = q+olat - 60.*lat
        yp = 60.*lat + xlat
        lat1 = np.arctan(rlatc*np.tan(yp*rad/60.))
        lat2 = np.arctan(rlatc*np.tan(olat*rad/60.))
        lat3 = (lat1+lat2)/2.
        clat1 = np.cos(lat3)
        bcl = bb*clat1
        if abs(bcl)>0.000001:
            p = x/(bb*clat1)
            lon = (p+olon)/60.
            xlon = p+olon - 60.*lon
            xlat = lat + xlat/60.
            xlon = lon + xlon/60.
            return xlat,xlon    # If able to conver back return values

    """
    If there was an error in conversion print values and break
    """
    print('subr. redist: \n')
    print('aa = %10.5f \n' % aa)
    print('bb = %10.5f \n' % bb)
    print('cos(lat1) = %10.5f \n' % clat1)

    raise Exception('division by zero run stops here')


def sdc2(xlat,xlon,i):
    """
    Call coordinate conversion subfunctions
    :::
    Convert coordinates of a point by short distance conversion.
    This function either passes the given x,y or lat,lon into:
        redist (i==1) - converts from xy to latlon
        or
        dist (i==-1) - converts from latlon to xy
    Parameter and return values depend on direction of conversion.
    :::
    PARAMETERS:
    xlat (float) ---- Latitude or X
    xlon (float) ---- Longitude or Y
    i (int) --- Switch indicating conversion direction
    :::
    RETURNS:
    x (float) --- Latitude or X
    y (float) --- Longitude or Y
    :::
    """
    if i!=1 and i!=-1:
        # Invalid switch token
        raise Exception('SDC: Specify conversion')
    if i==1:
        # Convert from x,y to lat,lon
        x,y, = redist(xlat,xlon)
    if i==-1:
        # Convert from lat,lon to x,y
        x,y = dist(xlat,xlon)

    return x,y

