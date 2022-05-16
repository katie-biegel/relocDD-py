# Import General Packages
import numpy as np
import os
from datetime import datetime,timedelta
import subprocess

# Other hypoDD functions
from hypoDD.hypoDD_files import hypoDD_input
from ph2dt.ph2dt_files import readphase
from utility.universal.misc import readstat
from utility.universal.raytrace import partials
from utility.synthetics.synthModel import synth_input

# Insert bin location to user's HYPOINVERSE2000 installation
hypoinvbin = '/Users/katie/Installations/hypoinv/bin/'

"""
The purpose of this script is use in synthetic model testing. This script contains 
functions that build the correct inputs for running hypoinverse-2000 as an 
intermediate catalog building step for relocation testing.
---
This script contains the utility functions:
    decdeg2dms:         Converts from decimal degrees to degree minutes,
    dms2decdeg:         Converts from degree minutes to decimal degrees,
    make_executable:    Changes file permissions to executable;
The hypoinverse subfunctions:
    velfile:            Writes the hypoinverse velocity file,
    stafile:            Writes the hypoinverse station file,
    phafile:            Writes the hypoinverse phase file,
    hypfile:            Writes the hypoinverse run file,
    runbash:            Writes an executes a bash script wrapper for hypoinverse,
    readcat:            Read back into memory the hypoinverse locations;
And the main hypoinverse wrapper:
    hypoinverse:        Main hypoinv run function
                        Calls the necessary functions to build input files (velfile,
                        stafile, phasefile, hypfile), to write and execute the bash
                        script (runbash), and to read back in the updated event 
                        locations (readcat).
    hypoinverse_wrapper:     Runs hypoinverse but to/from file.
---
"""

"""
Utility Functions
"""

def decdeg2dms(dd):
    """
    Convert from decimal degrees to degree minutes
    KB Note: Take and edited from stackoverflow example
    :::
    PARAMETERS:
    dd(float) ---- Decimal degree value
    :::
    RETURNS:
    degrees (int) ---- Degree value
    minutes (int) ---- Minute value
    :::
    """
    if dd>=0:
        degrees,minutes = divmod(dd*60,60)
    else:
        dd = abs(dd)
        degrees,minutes = divmod(dd*60,60)
        degrees = -degrees
    return (degrees,minutes)


def dms2decdeg(degrees,minutes):
    """
    Convert from degree minutes to decimal degrees
    :::
    PARAMETERS:
    degrees (int) ---- Degree value
    minutes (int) ---- Minute value
    :::
    RETURNS:
    dd(float) ---- Decimal degree value
    :::
    """
    dd = degrees
    decimal = minutes/60
    dd = dd+decimal

    return (dd)


def make_executable(path):
    """
    This function makes a file executable.
    KB Note: Function taken from stackoverflow example
    :::
    PARAMETERS:
    path (str) --- Path to file
    :::
    """
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2    # copy R bits to X
    os.chmod(path, mode)

    return None


"""
Input File Building Functions
"""

def velfile(mod_nl,mod_v,mod_top):
    """
    This function writes the velocity model input file
    for hypoinverse-2000
    :::
    This works for a layered vp model
    :::
    PARAMETERS:
    mod_nl (int) - number of velocity model layers
    mod_v (float array) - contains the velocity of each model layer 
    mod_top(float array) - contains the depth to top of each model layer
    :::
    """
    # Check for output directory
    if not os.path.isdir('synth_hypoINV'):
        os.mkdir('synth_hypoINV')

    """
    Write HYPOINV Velocity File
    """
    vpmodel = open('synth_hypoINV/synth_P.crh','w')
    # Opening line
    vpmodel.write('SYNTHETIC VP MODEL\n')
    # Add layers
    for ilay in range(0,mod_nl):
        vpmodel.write(' %3.2f %3.2f\n' % (mod_v[ilay],mod_top[ilay]))
    vpmodel.close() # Close files and return empty

    return None


def stafile(nsta,s_lab,s_lat,s_lon):
    """
    This function writes the station file for hypoinverse-2000 from
    numpy arrays containing the station information for the synthetic problem.
    :::  
    We're assuming all stations at the surface and no elevation needed to be
    taken into account in this problem.
    :::
    PARAMETERS:
    nsta (int) - number of stations
    s_lab (obj array) - numpy array containing strings of the station name
    s_lat (float array) - numpy array containing station latitudes
    s_lon (float array) - numpy array containing station longitudes
    :::
    """

    # Check for output directory
    if not os.path.isdir('synth_hypoINV'):
        os.mkdir('synth_hypoINV')

    """
    Write HYPOINV Station File
    """

    # Station data we have: nsta,s_lab,s_lat,s_lon
    stfile = open('synth_hypoINV/synth_stat.sta','w')
    # Convert from lat/lon decimal to degree minutes
    for ista in range(0,nsta):
        # Hypoinverse requires the degree minutes lat/lon format
        [latdeg,latmin] = decdeg2dms(s_lat[ista])
        [londeg,lonmin] = decdeg2dms(s_lon[ista])
        """
        Write each station to file
        ---
        This formatting can be found in the hypoinverse-2000 documentation.
        If you get issues when running hypoinverse, the problem is likely formatting
        of the input files.  Be sure to double-check here that formatting for the 
        station file is correct.
        """
        stfile.write('%s SY  EHZ  %02i %02.4f %03i %02.4f   000.0     0.00  0.00  0.00  0.00 0  0.00--EHZ\n' %
                     (s_lab[ista].ljust(5,' '),int(latdeg),latmin,int(londeg),lonmin))
    stfile.close() # Close file and return empty

    return None


def phafile(p_time,dist,az,ang,
            nsta,s_lab,s_lat,s_lon,
            nev,cuspid,lat,lon,depth,
            dates,times,mag,herr,verr,res,):

    """
    This function writes the phase file for hypoinverse-2000 from
    numpy arrays containing the event, station, and traveltime information
    in the format generated in hypoDD.
    :::
    This can be written in many formats: hypo71, hypo71 shadow, or other archive formats
    It was hard to find documentation for this but here is the shadow format guidlines:
        http://www.ncedc.org/ftp/pub/doc/ncsn/shadow2000.pdf
    And a vague description of the phase lines of the phase format can be found in the
    hypoinverse documentation.  I also used the test file in this situation as the baseline
    and simplified it as much as possible to this current working format.
    Typically the phase format is in this format for each event:
       Event Info Line (called Header Line)
       Phase Line 1 
       Phase Line 2
       ...
       Phase Line N
                   Terminator Line
    Multiple events can be combined into one phase file so you only need one per synthetic example.
    :::
    PARAMETERS:
    p_time[nev,2*nsta] (float array) --- Traveltimes
    dist[nsta,nev] (float array) --- Contains the event-station dist in km
    az[nsta,nev] (float array) --- Contains the event-station azimuth
    ang[nsta,nev] (float array) --- Contains the ray take-off angle
    nsta (int) --- number of stations
    s_lab (obj array) --- numpy array containing strings of the station name
    s_lat (float array) --- numpy array containing station latitudes
    s_lon (float array) --- numpy array containing station longitudes
    nev (int) --- number of events
    lat (float array) --- numpy array containing event latitudes
    lon (float array) --- numpy array containing event longitudes
    depth (float array) --- numpy array containing event depths
    dates (object array) --- numpy array containing event dates as datetime.date
    times (object array) --- numpy array containing event times as datetime.time
    mag (float array) --- numpy array containing the event magnitudes
    herr (float array) --- numpy array containing the horizontal event location errors
    verr (float array) --- numpy array containing the vertical event location errors
    res (float array) --- numpy array containing the traveltime residuals for each event
    :::
    """

    # Check for file directory
    if not os.path.isdir('synth_hypoINV'):
        os.mkdir('synth_hypoINV')

    # First open file
    phafile = open('synth_hypoINV/synth_pha.phs','w')
    strterm = ' '*62    # Placeholder empty string

    # Now loop over events and stations to write information into file:
    for iev in range(nev):
        # First write header line for each event which is:
        #   1) A summary line with location and event data
        #   2) IF IN SHADOW FORMAT: Includes reference time and archive tape info
        #   3) Can include optional events if it starts with "$2"
        [latdeg,latmin] = decdeg2dms(lat[iev])
        [londeg,lonmin] = decdeg2dms(lon[iev])

        # Event Location String
        dtimesec = times[iev].second + times[iev].microsecond/1000000

        """
        #################
        # Header String #
        #################
        Locinfo: 'YYYYMMDYHRMNSECSLT LTMNLON LNMNDEPTHMAGNUMAZIDST
        This first part of the header string includes:
        Col.    Len.    Format      
        1       4       I4          Year. *
        5       8       4I2         Month, day, hour and minute.
        13      4       F4.2        Origin time seconds.
        17      2       F2.0        Latitude (deg). First character must not be blank.
        19      1       A1          S for south, blank otherwise.
        20      4       F4.2        Latitude (min).
        24      3       F3.0        Longitude (deg).
        27      1       A1          E for east, blank otherwise.
        28      4       F4.2        Longitude (min).
        32      5       F5.2        Depth (km).
        37      3       F3.2        Magnitude from maximum S amplitude from NCSN stations *
        40      3       I3          Number of P & S times with final weights greater than 0.1.
        43      3       I3          Maximum azimuthal gap, degrees.
        46      3       F3.0        Distance to nearest station (km). 
        #################
        """
        locinfo = ('%4i%02i%02i%02i%02i%04i%02i %04i%03i %04i%5i100%3i%3i%3i' %
                   (dates[iev].year,dates[iev].month,dates[iev].day,
                    times[iev].hour,times[iev].minute,int(dtimesec*100),
                    int(latdeg),int(latmin*100),int(londeg),int(lonmin*100),
                    int(depth[iev]*100),int(2*nsta),
                    int(np.amax(az)),int(np.amin(dist))))

        """
        ###########################
        # Event Error Info String #
        ###########################

        errinfo:  'RMSTAZRDPPERRAZRDPPERRCDA###PERR##NUMHERRVERRFMS'
        This part of the header string includes:
        49      4       F4.2        RMS travel time residual.
        53      3       F3.0        Azimuth of largest principal error (deg E of N).
        56      2       F2.0        Dip of largest principal error (deg).
        58      4       F4.2        Size of largest principal error (km).
        62      3       F3.0        Azimuth of intermediate principal error.
        65      2       F2.0        Dip of intermediate principal error.
        67      4       F4.2        Size of intermediate principal error (km).
        71      3       F3.2        Coda duration magnitude from NCSN stations. *
        74      3       A3          Event location remark. (See table 7 below).
        77      4       F4.2        Size of smallest principal error (km).
        81      2       2A1         Auxiliary remarks (See note below).
        83      3       I3          Number of S times with weights greater than 0.1.
        86      4       F4.2        Horizontal error (km).
        90      4       F4.2        Vertical error (km).
        94      3       I3          Number of P first motions.
                                    (Not noted for this example so 0)
        ###########################
        """
        #errinfo = ('   0  0 0   0  0 0   0  0      0  %3i%4i%4i%3i' %
        #           (int(nsta),int(herr[iev]*100),int(verr[iev]*100),int(0)))
        errinfo = ('                                  %3i%4i%4i%3i' %
                   (int(nsta),int(herr[iev]*100),int(verr[iev]*100),int(0)))

        ##########################
        # KB NOTE 28 June 2021
        # Since my synthetic problem has no associated mag information I have removed all of this
        # Make sure not to mess with space.  If something needs to be re-added the spacing issue
        # needs to be of the utmost concern.
        # This fix made the whole thig work so I think it's ok.
        #########################


        """
        ####################
        # Data Info String #
        ####################

        Datainfo: 'SMGWDMGWDSADMG########NUM##MAGMGW#AMGAMWEVENTIDI10#PMGPMGW#ACDACDW##'
        This part of the header string includes:

        97      4       F4.1        Total of NCSN S-amplitude mag weights ~number of readings.*
        101     4       F4.1        Total of NCSN duration mag weights ~number of readings. *
        105     3       F3.2        Median-absolute-difference of NCSN S-amp magnitudes.
        108     3       F3.2        Median-absolute-difference of NCSN duration magnitudes.
        111     3       A3          3-letter code of crust and delay model. (See table 8 below).
        114     1       A1          Last authority for earthquake N=NCSC (USGS), B=UC Berkeley.
                                    (A T in this column is meaningless)
        115     1       A1          Most common P & S data source code. (See table 1 below).
        116     1       A1          Most common duration data source code. (See cols. 71-73)
        117     1       A1          Most common amplitude data source code.
        118     1       A1          Coda duration magnitude type code
        119     3       I3          Number of valid P & S readings (assigned weight > 0).
        122     1       A1          S-amplitude magnitude type code
        123     1       A1          "External" magnitude label or type code. Typically L for ML
                                    or W for MW. This information is not computed by Hypoinverse,
                                    but passed along, as computed by UCB.
        124     3       F3.2        "External" magnitude.
        127     3       F3.1        Total of "external" magnitude weights (~ number of readings).
        130     1       A1          Alternate amplitude magnitude label or type code (i.e. L for
                                    ML calculated by Hypoinverse from Wood Anderson amplitudes).
        131     3       F3.2        Alternate amplitude magnitude.
        134     3       F3.1        Total of the alternate amplitude mag weights ~no. of
                                    readings.
        137     10      I10         Event identification number
        147     1       A1          Preferred magnitude label code chosen from those available.
        148     3       F3.2        Preferred magnitude, chosen by the Hypoinverse PRE command.
        151     4       F4.1        Total of the preferred mag weights (~ number of readings). *
        155     1       A1          Alternate coda duration magnitude label or type code (i.e.Z).
        156     3       F3.2        Alternate coda duration magnitude.
        159     4       F4.1        Total of the alternate coda duration magnitude weights. *
        163     1       A1          QDDS version number of information. Starts at 0 for quick
                                    look reports. Incremented by one each time new information is
                                    added or revised: from quick location, final earthworm
                                    location with MD, ML added, etc.
        164     1       A1          “Origin instance” version number, distinguishes between
                                    different origins (hypocenters). It starts with 'a' ('0' for
                                    quick-look reports) and runs through the alphabet. When
                                    Berkeley has a final magnitude for each origin, the character
                                    is promoted to upper-case. 
        ####################
        """
        #datainfo = ('   0   0  0  0   TH   %3i L%3i%3i   0  0%10i                  \n' %
        #            (int(2*nsta),int(mag[iev]*100),int(nsta*10),cuspid[iev]))
        datainfo = ('                 TH   %3i L%3i%3i       %10i                  \n' %
                    (int(2*nsta),int(mag[iev]*100),int(nsta*10),cuspid[iev]))

        #########
        # KB NOTE 28 June 2021
        # At this point, I removed all the tailing information after the event ID in this line.  
        # That seemed to get rid of the line 199 hypphs reading problem that I was experiencing.
        # May need to be re-added depending on future synthetic problems.
        #########

        phafile.write(locinfo+errinfo+datainfo)

        """
        #####################
        # End Header String #
        #####################
        """

        """
        ######################
        # Start Phase String #
        ######################
        Then write each phase line to file:
           1) Phase lines with data for each individual station
           2) IF IN SHADOW FORMAT: Followed by coda duration parameters and seismogram recovery info
        The file format used below is from page 6 of the shadow document from the NCEDC
        """
        for ista in range(nsta):

            # Since times in the phase file are absolute times
            # Add the tt and noise togther
            pdelta = timedelta(seconds=p_time[iev,ista])
            sdelta = timedelta(seconds=p_time[iev,ista+nsta])

            # Add this to the event origin time
            dt = datetime.combine(dates[iev],times[iev])
            ptime = dt + pdelta
            stime = dt + sdelta

            # Now save the updated date and time information to
            # new variables
            date = ptime.date()
            ptime = ptime.time()
            stime = stime.time()

            # Format the time correctly
            pseconds = ptime.second + ptime.microsecond/1000000
            sseconds = stime.second + stime.microsecond/1000000

            """
            ###########################
            # Station Information Str #
            ###########################
            """
            stainfo = ('%5sSY  EHZ ' % (s_lab[ista].ljust(5,' ')))

            """
            ##############
            # P Phs Info #
            ##############
            """
            pinfo = (' P 0%4i%02i%02i%02i%02i%5i   0  0' % (date.year,date.month,date.day,
                                                            ptime.hour,ptime.minute,
                                                            int(pseconds*100)))

            """
            ##############
            # S Phs Info #
            ##############
            """
            sinfo = ('%5i S 0    0' % (int(sseconds*100)))

            """
            ##################
            # Amplitude Info #
            ##################
            """
            ampinfo = ('      0 0  0   0   0%4i%3i00  0     %3i  0  0   0   0H  -- 0EHZXX\n' % 
                       (int(dist[ista,iev]*10),int(ang[ista,iev]),int(az[ista,iev])))

            """
            #######################
            # Write Phase to File #
            #######################
            """
            phafile.write(stainfo+pinfo+sinfo+ampinfo)

        """
        #########################
        # Event Terminator Line #
        #########################
           1) Terminator ID number and optional trial hypocenter
           2) A shadow (usually blank)
        """
        phafile.write('%s%10i \n' % (strterm,cuspid[iev]))

    phafile.close() # Close file and return empty

    return None



def hypfile(nev,mod_nl,mod_ratio,mod_top,mod_v,mod_nl_reloc,mod_ratio_reloc,
            mod_top_reloc,mod_v_reloc):
    """
    This function write the hypoinverse-2000 control file or run file.
    This file contains information on the running.
    ###########
    PARAMETERS:
    nev (int) --- Number of events
    mod_nl (int) --- Number of layers in the raytracing velocity model
    mod_ratio (float) --- VPVS ratio in the raytracing velocity model
    mod_top (float array) --- Depth to top of layers in raytracing velocity model
    mod_v (float array) --- P velocity of layers in raytracing velocity model
    mod_nl_reloc (int) --- Number of layers in relocation velocity model
    mod_ratio_reloc (float) --- VPVS ratio in relocation velocity model
    mod_top_reloc (float array) --- Depth to top of layers in relocation velocity model
    mod_v_reloc (float array) --- P velocity of layers in relocation velocity model
    ###########
    """
    
    # First open file
    hypfile = open('synth_hypoINV/synth_run.hyp','w')

    # File header in comments (just explains the vmodel formats)
    # Helps keep files organized
    hypfile.write('* Simple Synthetic Model Location of %i Events \n' % nev)
    hypfile.write('* Using a %i layer vmodel for raytracing with \n' % mod_nl)
    for ilay in range(mod_nl):
        hypfile.write('* layer %i for raytracing: %2.2f %2.2f \n' % (ilay+1,mod_top[ilay],mod_v[ilay]))
    hypfile.write('* And a %i layer relocation vmodel with \n' % mod_nl)
    for ilay in range(mod_nl_reloc):
        hypfile.write('* layer %i for relocation: %2.2f %2.f \n' % (ilay+1,mod_top_reloc[ilay],mod_v_reloc[ilay]))

    # Hypoinverse-2000 Run Settings
    # The only thing that changes is the vpvs ratio
    hypfile.write('* Overall Run Settings \n')
    hypfile.write('LET 5 0 0 0 0 \n')       # Only need to match station code (since network and 
                                            # channel codes are not necessary)
    hypfile.write('POS %2.2f \n' % mod_ratio_reloc)     # Set P and S Ratio
    hypfile.write('ERR %f \n' % 0.01)       # Assumed timing errors
    hypfile.write('RMS 4 .16 1.5 3 \n')     # Default residual downweighting
    hypfile.write('NET 0 \n')               # No defined network
    hypfile.write('DI1 100 50 1 3 \n')      # Default initial distance downweighting (No 
                                            # initial weighting)
    hypfile.write('DIS 4 50 1 3 \n')        # Default distance downweighting (NCSN Default)
    hypfile.write('WET 1.0 0.75 0.5 0.25 \n')           # Default traditional weighting factors
    # Then Output Format
    hypfile.write('ERF T \n')               # Output error messages to terminal
    hypfile.write('TOP F \n')               # No page ejects for event ouputs
    hypfile.write('REP F F \n')             # Turns off terminal window outputs
                                            # **** KB Note: I've been changing this back and forth
                                            #               Sometimes it's good to see but can be a
                                            #               lot if you're running a large amount of 
                                            #               
    hypfile.write('KPR 0 \n')               # Print out only the final locations
    hypfile.write('H71 2 1 3 \n')           # Set file formats: 2 for output hypo71, 1 for 
                                            # input hypoinverse, 3 for 4 digit years in station 
                                            # file

    # Set file locations for:
    # Input Files
    hypfile.write('* Input File Formats \n')
    # Then Station Data
    hypfile.write('STA %s \n' % 'synth_stat.sta')  # Station Input File
    # Then Crustal Model
    hypfile.write('CRH 1 %s \n' % 'synth_P.crh')   # P Vel Model set to 1
    ####################
    # THESE LINES NOT NECESSARY IF THERE'S ONLY ONE VPVS RATIO
    #CRH 2 'synth_S.crh'     # S Vel Model set to 2
    #SAL 1 2                 # Link P and S Velocity Models Together
    ####################
    # Then Default Phase File
    hypfile.write('PHS %s \n' % 'synth_pha.phs')   # Identify Phase File
    hypfile.write('FIL \n')                        # Determine Phase Format Automatically
    
    # Outputs and Run Files
    hypfile.write('* Ouput Files and Run The Locations\n')
    # Then Default Ouput Files
    hypfile.write('PRT %s \n' % 'synth.prt')       # Print file
    hypfile.write('SUM %s \n' % 'synth.sum')       # Event Summary File
    # Then Loc
    hypfile.write('LOC\n')                         # Final command, runs relocations
    hypfile.close() # Close and return

    return None


def runbash(hypoinvwd):
    """
    This function creates and runs the bash file for hypoinverse-2000 from inside python
    ###########
    PARAMETERS:
    hypoinvwd (str) --- Path to working directory for hypoinverse-2000
    ###########
    """

    """
    #####################
    # Build Bash Script #
    #####################
    """
    bashfile = open('synth_hypoINV/hypoinv.sh','w')

    # Record Current Directory
    bashfile.write('#!/bin/bash\n')
    bashfile.write('cwd=$(pwd)\n')
    bashfile.write('pwd\n')
    # Change to hypoinv directory
    bashfile.write('cd synth_hypoINV\n')
    bashfile.write('pwd\n')
    # Run hypoinverse
    bashfile.write('OBJPATH="%s"\n' % hypoinvwd)
    bashfile.write('PROCESS=$OBJPATH./hyp1.40\n')
    bashfile.write('echo "@synth_run.hyp" | $PROCESS\n')
    #bashfile.write('echo "ctrl-C"\n')
    # Hypoinverse only closed by ctrl-C or ctrl-Z kill switch
    # This kills process from inside bash
    #bashfile.write('PID=$!\n')
    #bashfile.write('kill -INT $PID\n')

    # Change Back to Current Directory
    bashfile.write('cd $cwd\n')
    bashfile.write('pwd\n')
    # Close Bashfile
    bashfile.close()

    """
    ############################
    # Make bashfile executable #
    ############################
    """
    make_executable('synth_hypoINV/hypoinv.sh')

    """
    ###############################
    # Run hypoinverse from python #
    ###############################
    """
    datet = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    print('Run Bash Hypoinverse %s' % datet)
    subprocess.call('synth_hypoINV/hypoinv.sh')
    datet = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    print('Hypoinverse Finished %s' % datet)

    return None # Return empty


def readcat(dates,times,lat,lon,depth,mag,
            res,herr,verr,cuspid):
    """
    Reads in and updates the event arrays with the hypoinverse catalog locations
    :::
    PARAMETERS:
    dates (object array) --- Numpy object array containing the event datetime.date objects
    times (object array) --- Numpy object array containing the event datetime.time objects
    lat (float array) --- Numpy array containing the event latitudes
    lon (float array) --- Numpy array containing the event longitudes
    depth (float array) --- Numpy array containing the event depths
    mag (float array) --- Numpy array containing the event magnitudes
    res (float array) --- Numpy array containing the event residuals
    herr (float array) --- Numpy array containing the event horizontal location errors
    verr (float array) --- Numpy array containing the event vertical location errors
    cuspid (int array) --- Numpy array containing the event IDs
    :::
    """

    # Open and read in hypoinverse-2000 output
    evcatfile = open('synth_hypoINV/synth.sum','r')
    events = evcatfile.readlines()

    for eindx,event in enumerate(events):
        event = event.split(' ')
        event = list(filter(None,event))

        # Convert event date/times to correct formats
        year = int(event[0][0:4])
        mon = int(event[0][4:6])
        day = int(event[0][6::])
        hr = int(event[1][0:2])
        mn = int(event[1][2::])
        sec = int(float(event[2]))
        microsecond = int((float(event[2])-int(float(event[2])))*1000000)

        dt = datetime(year,mon,day,hr,mn,sec,microsecond)
        #datetimes[eindx] = dt
        dates[eindx] = dt.date()
        times[eindx] = dt.time() 

        # Convert event latitudes and longitudes back into decimal
        # degrees from degree minutes
        latdeg = int(event[3])
        latmin = float(event[4])
        lat[eindx] = dms2decdeg(latdeg,latmin)
        longitude = event[5].split('W')
        londeg = int(longitude[0])
        lonmin = float(longitude[1])
        lon[eindx] = dms2decdeg(londeg,lonmin)

        # Record the other variables as well
        depth[eindx] = float(event[6])
        mag[eindx] = float(event[7][0:4])
        # Typically the res and location errors will be zero before the
        # catalog location stage.  These values are now recorded.
        res[eindx] = float(event[10])
        herr[eindx] = float(event[11])
        verr[eindx] = float(event[12])

        cuspid[eindx] = int(event[-2])

    evcatfile.close()
    # Close and return event information

    return dates,times,lat,lon,depth,mag, \
           res,herr,verr,cuspid


"""
Hypoinverse Wrapper Function
"""
def hypoinverse(p_time,dist,az,ang,nsta,s_lab,s_lat,s_lon,
                rt_nl,rt_ratio,rt_top,rt_v,
                reloc_nl,reloc_ratio,reloc_top,reloc_v,
                nev,cuspid,lat,lon,depth,
                dates,times,mag,herr,verr,res):
    """
    Hypoinverse catalog building step
    :::
    Parameters:
    p_time (2D float array (nev,2*nsta)) --- Phase traveltimes
    dist (2D float array (nsta,nev)) --- Event-station sep. dist. in km
    az (2D float array (nsta,nev)) --- Event-station azimuths
    ang (2D float array (nsta,nev)) --- Take-off angles
    nsta (int) --- No. of stations
    s_lab (object array (nsta)) --- Station codes
    s_lat (float array (nsta)) --- Station latitudes
    s_lon (float array (nsta)) --- Station longitudes
    rt_nl (int) --- Raytracing vmodel no. of layers
    rt_ratio (float) --- Raytracing vmodel VPVS ratio
    rt_top (float array) --- Raytracing vmodel layer depths
    rt_v (float array) --- Raytracing vmodel layer P velocities
    reloc_nl (int) --- Relocation vmodel no. of layers
    reloc_ratio (float) --- Relocation vmodel VPVS ratio
    reloc_top (float array) --- Relocation vmodel layer depths
    reloc_v (float array) --- Relocation vmodel layer P velocities
    nev (int) --- No. of events
    cuspid (int array) --- Event IDs
    lat (float array) --- Event latitudes
    lon (float array) --- Event longitudes
    depth (float array) --- Event depths
    dates (object array) --- Event datetime.date objects
    times (object array) --- Event datetime.time objects
    mag (float array) --- Event magnitudes
    herr (float array) --- Horizontal event location errors
    verr (float array) --- Vertical event location errors
    res (float array) --- Event traveltime residuals
    :::
    Returns:
    dates (object array) --- Event datetime.date objects
    times (object array) --- Event datetime.time objects
    lat (float array) --- Event latitudes
    lon (float array) --- Event longitudes
    depth (float array) --- Event depths
    mag (float array) --- Event magnitudes
    res (float array) --- RMS of event traveltime residuals
    herr (float array) --- Horizontal event errors
    verr (float array) --- Vertical event errors
    cuspid (int array) --- Event ID arrays
    :::
    """

    """
    Write files for hypoINV using the vmodel for relocation
    """
    # Use this relocation vmodel to write the vmodel files for hypoinv
    velfile(reloc_nl,reloc_v,reloc_top)
    # Write station file
    stafile(nsta,s_lab,s_lat,s_lon)
    # Write phase file
    phafile(p_time,dist,az,ang,
            nsta,s_lab,s_lat,s_lon,# Set first iteration start to 0
            nev,cuspid,lat,lon,depth,
            dates,times,mag,herr,verr,res)
    # Write run file
    hypfile(nev,rt_nl,rt_ratio,rt_top,rt_v,reloc_nl,reloc_ratio,reloc_top,reloc_v)

    """
    Run hypoinverse 
    ---
    Makes bash script then executes it
    """
    hypoinvwd = hypoinvbin
    runbash(hypoinvwd)

    """
    Update Event Arrays
    ---
    Update lat,lon,depth, and maybe dates/times from 
    the output of hypoinverse
    Read in .sum file from hypoinverse and update event arrays
    """
    [date,times,lat,lon,depth,mag,
     res,herr,verr,cuspid] = readcat(dates,times,lat,lon,depth,mag,
                                     res,herr,verr,cuspid)

    """
    Return updated event arrays
    """
    return [dates,times,lat,lon,depth,mag,
            res,herr,verr,cuspid]


"""
To File and Not to Arrays
"""
def hypoinverse_wrapper(log,fileout=0,hinput='hypoDD.inp',sinput='synth.inp',data=[]):
    """
    Prep for hypoinverse step
    :::
    Parameters:
    log (file object) --- Log file
    hinput (string) --- HypoDD.inp file location
    sinput (string) --- Synth.inp file location
    :::
    """

    """
    Read in synth.inp
    """
    [statfile,nev,depth,latmid,lonmid,
     stepsize,latstep,lonstep,
     rt_nl,rt_ratio,rt_top,rt_v] = synth_input(log,sinput)

    """
    Read in hypoDD.inp
    """
    [fn_cc,fn_ct,fn_sta,fn_eve,
     fn_loc,fn_reloc,fn_res,fn_stares,fn_srcpar,
     idata,iphase,minobs_cc,minobs_ct,
     amaxres_cross,amaxres_net,amaxdcc,amaxdct,
     maxdist,awt_ccp,awt_ccs,awt_ctp,awt_cts,
     adamp,istart,maxiter,isolv,
     niter,aiter,
     reloc_nl,reloc_ratio,reloc_v,reloc_top,
     iclust,ncusp,icusp] = hypoDD_input(log,hinput,fileout=fileout)


    if fileout==0:
        """
        Read in events and stations
        """
        [nsta,s_lab,s_lat,s_lon] = readstat(log,statfile)

        """
        Read in phases
        """
        [nev,lat,lon,depth,cuspid,dates,times,mag,herr,verr,res,
         npha,nobs_ct,p_pha,p_sta,p_time,p_wghtr] = readphase(log,'phase.dat',0.,minobs_ct,nsta)
        [ttp,tts,xp,yp,zp,
         dist,az,ang] = partials(nev,cuspid,lat,lon,depth, 
                                 nsta,s_lab,s_lat,s_lon, 
                                 rt_nl,rt_ratio,rt_v,rt_top,return_all=True)
    elif fileout==1:
        """
        Unpack data arrays
        """
        [nsta,s_lab,s_lat,s_lon,
         nev,lat,lon,depth,cuspid,dates,times,mag,herr,verr,res,
         p_time,dist,az,ang] = data

    """
    Call hypoinverse wrapper
    """
    [dates,times,lat,lon,depth,mag,
     res,herr,verr,cuspid] = hypoinverse(p_time,dist,az,ang,
                                         nsta,s_lab,s_lat,s_lon,
                                         rt_nl,rt_ratio,rt_top,rt_v,
                                         reloc_nl,reloc_ratio,reloc_top,reloc_v,
                                         nev,cuspid,lat,lon,depth,
                                         dates,times,mag,herr,verr,res)

    if fileout==0:
        """
        Update event file (event.dat) if fileout
        """
        evfile = open(fn_eve,'w')
        for i in range(len(cuspid)):
            evdate = '%04i%02i%02i' % (dates[i].year,dates[i].month,dates[i].day)
            evtime = '%02i%02i%02i%02i' % (times[i].hour,times[i].minute,times[i].second,int(times[i].microsecond/10000))
            evfile.write('%s %s %5.6f %5.6f %5.6f %1.2f %6.2f %6.2f %6.2f %8i \n' %
                         (evdate,evtime,lat[i],lon[i],depth[i],mag[i],herr[i],verr[i],res[i],cuspid[i]))
        evfile.close()

        """
        Update phase file
        """
        phfile = open('phase.dat','r')
        phases = phfile.readlines()
        phfile.close()
        phfile = open('phase.dat','w')
        for ph in phases:
            if ph[0]=='#':
                ph = ph.split(' ')
                ph = list(filter(None,ph))
                i = int(ph[-2])-1
                phfile.write('# %4i %2i %2i %2i %2i %2.2f %5.6f %5.6f %5.6f %1.2f %6.2f %6.2f %6.2f %8i \n' %
                             (dates[i].year,dates[i].month,dates[i].day,
                              times[i].hour,times[i].minute,(times[i].second + (times[i].microsecond/1000000)),
                              lat[i],lon[i],depth[i],mag[i],herr[i],verr[i],res[i],cuspid[i]))
            else:
                phfile.write(ph)
        phfile.close()
        
        return None
    elif fileout==1: 
        """
        Else return updated event arrays
        """   
        return [dates,times,lat,lon,depth,mag,res,herr,verr,cuspid]

    raise Exception('End of function hypoinverse wrapper.')
