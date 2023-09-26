#!/usr/bin/env python3
import numpy as np
from datetime import datetime
import sys
import os
import shutil

# Import run functions
from runOptions.runFunctions import run_input,readinfiles,readinnoise,statsOut
from runOptions.bootstrap import bootstrap 
from runOptions.plotting import resplot,locplots,evpairsplot

# Import ph2dt functions
from ph2dt.ph2dt_files import ph2dt_input
from ph2dt.ph2dt import ph2dt

# Import hypoDD functions
from hypoDD.hypoDD_files import hypoDD_input
from hypoDD.hypoDD import hypoDD

# Import utility functions
from utility.synthetics.synthModel import synth_generate
from utility.synthetics.noise import generate_noise
from utility.synthetics.hypoinv import hypoinverse_wrapper


"""

:::
From relocDD-py VERSION 1.0
:::
relocDD-py includes multiple method implementations of double-difference (DD)
relocation including event-pair (hypoDD), station-pair, and double-pair DD methods
:::
The main code workflow and architecture is based on the Implementation of HypoDD Version 1.3 -
Original Author: Felix Waldhauser, felixw@ldeo.columbia.edu
Original Code in Fortran can be Found: https://www.ldeo.columbia.edu/~felixw/hypoDD.html
See hypoDD user guide for a description of many of the parameters seen in this code.
Citations required if you use this code:
    Waldhauser F. and W.L. Ellsworth, A double-difference earthquake location algorithm: 
        Method and application to the northern Hayward fault, Bull. Seism. Soc. Am., 90, 
        1353-1368, 2000.
    Waldhauser, F., HypoDD: A computer program to compute double-difference earthquake 
        locations, USGS Open File Rep., 01-113, 2001.
:::
Citations for the additional DD methods (station and double-pair) are:
    Zhang, H., Nadeau, R.M., and M.N., Locating nonvolcanic tremors beneath the San 
        Andreas Fault using a station-pair double-difference location method, 
        Geophys. Res. Lett., 37, L13304, 2010.
    Guo, H. and H. Zhang, Development of double-pair double difference earthquake 
        location algorithm for improving earthquake locations, Geophys. J. Int. 208, 
        333-348, 2017.
Additional citations for some subroutines can be found in those specific functions.
:::
Author Python Version of RelocDD: K. Biegel, katherine.biegel@ucalgary.ca
Python Version and Updates can be Found: 
Updated User Manual: 
The updated user manual is specific to the relocDD-py base implementation.
Citations for relocDD-py specifically:

:::
"""


"""
---
Run.py serves as an example run script for relocDD-py.
---
This script can be run in current form or used to implement your own run structure.
---
"""
def main(finput,yesph2dt=False,yeshypoDD=False,yesplot=False,yesboot=False):
    """
    :::
    Main run function.  
    :::
    Depending on run.inp variables
    this will call everything from the classic hypoDD run to synthetic 
    modelling to bootstrapping.
    :::
    There are two ways to run main() either by:
        python main.py [/path/to/run.inp] [yesph2dt] [yeshypoDD] [yesplot] [yesboot]
    Or by running the following in another script:
        from run import main
        main(finput,yesph2dt,yeshypoDD,yesplot,yesboot)
    :::
    Parameters:
        finput (str)        --- Path to run input file
        yesph2dt (bool)     --- Switch for diff. time catalog building
        yeshypoDD (bool)    --- Switch for relocation
        yesplot (bool)      --- Switch for plotting outputs
        yesboot (bool)      --- Switch for bootstrapping
    :::
    Returns:
        NONE
    :::
    """


    """
    :::
    Read in run input file
    :::
    Path to run.inp must be specified in the 
    function call statement
    :::
    """
    [inputfol,datfol,outfol,reloctype,fileout,makedata,hypoinv,
     noiseswitch,noisediff,stdcc,stdct,nboot,nplot] = run_input(finput)

    """
    :::
    Declare Folders for Outputs if Necessary
    :::
    """
    # Method specific outfolder
    os.makedirs(outfol,exist_ok='True')
    if reloctype==1:
        outfol = os.path.join(outfol,'EDD')
        if os.path.isdir(outfol):
            shutil.rmtree(outfol)
        os.makedirs(outfol)
    elif reloctype==2:
        outfol = os.path.join(outfol,'SDD')
        if os.path.isdir(outfol):
            shutil.rmtree(outfol)
        os.makedirs(outfol)
    elif reloctype==3:
        outfol = os.path.join(outfol,'DDD')
        if os.path.isdir(outfol):
            shutil.rmtree(outfol)
        os.makedirs(outfol)
    # Output file folders
    tradout = os.path.join(outfol,'tradouts')
    os.makedirs(tradout,exist_ok='True')

    """
    :::
    Open general log file
    :::
    """
    if yesph2dt and not yeshypoDD:
        log = open(os.path.join(outfol,'ph2dt.log'),'w')
    elif yeshypoDD and not yesph2dt:
        log = open(os.path.join(outfol,'hypoDD.log'),'w')
    elif yesph2dt and yeshypoDD:
        log = open(os.path.join(outfol,'combo.log'),'w')
    else:
        log = open(os.path.join(outfol,'out.log'),'w')
    datet = str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    log.write('\n\n(main) Beginning hypoDD_py run at %s. \n\n' % datet)

    """
    :::
    Read in additional input files based on inputfol
    :::
    """
    # Set input file path names
    pinput = os.path.join(inputfol,'ph2dt.inp')
    hinput = os.path.join(inputfol,'hypoDD.inp')

    # Read in ph2dt input file
    pinputs = ph2dt_input(log,pinput,dd_version=reloctype)
    log.write('\nPh2dt.inp read into memory.\n')
    # Update path strings
    pinputs[0] = os.path.join(datfol,pinputs[0])
    pinputs[1] = os.path.join(datfol,pinputs[1])
    log.write('\nData file path updated for station and phase files.\n')

    # Read in hypoDD input files
    hinputs = hypoDD_input(log,hinput)
    log.write('\nPh2dt.inp read into memory.\n')
    # Update data path strings
    hinputs[0] = os.path.join(datfol,hinputs[0])
    hinputs[1] = os.path.join(datfol,hinputs[1])
    hinputs[2] = os.path.join(datfol,hinputs[2])
    hinputs[3] = os.path.join(datfol,hinputs[3])
    log.write('\nData file path updated for cc, ct, event, and station files.\n')
    # Update hypoDD output strings 
    hinputs[4] = os.path.join(tradout,hinputs[4])
    hinputs[5] = os.path.join(tradout,hinputs[5])
    hinputs[6] = os.path.join(tradout,hinputs[6])
    hinputs[7] = os.path.join(tradout,hinputs[7])  
    hinputs[8] = os.path.join(tradout,hinputs[8])
    log.write('\nOutput file path updated for hypoDD loc, reloc, sta, res, and src files.\n')

    """
    :::
    Set random seed
    :::
    Set once at the beginning
    Keeps noise consistent throughout trials
    Used in bootstrapping and also synthetic 
    cases
    :::
    """
    np.random.seed(0)
    log.write('\n(main) Set noise seed to 0. \n')

    
    """
    :::
    Here is where the main run loop occurs.
    :::
    """

    if not yesboot: # Single-run version (not bootstrapping)

        """
        :::
        Generate Synthetic Data If Needed
        :::
        """
        synthmod = []
        if makedata==1:
            # Synthetic Input File
            sinput = os.path.join(inputfol,'synth.inp')

            # Generate Synthetic Data
            log.write('Build synthetic model.\n\n')
            if fileout==0:
                print('Building synthetic model...')
            synthmod = synth_generate(log,sinput,datfol,fileout)
            if fileout==0:
                print('Phase.dat file generated.')

        """
        :::
        Run ph2dt
        :::
        """
        catalog = []
        if yesph2dt:
            log.write('\n\nRun ph2dt.\n\n')
            if fileout==0:
                print('\n\nRun ph2dt.')
            
            tempstart = datetime.now()
            catalog = ph2dt(log,pinputs,makedata=makedata,datfol=datfol,outfol=outfol,
                            reloctype=reloctype,fileout=fileout,idata=hinputs[9],
                            datavals=synthmod)
            
            if fileout==0:
                print('\n\nPh2dt completed.')
                print('Time to run ph2dt: ',datetime.now()-tempstart)
            log.write('\n\nPh2dt completed.')
            log.write('Time to run ph2dt: %i\n' % (datetime.now()-tempstart).seconds)

        # """
        # ---
        # Updates for synthetic models
        # ---
        # """
        # if makedata==1:
        #     """
        #     ---
        #     Add noise to datafiles
        #     ---
        #     """
        #     if noiseswitch==1:
        #         log.write('Add Noise.\n\n')
        #         if fileout==0:
        #             print('Add Noise.')
        #         noise = generate_noise(0,reloctype,noisediff,stdcc,stdct,data=['dt.cc','dt.ct'])
        #         if fileout==0:
        #             print('Add Noise completed.')

        #     """
        #     ---
        #     Run hypoinverse catalog step
        #     ---
        #     """
        #     if hypoinv==1:
        #         log.write('Updated Event Locations using HypoINVERSE-2000.\n\n')
        #         if fileout==0:
        #             print('Updated Event Locations using HypoINVERSE-2000.')
        #         cat = hypoinverse_wrapper(log,0,hinput,'synth.inp',data=[])
        #         if fileout==0:
        #             print('Catalog Locations completed.')

        """
        :::
        Run hypoDD (relocation)
        :::
        """
        if yeshypoDD:
            log.write('\n\nRun relocation.\n\n')
            if fileout==0:
                print('\n\nRun relocation.')

            tempstart = datetime.now()
            outputs = hypoDD(log,hinputs,outfol=outfol,reloctype=reloctype,fileout=0,
                             hypoDDdata=catalog)

            log.write('\n\nRelocation completed.')
            log.write('Time to run relocation: %i\n' % (datetime.now()-tempstart).seconds)

    elif yesboot:
        """
        :::
        Run bootstrapping
        :::
        Bootstrapping requires fileout==1 for limited 
        file/terminal input/outputs.
        :::
        """
        log.write('\n\nStarting bootstrapping... \n\n')
        print('\n\nStarting bootstrapping... \n')

        plotting=False
        plotstring=''

        """
        Run bootstrapping loop
        """
        tempstart = datetime.now()
        bootstrap(log,reloctype,nboot,pinputs,hinputs,plotting,plotstring,nplot)
        
        print('\nTime to run bootstrap: ',datetime.now()-tempstart)
        log.write('\nTime to bootstrap: %i\n' % (datetime.now()-tempstart).seconds)

    # niter = 8
    # nclust = 4
    # """
    # If plotting:
    # """
    # if plotting==1:
    #     """"
    #     Final plots for relocations
    #     """

    #     """
    #     Read in or organize data
    #     """
    #     if fileout==0:
    #         tempstart = datetime.now()
            
    #         """
    #         Read in from files
    #         """
    #         if reloctype==2:
    #             ncust=1
    #         [x,cal1,cal4,dtdt1,dtdt4,del1,del4,
    #          abs0,abs1,abs4,tru1] = readinfiles(noiseswitch,nclust,niter)
            
    #         print('\nTime to run readinfiles (run.py:main): ',datetime.now()-tempstart)
    #         log.write('\nTime to readinfiles (run.py:main): %i\n' % (datetime.now()-tempstart).seconds)

    #         ndt = len(dtdt1)
    #         ncc = int(ndt/2)

    #     # if fileout==2:
    #     #     tempstart = datetime.now()
    #     #     x = readinnoise(noiseswitch)
    #     #     print('\nTime to run readinnoise (run.py:main): ',datetime.now()-tempstart)
    #     #     log.write('\nTime to readinnoise (run.py:main): %i\n' % (datetime.now()-tempstart).seconds)

    #     #     cal1 = calstart
    #     #     cal4 = calend
    #     #     dtdt1 = dtdtstart
    #     #     dtdt4 = dtdtend

    #     #     del1 = np.transpose(locdel[:,:,0])
    #     #     del2 = np.transpose(locdel[:,:,1])
    #     #     del3 = np.transpose(locdel[:,:,2])
    #     #     del4 = np.transpose(locdel[:,:,3])

    #     #     abs0 = np.transpose(locabs[:,:,0])
    #     #     abs1 = np.transpose(locabs[:,:,1])
    #     #     abs2 = np.transpose(locabs[:,:,2])
    #     #     abs3 = np.transpose(locabs[:,:,3])
    #     #     abs4 = np.transpose(locabs[:,:,4])

    #     #     tru1 = np.transpose(loctru)

    #     #     ndt = len(dtdt1)
    #     #     ncc = int(ndt/2)

    #     # elif fileout==1:
    #     #     """
    #     #     Transform Arrays
    #     #     """
    #     #     cal1 = calstart
    #     #     cal4 = calend
    #     #     dtdt1 = dtdtstart
    #     #     dtdt4 = dtdtend

    #     #     del1 = np.transpose(locdel[:,:,0])
    #     #     del2 = np.transpose(locdel[:,:,1])
    #     #     del3 = np.transpose(locdel[:,:,2])
    #     #     del4 = np.transpose(locdel[:,:,3])

    #     #     abs0 = np.transpose(locabs[:,:,0])
    #     #     abs1 = np.transpose(locabs[:,:,1])
    #     #     abs2 = np.transpose(locabs[:,:,2])
    #     #     abs3 = np.transpose(locabs[:,:,3])
    #     #     abs4 = np.transpose(locabs[:,:,4])

    #     #     tru1 = np.transpose(loctru)

    #     """
    #     Ouput Statistics:
    #     """
    #     tempstart = datetime.now()
    #     [ed_tru,ed_ini,ed_rel] = statsOut(makedata,noiseswitch,reloctype,
    #                                       noisediff,stdcc,stdct,
    #                                       x,ncc,ndt,dtdt1,dtdt4,cal1,cal4,
    #                                       tru1,abs0,abs4,hypoinv)
    #     print('\nTime to run statsOut (run.py:main): ',datetime.now()-tempstart)
    #     log.write('\nTime to statsOut (run.py:main): %i\n' % (datetime.now()-tempstart).seconds)
    #     ####################
    #     # Relocation Plots #
    #     ####################
    #     nev = (abs0.shape)[0]

    #     if not os.path.isdir(plotstring):
    #         os.mkdir(plotstring)

    #     """
    #     Residual Histogram Plot
    #     """
    #     if yesplot:
    #         tempstart = datetime.now()
    #         resplot(plotstring,noiseswitch,reloctype,
    #                 x,stdcc,stdct,dtdt1,dtdt4,cal1,cal4,noCC=True)
    #         print('\nTime to run resplot (run.py:main): ',datetime.now()-tempstart)
    #         log.write('\nTime to resplot (run.py:main): %i\n' % (datetime.now()-tempstart).seconds)

    #         """
    #         Absolute Location Plots
    #         """
    #         tempstart = datetime.now()
    #         locplots(plotstring,reloctype,nev,
    #                  tru1,abs0,abs1,abs4,
    #                  ed_tru,ed_ini,ed_rel,hypoinv)
    #         print('\nTime to run locplots (run.py:main): ',datetime.now()-tempstart)
    #         log.write('\nTime to locplots (run.py:main): %i\n' % (datetime.now()-tempstart).seconds)

    #         """
    #         Event Pairs
    #         """
    #         tempstart = datetime.now()
    #         #if (reloctype==1 or reloctype==3) and fileout==0:
    #         evpairsplot(nev,plotstring,reloctype,tru1,abs0,abs4,hypoinv)
    #         print('\nTime to run evpairsplot (run.py:main): ',datetime.now()-tempstart)
    #         log.write('\nTime to evpairsplot (run.py:main): %i\n' % (datetime.now()-tempstart).seconds)

    return None


tstart = datetime.now()
# Get yesph2dt switch
print('\nBeginnging main (run.py) at ',tstart,' ...')
if len(sys.argv)==4:
    fstr = str(sys.argv[1])
    yesph2dt = bool(int(sys.argv[2]))
    yeshypoDD = bool(int(sys.argv[3]))
    main(fstr,yesph2dt=yesph2dt,yeshypoDD=yeshypoDD)
elif len(sys.argv)==5:
    fstr = str(sys.argv[1])
    yesph2dt = bool(int(sys.argv[2]))
    yeshypoDD = bool(int(sys.argv[3]))
    yesplot = bool(int(sys.argv[4]))
    main(fstr,yesph2dt=yesph2dt,yeshypoDD=yeshypoDD,yesplot=yesplot)
elif len(sys.argv)==6:
    fstr = str(sys.argv[1])
    yesph2dt = bool(int(sys.argv[2]))
    yeshypoDD = bool(int(sys.argv[3]))
    yesplot = bool(int(sys.argv[4]))
    yesboot = int(sys.argv[5])
    main(fstr,yesph2dt=yesph2dt,yeshypoDD=yeshypoDD,yesplot=yesplot,yesboot=yesboot)
else:
    main('run.inp')
print('\n\nTime to run main (run.py): ',datetime.now()-tstart)

sys.exit()
