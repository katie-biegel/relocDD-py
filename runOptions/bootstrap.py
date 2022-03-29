#!/usr/bin/env python
import os
import sys
import numpy as np

# Import needed hypoDD-py functions
#from ph2dt.ph2dt import ph2dt
#from ph2dt_subfunc import ph2dt_prep
#from synth_model import generate_noise
#from hypoinv_wrapper import hypoinverse_catalog
from hypoDD.hypoDD import hypoDD


"""
This script contains the functions needed to run bootstrapping.
----
The built-in bootstrapping requires the limited file i/o run option.
----
This script contains the functions:
    savefiles:      Save files from bootstrapping iterations
    plotting:       Plotting for bootstrap iterations w/ plotting
    bootstrap:      Run bootstrapping
"""

def savefiles(bsfolder,iboot,nev,calstart,calend,dtdtstart,dtdtend,
              loctru,locdel,locabs,x):    
    """
    Save bootstrap outputs to file
    :::
    Parameters:
    bsfolder (str) --- Bootstrap output folder
    iboot (int) --- Bootstrap iteration
    nev (int) --- Number of events
    calstart[ndt] (float array) --- Calculated diff. times starting iteration
    calend[ndt] (float array) --- Calculated diff. times final iteration
    dtdtstart[ndt] (float array) --- Measured diff. times starting iteration
    dtdtend[ndt] (float array) --- Measured diff. times final iteration
    loctru[nev,4] (float array) --- True locations (x,y,z,t)
    locdel[nev,4,niter] (float array) --- Delta model parameters (x,y,z,t) per iteration
    locabs[nev,4,niter] (float array) --- Abs. model parameters (x,y,z,t) per iteration
    x[ndt](float array) --- Noise array
    :::
    Return:
    None
    :::
    """

    """
    New folder for each iteration
    """
    ibootstr = bsfolder + '/' + str(iboot+1)
    if not os.path.isdir(ibootstr):
        os.mkdir(ibootstr)

    """
    Save Files
    """
    # Save traveltimes residuals to file
    np.savetxt('%s/tot_calstart.txt' % ibootstr,calstart)
    np.savetxt('%s/tot_calend.txt' % ibootstr,calend)
    np.savetxt('%s/tot_dtdtstart.txt' % ibootstr,dtdtstart)
    np.savetxt('%s/tot_dtdtend.txt' % ibootstr,dtdtend)
    # Save true locations
    np.savetxt('%s/tot_loctru.txt' % ibootstr,loctru)
    # Save nosie
    np.savetxt('%s/noise.txt' % ibootstr,x)
    # Loop over and save location updates in event specific files
    # *** Necessary since there are multiple inversion iterations
    for iev in range(nev):
        np.savetxt('%s/tot_locdel_%i.txt' % (ibootstr,iev+1),locdel[:,iev,:])
        np.savetxt('%s/tot_locabs_%i.txt' % (ibootstr,iev+1),locabs[:,iev,:])

    return None


def plotting(plotstrin,cal1,cal4,dtdt1,dtdt4,locdel,locabs,loctru,nev,niter=4):
    """
    Generate and save plots for bootstrapping
    :::
    Parameters:
    plotstrin (str) --- Plotting folder string location
    cal1[ndt] (float array) --- Starting cal. dt times
    cal4[ndt] (float array) --- Last iteration cal. dt times
    dtdt1[ndt] (float array) --- Starting meas. dt times
    dtdt4[ndt] (float array) --- Last iteration meas. dt times
    locdel[4,nev,niter] (float array) --- Change in model parameter each iteration
    locabs[4,nev,niter] (float array) --- Abs. model parameters each iteration
    loctru[4,nev] (float array) --- True event locations
    niter (int) --- Number of relocation inversion iterations
    :::
    Returns:
    plotcount (int) ---- Set to 0, returns reset plotcounter 
    :::
    """

    """
    Begin Plotting
    """
    # Build plotting folder
    plotstr = bsfolder + '/' + plotstrin
    if not os.path.isdir(plotstr):
        os.mkdir(plotstr)

    # Plot 1: Map view absolute locations with iteration
    plt.figure()
    plt.title('Abs Locs in XY')
    for i in range(nev):
        plt.plot(loctru[0,i],loctru[1,i],'ok') # Need to code in plotting all iterations
        for j in range(niter):
            plt.plot(locabs[i,0],abs1[i,1],'xk')
            #plt.plot(locabs[i,0],abs2[i,1],'xr')
            #plt.plot(locabs[i,0],abs3[i,1],'xg')
            plt.plot(locabs[i,0],abs4[i,1],'ob')
    plt.axis('equal')
    plt.savefig('%s/fig_xyabs_%i.png' % (plotstr,iboot+1))

    # Plot 2: Depth view absolute locations with iteration
    plt.figure()
    plt.title('Abs locs in depth')
    for i in range(nev):
        #plt.plot(tru1[i,1],tru1[i,2],'ok')
        plt.plot(abs1[i,1],abs1[i,2],'xk')
        #plt.plot(abs2[i,1],abs2[i,2],'xr')
        #plt.plot(abs3[i,1],abs3[i,2],'xg')
        plt.plot(abs4[i,1],abs4[i,2],'ob')
    plt.gca().invert_yaxis()
    plt.axis('equal')
    plt.savefig('%s/fig_depabs_%i.png' % (plotstr,iboot+1))

    # Change to true locations map view
    plt.figure()
    plt.title('Subtract true locations in XY')
    for i in range(nev):
        plt.plot(abs1[i,0]-tru1[i,0],abs1[i,1]-tru1[i,1],'xk')
        #plt.plot(abs2[i,0]-tru1[i,0],abs2[i,1]-tru1[i,1],'xr')
        #plt.plot(abs3[i,0]-tru1[i,0],abs3[i,1]-tru1[i,1],'xg')
        plt.plot(abs4[i,0]-tru1[i,0],abs4[i,1]-tru1[i,1],'ob')
    plt.grid(True)
    plt.axis('equal')
    plt.savefig('%s/fig_xy_minustrue_%i.png' % (plotstr,iboot+1))

    # Change to depth locations depth view
    plt.figure()
    plt.title('Subtract true locations in depth')
    for i in range(nev):
        plt.plot(abs1[i,1]-tru1[i,1],abs1[i,2]-tru1[i,2],'xk')
        #plt.plot(abs2[i,1]-tru1[i,1],abs2[i,2]-tru1[i,2],'xr')
        #plt.plot(abs3[i,1]-tru1[i,1],abs3[i,2]-tru1[i,2],'xg')
        plt.plot(abs4[i,1]-tru1[i,1],abs4[i,2]-tru1[i,2],'ob')
    plt.grid(True)
    plt.axis('equal')
    plt.savefig('%s/fig_dep_minustrue_%i.png' % (plotstr,iboot+1))

    # Map view deltas x and y per iteration
    #plt.figure()
    #plt.title('Del Locs in XY')
    #for i in range(nev):
    #    plt.plot(del1[i,0],del1[i,1],'xk')
    #    plt.plot(del2[i,0],del2[i,1],'xr')
    #    plt.plot(del3[i,0],del3[i,1],'xg')
    #    plt.plot(del4[i,0],del4[i,1],'ob')
    #plt.axis('equal')
    #plt.savefig('%s/fig_xydel_%i.png' % (plotstr,iboot+1))

    # Depth view deltas y and z per iteration
    #plt.figure()
    #plt.title('Del locs in depth')
    #for i in range(nev):
    #    plt.plot(del1[i,1],del1[i,2],'xk')
    #    plt.plot(del2[i,1],del2[i,2],'xr')
    #    plt.plot(del3[i,1],del3[i,2],'xg')
    #    plt.plot(del4[i,1],del4[i,2],'ob')
    #plt.gca().invert_yaxis()
    #plt.axis('equal')
    #plt.savefig('%s/fig_depdel_%i.png' % (plotstr,iboot+1))

    # Historgram Plots side by side
    #res = dtdt4-cal4
    #res_start = dtdt1-cal1

    #NDAT=int((x.shape)[0]/2)
    # fig, axs = plt.subplots(nrows=1, ncols=2)
    # fig.suptitle('Noise vs Final Residuals')
    # axs[0].set_title('CC Data')
    # axs[0].hist(x[:NDAT], 100,color='k', density=False, histtype='step',alpha=0.75)
    # axs[0].hist(res[:NDAT], 100,color='r', density=False, histtype='step',alpha=0.75)
    # axs[1].set_title('CT Data')
    # axs[1].hist(x[NDAT:], 100,color='k', density=False, histtype='step',alpha=0.75)
    # axs[1].hist(res[NDAT:], 100,color='r', density=False, histtype='step',alpha=0.75)
    # plt.savefig('%s/fig_resnoise_hist_%i.png' % (plotstr,iboot+1))

    # plt.close('all')
    
    """
    Return reset plot count
    """
    plotcount=0
    return plotcount


def bootstrap(log,reloctype,nboot,pinput,hinput,iplot,plotstrin,nplot,imakedata=0,
              noiseswitch=0,noisediff=0,hypoinv=0,stdcc=0.0,stdct=0.0,bstart=0):
    """
    Run bootstrapping loop for ph2dt and hypoDD run (disk i/o turned off)
    :::
    Parameters:
    log (file object) --- Log file
    imakedata (int) --- Make synthetic model switch (set in run.inp)
    pinput (str) --- Ph2dt input file name (set in run.inp)
    hinput (str) --- HypoDD input file name (set in run.inp)
    nplot (int) --- No. of plotting iterations (set in run.inp)
    noisediff (int) --- Noise difference switch (set in run.inp)
    iplot (int) --- Plotting switch (set in run.inp)
    plotstrin (str) --- Plotting folder location (set in run.inp)
    hypoinv (int) --- Hypoinverse step switch (set in run.inp)
    nboot (int) --- Number of bootstrap iterations to run (set in run.inp)
    :::
    Return:
    None
    :::
    """

    """
    Set/make bootstrap output folder
    """
    bsfolder = str('bootstrap_output_%i' % nboot)
    if not os.path.isdir(bsfolder):
        os.mkdir(bsfolder)
    
    """
    Set plotting counter
    """
    plotcount=int(0)
    
    """
    Require limited files outputs and print statements for bootstrapping
    """
    #diskio=1

    """
    Start bootstap loop
    """
    print(bstart)
    if bstart!=0:
        for i in range(bstart):
            tmp = np.random.normal(loc=0.,scale=0.005,size=1000)
    #import pdb; pdb.set_trace()
    for iboot in range(bstart,nboot):
        log.write('Bootstrap iteration %i out of %i \n' % (iboot+1,nboot))
        print('Bootstrap iteration %i out of %i \n' % (iboot+1,nboot))
        plotcount+=1

        """
        Generate Synthetic Model or Read in Phases
        #############
        Pulls this step completely out of Ph2DT and into a separate function.
        Means we now pass the phase/event data into the Ph2DT function.
        """
        #log.write('Prepare for ph2dt function.  Call ph2dt_prep. \n')
        #[ph2dtinputs,nsta,s_lab,s_lat,s_lon,nev,lat,lon,depth,
        # cuspid,date,etime,mag,herr,verr,res,npha,nobs_ct,
        # p_pha,p_sta,p_time,p_wghtr,dist,az,ang,
        # rt_nl,rt_ratio,rt_top,rt_v,hypoDDinputs] = ph2dt_prep(log,imakedata,diskio,pinput,hinput)

        """
        Run ph2dt
        """
        #ph2dtlist = [ph2dtinputs[0],ph2dtinputs[1],ph2dtinputs[2],ph2dtinputs[3],
        #             ph2dtinputs[4],ph2dtinputs[5],nsta,s_lab,s_lat,s_lon,
        #             nev,lat,lon,depth,cuspid,npha,nobs_ct,p_pha,p_sta,p_time,p_wghtr]

        #[dt_sta,dt_dt,dt_qual,dt_c1,dt_c2,dt_idx,dt_ista,dt_ic1,dt_ic2,dt_offs,
        # ndt,nccp,nccs,nctp,ncts] = ph2dt(log,ph2dtlist,hypoDDinputs[9],imakedata,diskio)

        #ncc = nccp+nccs
        #nct = nctp+ncts

        """
        Add noise to dt_dt data
        """
        #if noiseswitch==1:
        #    [x,dt_dt,p_time] = generate_noise(noisediff,nsta,ncc,nct,stdcc,stdct,p_time,
        #                                      dt_dt,dt_ic1,dt_ic2,dt_ista,dt_idx)

        """
        Run hypoinverse catalog step if called
        """
        #print(lat)
        #if hypoinv==1:
        #    [dates,times,lat,lon,depth,mag,
        #     res,herr,verr,cuspid] = hypoinverse_catalog(log,diskio,p_time,dist,az,ang,
        #                                                 nsta,s_lab,s_lat,s_lon,
        #                                                 nev,cuspid,lat,lon,depth,
        #                                                 date,etime,mag,herr,verr,res,
        #                                                 rt_nl,rt_ratio,rt_top,rt_v,
        #                                                 hypoDDinputs[28],hypoDDinputs[29],
        #                                                 hypoDDinputs[30],hypoDDinputs[31])
        #print(lat)
        #sys.exit()

        """
        Run hypoDD
        """
        # Set input data into list
        #ph2data = [date,rtime,cuspid,lat,lon,depth,
        #           mag,herr,zerr,res,
        #           s_lab,s_lat,s_lon,
        #           dt_sta,dt_dt,dt_qual,dt_c1,dt_c2,dt_idx,
        #           dt_ista,dt_ic1,dt_ic2,dt_offs,
        #           nev,nsta,ndt,nccp,nccs,nctp,ncts]

        # Run hypoDD
        ibootstr = bsfolder + '/' + str(iboot+1) + '/'
        folder1 = '%stradout/' % ibootstr
        folder2 = '%stxtoutputs/' % ibootstr
        if not os.path.isdir(ibootstr):
            os.mkdir(ibootstr)
        if not os.path.isdir(folder1):
            os.mkdir(folder1)
        if not os.path.isdir(folder2):
            os.mkdir(folder2)

        #import pdb; pdb.set_trace()
        #try:
        hypoDD(log,hinput,reloctype=reloctype,iboot=ibootstr)
        #except:
        #    print('Error on bootstrap iteration: %i.' % iboot)
        #    print('Error exit from hypoDD.')
        #    log.write('Error on bootstrap iteration: %i. \n' % iboot)
        #    log.write('Error exit from hypoDD. \n')
        #    continue

        """        
        Track and save outputs
        ---
        Need to track locs and residuals for each iteration 
        Pass out of HypoDD and save to file
        """

        #savefiles(bsfolder,iboot,nev,calstart,calend,dtdtstart,dtdtend,
        #          loctru,locdel,locabs,x)
        #sys.exit()
        """
        Plotting
        ---
        If this is a plotting iteration and plotting
        is turned on
        """
        #if plotcount==nplot and iplot==1:
        #    print('Plotting to file for iteration %i' % int(iboot+1))
        #    print('Plots saved to folder %s' % (plotstrin))

        #    plotcount = plotting(plotstrin,calstart,calend,dtdt_start,dtdt_end,
        #                         locdel,locabs,loctru)

    return None





