import numpy as np
import datetime
import sys
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import interpolate

from runOptions.runFunctions import run_input
from hypoDD.hypoDD_functions import cluster,dataprep,trialsrc,dtres,weighting,resstat
from hypoDD.hypoDD_functions import skip,sigcoherency,statres,eventstats
from hypoDD.inversion import lsqr
from hypoDD.hypoDD_files import hypoDD_input, terminaloutputs
from utility.universal.geodetics import setorg,sdc2
from utility.universal.raytrace import partials


if len(sys.argv)==5:
    rinput = str(sys.argv[1])
    minlam = float(sys.argv[2])
    maxlam = float(sys.argv[3])
    nlam = int(sys.argv[4])
else:
    rinput = str(input('Input File: '))
    minlam = float(input('Min Damp to Test: '))
    maxlam = float(input('Max Damp to Test: '))
    nlam = int(input('No. Damp to Test: '))


"""
Read in default input file
"""
[inputfol,datfol,outfol,reloctype,fileout,makedata,hypoinv,
 noiseswitch,noisediff,stdcc,stdct,nboot,nplot] = run_input(rinput)
hinput = os.path.join(inputfol,'hypoDD.inp')
pinput = os.path.join(inputfol,'ph2dt.inp')

if reloctype==1:
    outfol = os.path.join(outfol,'EDD')
elif reloctype==2:
    outfol = os.path.join(outfol,'SDD')
elif reloctype==3:
    outfol = os.path.join(outfol,'DDD')
else:
    outfol = os.path.join(outfol,'EDD')

log = open(os.path.join(outfol,'damping.log'),'w')
log.write('Damping tests for: \n\n')
log.write('Min. damping value tested: %f \n' % (10.**minlam))
log.write('Max. damping value tested: %f \n' % (10.**maxlam))
log.write('No. damping value tested: %i \n' % (nlam))

[fn_cc,fn_ct,fn_sta,fn_eve,
 fn_loc,fn_reloc,fn_res,fn_stares,fn_srcpar,
 idata,iphase,minobs_cc,minobs_ct,
 amaxres_cross,amaxres_net,amaxdcc,amaxdct,maxdist,
 awt_ccp,awt_ccs,awt_ctp,awt_cts,adamp,
 istart,maxiter,isolv,niter,aiter,
 mod_nl,mod_ratio,mod_v,mod_top,
 iclust,ncusp,icusp] = hypoDD_input(log,hinput,fileout)

"""
Read in Data
"""
fn_cc = os.path.join(datfol,fn_cc)
fn_ct = os.path.join(datfol,fn_ct)
fn_sta = os.path.join(datfol,fn_sta)
fn_eve = os.path.join(datfol,fn_eve)
txtfol = os.path.join(outfol,'txtoutputs')

[ev_date,ev_time,ev_cusp,ev_lat,ev_lon,ev_dep,ev_mag,
 ev_herr,ev_zerr,ev_res,sta_lab,sta_lat,sta_lon,
 dt_sta1,dt_sta2,dt_dt,dt_qual,dt_c1,dt_c2,dt_idx,
 dt_ista1,dt_ista2,dt_ic1,dt_ic2,dt_offse,dt_offss,
 nev,nsta,ndt,ncc,nct,nccp,nccs,nctp,
 ncts]= dataprep(log,reloctype,fileout,fn_cc,fn_ct,
                 fn_sta,fn_eve,idata,iphase,ncusp,
                 icusp,maxdist,amaxdct[0],amaxdcc[0])

"""
Declare Open Variables
"""
tmpr1 = float(0.)
tmpr2 = float(0.)
minwght = float(0.00001)
rms_ccold = float(0.)
rms_ctold = float(0.)
rms_cc0old = float(0.)
rms_ct0old = float(0.)
ineg = int(0)
rms_cc = float(0.)
rms_ct = float(0.)
rms_cc0 = float(0.)
rms_ct0 = float(0.)
xav = float(0.)
yav = float(0.)
zav = float(0.)
tav = float(0.)
sav = float(0.)
alat = float(0.)
alon = float(0.)
adep = float(0.)

""""
Initial Clustering
"""
if (idata==1 and minobs_cc==0) or (idata==2 and minobs_ct==0) or (idata==3 and minobs_ct+minobs_cc==0):
    # No clustering (one cluster includes all events)
    nclust = 1
    clust = np.zeros((1,nev+1))
    clust[0,0] = nev
    for i in range(0,nev):
        clust[0,i+1] = ev_cusp[i]
    log.write('(hypoDD) No clustering performed. \n')
    if fileout==0:
        print('(hypoDD) No clustering performed.')
else: # Clustering
    if reloctype==2:
        nclust = 1
        clust = np.zeros((1,nev+1))
        clust[0,0] = nev
        for i in range(0,nev):
            clust[0,i+1] = ev_cusp[i]

        nclust = 1
        clust = np.zeros((nclust,nev+1),dtype='int')
        for i in range(1):
            evtoclust = open(os.path.join(outfol,('txtoutputs/clust_%i.dat' % i)),'r')
            clustlines = evtoclust.readlines()
            evtoclust.close()

            count = 1
            for line in clustlines:
                line = line.strip()
                line = line.split(' ')
                line = list(filter(None,line))

                clust[i,count:count+len(line)] = line[0:len(line)]
                count += len(line)
            clust[i,0] = count-1
    else:
        [clust,noclust,nclust] = cluster(log,txtfol,nev,ndt,idata,minobs_cc,
                                         minobs_ct,dt_c1,dt_c2,ev_cusp)

"""
Initalise loop over iterations
"""
jiter=int(0)
if iclust!=0:
    if iclust<0 or iclust>nclust:
        raise Exception('(hypoDD) Error: invalid cluster number %5i. Must be between 1 and nclust (%5i)' % (iclust,nclust))
    ibeg=iclust-1
    iend=iclust
else:
    ibeg=int(0)
    iend=nclust

# Initialise airquake array (cuspid's for airquakes later stored here)
amcusp = np.zeros(1000,dtype='int')

"""
Initialize variables
"""
exav=int(0)
eyav=int(0)
ezav=int(0)
etav=int(0)
esav=int(0)
dxav=int(0)
dyav=int(0)
dzav=int(0)
dtav=int(0)
dsav=int(0)
# Initialize arrays for storing station residuals statistics
sta_np = np.zeros(nsta,dtype='int')
sta_ns = np.zeros(nsta,dtype='int')
sta_nnp = np.zeros(nsta,dtype='int')
sta_nns = np.zeros(nsta,dtype='int')
sta_rmsc = np.zeros(nsta,dtype='float')
sta_rmsn = np.zeros(nsta,dtype='float')

"""
Big cluster loop
"""
adamp_tmp = np.logspace(minlam,maxlam,nlam)
dampnorms = np.zeros((len(adamp_tmp),3),dtype='float')

for dampindx,damptmp in enumerate(adamp_tmp):
    print('\n\nDamp Loop %i --- Testing Damping Value: %f' % (dampindx,damptmp))
    avgresnorm = []
    avgsolnorm = []
    maxiter=3
    icl=0
    time_cluststart = datetime.datetime.now()
    datet = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    log.write('\n\n(hypoDD) RELOCATION OF CLUSTER: %2i     %s \n\n\n' % (icl+1,datet))
    print('\n(hypoDD) RELOCATION OF CLUSTER: %2i     %s' % (icl+1,datet))

    """
    Get data for each cluster if clustering invoked
    -----
    Reads in only correct data for each cluster
    """
    if (nclust!=1) and (minobs_cc>0 or minobs_ct>0):
        ncusp = clust[icl,0]
        icusp = clust[icl,1:ncusp+1] #np.zeros(ncusp,dtype='int')
        
        [ev_date,ev_time,ev_cusp,ev_lat,ev_lon,ev_dep,ev_mag,
         ev_herr,ev_zerr,ev_res,sta_lab,sta_lat,sta_lon,
         dt_sta1,dt_sta2,dt_dt,dt_qual,dt_c1,dt_c2,dt_idx,
         dt_ista1,dt_ista2,dt_ic1,dt_ic2,dt_offse,dt_offss,
         nev,nsta,ndt,ncc,nct,nccp,nccs,nctp,
         ncts]= dataprep(log,reloctype,fileout,fn_cc,fn_ct,
                         fn_sta,fn_eve,idata,iphase,ncusp,icusp,
                         maxdist,amaxdct[0],amaxdcc[0])
    # Recount all data (to update for cluster subset)
    nccold = ncc
    nctold = nct
    nevold = nev

    """
    Get cluster centroid
    """
    sdc0_lat = 0.
    sdc0_lon = 0.
    sdc0_dep = 0.
    sdc0_lat = np.sum(ev_lat)
    sdc0_lon = np.sum(ev_lon)
    sdc0_dep = np.sum(ev_dep)
    sdc0_lat = sdc0_lat/nev
    sdc0_lon = sdc0_lon/nev
    sdc0_dep = sdc0_dep/nev
    log.write('Cluster centroid at: %10.6f  %11.6f  %9.6f \n' % (sdc0_lat,sdc0_lon,sdc0_dep))

    """"
    Set up cartesian coordinates from cluster centroid epicenter lat and lon
    """
    setorg(sdc0_lat,sdc0_lon)

    """
    Convert all events to cartesian
    """
    ev_x = np.zeros(nev)
    ev_y = np.zeros(nev)
    ev_z = np.zeros(nev)
    for i in range(0,nev):
        lat = ev_lat[i]
        lon = ev_lon[i]
        [x,y] = sdc2(lat,lon,-1)
        ev_x[i] = x*1000
        ev_y[i] = y*1000
        ev_z[i] = (ev_dep[i]-sdc0_dep)*1000
    log.write('# Events: %5i \n' % nev)

    """
    Get initial trial sources
    """
    [nsrc,src_cusp,src_lat0,src_lon0,src_dep0,
     src_x0,src_y0,src_z0,src_t0,
     src_lat,src_lon,src_dep,
     src_x,src_y,src_z,src_t,
     src_xi,src_yi,src_zi,src_ti] = trialsrc(istart,sdc0_lat,sdc0_lon,sdc0_dep,
                                             nev,ev_cusp,ev_lat,ev_lon,ev_dep)
    src_s = np.zeros(nev,dtype='float')
    log.write('(hypoDD) # Initial trial sources: %6i \n' % nsrc)

    """
    Loop over iterations starts here:
    """
    for i in range(1):
        #import pdb; pdb.set_trace()
        #for i in range(0,niter):
        aiter[0] = aiter[0] - jiter     # Aiter from input file
                                     # Number of iterations per set
        maxiter = maxiter - jiter 
        kiter = 0   # counter for iter with data skipping
        jiter = 0   # counter for iter with no updating (air quakes)
        mbad = 0    # counter for air quakes
        normvar=np.zeros((int(10*maxiter),3))
        iteri = 0   # start on iteration 1
        time_preit = datetime.datetime.now()

        """
        Get weighting parameters for this iteration:
        """
        maxres_cross = amaxres_cross[0] # Max. residual threshold (corr)
        maxres_net = amaxres_net[0]     # Max. residual threshold (cat)
        maxdcc = amaxdcc[0]             # Max. interevent sep for linked pairs (corr)
        maxdct = amaxdct[0]             # Max. interevent sep for linked pairs (cat)
        wt_ccp = awt_ccp[0]             # Weight multiplier for cc_p (corr)
        wt_ccs = awt_ccs[0]             # Weight multiplier for cc_s (corr)
        wt_ctp = awt_ctp[0]             # Weight multiplier for ct_p (cat)
        wt_cts = awt_cts[0]             # Weight multiplier for ct_s (cat)
        damp = damptmp

        # Write inversion parameters to log
        log.write('Weighting parameters for this iteration: \n')
        log.write('wt_ccp= %7.4f  wt_ccs= %7.4f  maxr_cc= %7.4f  maxd_cc= %7.4f \n' % 
                  (wt_ccp,wt_ccs,maxres_cross,maxdcc))
        log.write('wt_ctp= %7.4f  wt_cts= %7.4f  maxr_ct= %7.4f  maxd_ct= %7.4f  damp= %7.4f \n\n' % 
                  (wt_ctp,wt_cts,maxres_net,maxdct,damp))

        """
        Calculate travel times and slowness vectors
        """
        [tmp_ttp,tmp_tts,
         tmp_xp,tmp_yp,tmp_zp] = partials(nsrc,src_cusp,src_lat,src_lon,src_dep,
                                          nsta,sta_lab,sta_lat,sta_lon,
                                          mod_nl,mod_ratio,mod_v,mod_top,fn_srcpar)
        """
        Calculate double difference vector and theor dt
        """
        [dt_cal,dt_res,tt] = dtres(reloctype,ndt,nsrc,dt_dt,dt_idx,
                                   dt_ista1,dt_ista2,dt_ic1,dt_ic2,
                                   src_t,tmp_ttp,tmp_tts)
        """
        Get a priori weights and reweight residuals
        """
        [ineg,dt_wt,kiter] = weighting(log,reloctype,ndt,mbad,amcusp,idata,kiter,ineg,
                                       maxres_cross,maxres_net,maxdcc,maxdct,minwght,
                                       wt_ccp,wt_ccs,wt_ctp,wt_cts,
                                       dt_c1,dt_c2,dt_idx,dt_qual,dt_res,
                                       dt_offse,dt_offss)
        """
        Skip outliers and/or air quakes
        """
        if ineg>0 or reloctype==2:
            [ndt,nev,nsrc,nsta,
             ev_cusp,ev_date,ev_time,ev_mag, 
             ev_lat,ev_lon,ev_dep,
             ev_x,ev_y,ev_z,
             ev_herr,ev_zerr,ev_res,
             src_cusp,src_lat,src_lon,src_dep,
             src_lat0,src_lon0,
             src_x,src_y,src_z,src_t,
             src_x0,src_y0,src_z0,src_t0,
             src_xi,src_yi,src_zi,src_ti,
             sta_lab,sta_lat,sta_lon,
             sta_rmsc,sta_rmsn,
             sta_np,sta_ns,sta_nnp,sta_nns,
             dt_sta1,dt_sta2,dt_c1,dt_c2,
             dt_idx,dt_dt,dt_qual,dt_cal,
             dt_ista1,dt_ista2,dt_ic1,dt_ic2,
             dt_res,dt_wt,dt_offse,dt_offss,
             tmp_ttp,tmp_tts,tmp_xp,tmp_yp,tmp_zp,
             nct,ncc] = skip(log,reloctype,kiter,minwght,ndt,nev,nsrc,nsta,
                             ev_cusp,ev_date,ev_time,ev_mag,ev_lat,ev_lon,ev_dep,
                             ev_x,ev_y,ev_z,ev_herr,ev_zerr,ev_res,
                             src_cusp,src_lat,src_lon,src_dep,src_lat0,src_lon0,
                             src_x,src_y,src_z,src_t,src_x0,src_y0,src_z0,src_t0,
                             src_xi,src_yi,src_zi,src_ti,
                             sta_lab,sta_lat,sta_lon,sta_rmsc,sta_rmsn,
                             sta_np,sta_ns,sta_nnp,sta_nns,
                             dt_sta1,dt_sta2,dt_c1,dt_c2,dt_idx,dt_dt,dt_qual,dt_cal,
                             dt_ista1,dt_ista2,dt_ic1,dt_ic2,dt_res,dt_wt,dt_offse,dt_offss,
                             tmp_ttp,tmp_tts,tmp_xp,tmp_yp,tmp_zp,nct,ncc,amcusp,mbad)
        """
        Get initial residual statistics (avrg,rms,var...)
        """
        [rms_cc,rms_ct,rms_cc0,rms_ct0,rms_ccold,rms_ctold,
         rms_cc0old,rms_ct0old,resvar1] = resstat(log,reloctype,idata,ndt,nev,dt_res,dt_wt,dt_idx,
                                                  rms_cc,rms_ct,rms_cc0,rms_ct0,rms_ccold,
                                                  rms_ctold,rms_cc0old,rms_ct0old)
        """
        Call inversion
        """
        [src_cusp,src_dx,src_dy,src_dz,src_dt,src_ds, 
         src_ex,src_ey,src_ez,src_et,src_es,
         exav,eyav,ezav,etav,esav,
         dxav,dyav,dzav,dtav,dsav,
         rms_cc,rms_ct,rms_cc0,rms_ct0,
         rms_ccold,rms_ctold,rms_cc0old,rms_ct0old,
         acond,normvar,sam] = lsqr(log,reloctype,iteri,ndt,nev,nsrc,damp,mod_ratio, 
                                     idata,ev_cusp,src_cusp,dt_res,dt_wt,
                                     dt_ista1,dt_ista2,dt_ic1,dt_ic2,
                                     exav,eyav,ezav,etav,esav,
                                     dxav,dyav,dzav,dtav,dsav,
                                     rms_cc,rms_ct,rms_cc0,rms_ct0,
                                     rms_ccold,rms_ctold,rms_cc0old,rms_ct0old,
                                     tmp_xp,tmp_yp,tmp_zp,dt_idx,normvar)
        print('R1Norm of ',sam[0],' for damp value of ',damptmp)
        print('Xnorm of ',sam[1])
        print('R2Norm of ',sam[2])
        dampnorms[dampindx,0] += sam[0]
        dampnorms[dampindx,1] += sam[1]
        dampnorms[dampindx,2] += sam[2]

        #import pdb; pdb.set_trace()
        """
        Check for air quakes
        -----
        Air quakes are relocated events with a negative depth
        Negative depths indicate a depth above 0.0
        Air quake IDs are saved to amcusp and then removed before inversion
        in the next iteration
        """
        mbad = 0
        k = 0
        for i in range(nsrc):
            if (src_dep[i]+src_dz[i]/1000.)<0.: # If depth negative
                log.write('(hypoDD) >>> Warning: negative depth - %12i \n' % ev_cusp[i])
                amcusp[k]=ev_cusp[i]  # Save event ID
                k += 1
                if k>1000:  # Trigger exception if there are too many airquakes.
                            # Indicates poor data quality or poorly constrained inversion              
                    raise Exception('(hypoDD) >>> More than 1000 air quakes.')
        mbad = k # Number of neg depth events

        """
        Update iteration numbers:
        ----
        If there are airquakes the iteration DOES NOT COUNT
        If there are airquakes the event locations will not be updated (iskip)
        """
        iskip=0
        if mbad>0:
            for i in range(niter):
                aiter[i] = aiter[i] + 1
            jiter+=1 # iteration with no update
            maxiter+=1
            normvar = np.append(normvar,np.zeros((1,3)),axis=0)

            log.write('(hypoDD) Number of air quakes (AQ) = %i \n' % mbad)
            if (nsrc-mbad)<=1:
                log.write('(hypoDD) Warning: number of non-airquakes < 2.  Skipping cluster. \n')
                print('(hypoDD) Warning: number of non-airquakes < 2.  Skipping cluster. \n')
                continue
            iskip=1

        """
        Update source parameters:
        """
        if iskip==0:
            xav = 0 # mean centroid shift
            yav = 0
            zav = 0
            tav = 0
            sav = 0
            alon = 0
            alat = 0
            adep = 0

            if nsrc==1:
                nsrc=nev

            # Save iteration updates to file
            # # Open all files
            # if fileout==0:
            #     fds = "txtoutputs/delta_source_"+str(icl+1)+'_'+str(iteri+1-jiter)+".txt"
            #     delta_source = open(fds,'w')
            #     fas = "txtoutputs/abs_source_"+str(icl+1)+'_'+str(iteri+1-jiter)+".txt"
            #     abs_source = open(fas,'w')
            #     fini = "txtoutputs/abs_source_"+str(icl+1)+".txt"
            #     abs_source0 = open(fini,'w')
            #     fts = "txtoutputs/tru_source_"+str(icl+1)+'_'+str(iteri+1-jiter)+".txt"
            #     tru_source = open(fts,'w')

            # Loop over events
            for i in range(0,nsrc):
                src_cusp[i] = ev_cusp[i]
                # Update absolute source parameters
                src_x[i] = src_x[i] + src_dx[i]
                src_y[i] = src_y[i] + src_dy[i]
                src_z[i] = src_z[i] + src_dz[i]
                if reloctype==1:
                    src_t[i] = src_t[i] + src_dt[i]
                if reloctype==2:
                    src_s[i] = src_s[i] + src_ds[i]
                # Save updates and file locations to file
                # if fileout==0:
                #     delta_source.write('%10.6f %10.6f %10.6f %10.6f %10.6f \n' % (src_dx[i],src_dy[i],src_dz[i],src_dt[i],src_ds[i]))
                #     abs_source.write('%10.6f %10.6f %10.6f %10.6f %10.6f \n' % (src_x[i],src_y[i],src_z[i],src_t[i],src_s[i]))
                #     abs_source0.write('%10.6f %10.6f %10.6f %10.6f %10.6f \n' % (src_x0[i],src_y0[i],src_z0[i],src_t0[i],src_t0[i]))
                #     tru_source.write('%10.6f %10.6f %10.6f %10.6f \n' % (src_xi[i],src_yi[i],src_zi[i],src_ti[i]))
                # Save delta and abs locations into arrays
                # locdel[:,i,iteri-jiter] = np.array([src_dx[i],src_dy[i],src_dz[i],src_dt[i],src_ds[i]])
                # locabs[:,i,iteri-jiter+1] = np.array([src_x[i],src_y[i],src_z[i],src_t[i],src_s[i]])
                # Also save true locations
                # locabs[:,i,0] = np.array([src_x0[i],src_y0[i],src_z0[i],src_t0[i],src_t0[i]])
                # loctru[:,i] = np.array([src_xi[i],src_yi[i],src_zi[i],src_ti[i]])
                # Update absolute source locations both in
                # cartesian coordinates and in lat/lon
                src_dep[i] = src_dep[i] + src_dz[i]/1000
                [lat,lon] = sdc2(src_x[i]/1000,src_y[i]/1000,1)
                src_lon[i] = lon
                src_lat[i] = lat
                alon = lon+alon
                alat = lat+alat
                adep = adep+src_dep[i]
                # Calculate the centroid shift
                ################## KB Note: Pull this out of loop maybe 
                xav = xav + (src_x[i] - src_x0[i])
                yav = yav + (src_y[i] - src_y0[i])
                zav = zav + (src_z[i] - src_z0[i])
                tav = tav + (src_t[i] - src_t0[i])
                sav = sav + (src_s[i])
            #import pdb; pdb.set_trace()
            xav = xav/nsrc
            yav = yav/nsrc
            zav = zav/nsrc
            tav = tav/nsrc
            sav = sav/nsrc
            alon = alon/nsrc
            alat = alat/nsrc
            adep = adep/nsrc

            # Document centroid shift
            log.write('\n(hypoDD) Cluster centroid at: %10.6f  %11.6f  %9.6f \n' % 
                      (alat,alon,adep))
            log.write('(hypoDD) Mean centroid (origin) shift in x, y, z, t [m,ms]: %7.1f, %7.1f, %7.1f, %7.1f \n' % 
                      (xav,yav,zav,tav))
            # Close iteration files
            # if fileout==0:
            #     delta_source.close()
            #     abs_source.close()
            #     tru_source.close()

            """
            Get interpair distance for each observation and average signal coherency
            """
            [ncc,nct,cohav,picav,
             dt_offse] = sigcoherency(log,reloctype,nct,ncc,ndt,idata,src_x,src_y,src_z,
                                      dt_ic1,dt_ic2,dt_offse,dt_qual,dt_idx)

            """
            Get number of observations and mean residual at each station
            ----
            Resstat calculates residuals statistics per event
            This loop calculates residual statistics per station
            """
            [sta_np,sta_ns,sta_nnp,sta_nns,
             sta_rmsc,sta_rmsn,tmpr1,tmpr2] = statres(log,nsta,ndt,idata,reloctype,
                                                      sta_lab,dt_ista1,dt_ista2,
                                                      dt_idx,dt_res)


    dampnorms[dampindx,0] = dampnorms[dampindx,0]#/3.
    dampnorms[dampindx,0] = dampnorms[dampindx,0]#/3.
    dampnorms[dampindx,0] = dampnorms[dampindx,0]#/3.


# Save Outputs in Case Plotting Breaks
np.savetxt(os.path.join(outfol,'DampTested.dat'),adamp_tmp)
np.savetxt(os.path.join(outfol,'R1NormsDamping.dat'),dampnorms[:,0])
np.savetxt(os.path.join(outfol,'XNormsDamping.dat'),dampnorms[:,1])
np.savetxt(os.path.join(outfol,'R2NormsDamping.dat'),dampnorms[:,2])

# Spline Fit The L-Curve
cutoff = len(adamp_tmp)
x1 = dampnorms[:cutoff,0]
y1 = dampnorms[:cutoff,1]

logadamp = np.log10(adamp_tmp[1:cutoff])
logx1 = np.log10(dampnorms[1:cutoff,0])
logy1 = np.log10(dampnorms[1:cutoff,1])
dampnew = np.linspace(0,adamp_tmp[cutoff-1],1000)
#dampnew = dampnew[1:]
logdamp = np.log10(dampnew[1:])

#fresnorm = interpolate.UnivariateSpline(logadamp,logx1,k=5)
#fsolnorm = interpolate.UnivariateSpline(logadamp,logy1,k=5)
fresnorm = interpolate.UnivariateSpline(adamp_tmp,x1,k=5)
fsolnorm = interpolate.UnivariateSpline(adamp_tmp,y1,k=5)
#x = fresnorm(logdamp[34:])
#y = fsolnorm(logdamp[34:])
x = fresnorm(dampnew)
y = fsolnorm(dampnew)

# Plot L-Curve
fig,ax = plt.subplots(1,2,figsize=(11,5))
ax[0].plot(logx1,logy1,'ko-')
ax[0].plot(np.log10(x),np.log10(y),'r--')
ax[0].set_xlabel('$log_{10}(||Am-d||_2)$')
ax[0].set_ylabel('$log_{10}(||m||_2)$')

# Curvature of the Log-Scaled L-Curve
#
# K = 2 * (xhat'yhat" - xhat"yhat') / ((xhat')**2 + (yhat')**2)**(3/2.)
#
# xhat'  = x'/x
# yhat'  = y'/y
# xhat" = (x"x -(x')**2)/x**2
# yhat" = (y"y -(y')**2)/y**2
#
# Take derivatives
#lfresnorm = interpolate.UnivariateSpline(logadamp,x1,k=5)
#lfsolnorm = interpolate.UnivariateSpline(logadamp,y1,k=5)
fderx1 = fresnorm.derivative(n=1)
fderx2 = fresnorm.derivative(n=2)
fdery1 = fsolnorm.derivative(n=1)
fdery2 = fsolnorm.derivative(n=2)
#x = fresnorm(logdamp)
#y = fsolnorm(logdamp)
#dx1 = fderx1(logdamp)
#dx2 = fderx2(logdamp)
#dy1 = fdery1(logdamp)
#dy2 = fdery2(logdamp) 
dx1 = fderx1(dampnew)/x
dx2 = (fderx2(dampnew)*x - fderx1(dampnew)**2)/x**2
dy1 = fdery1(dampnew)/y
dy2 = (fdery2(dampnew)*y - fdery1(dampnew)**2)/y**2

top = (dx1*dy2 - dx2*dy1)
bottom = ((dx1)**2 + (dy1)**2)**(3/2.)
K = 2*top/bottom

# Lambda Selection
print('Max curvature: ',np.max(K))
print('Damp at max curvature: ',(dampnew[np.argmax(K)]))
#log.write('\n\nMax curvature: %f \n' % np.max(K))
#log.write('Damp at max curvature: %f \n' % dampnew[np.argmax(K)])

ax[1].plot(np.log10(dampnew[1:]),K[1:],'k-')
ax[1].axvline(np.log10(dampnew[np.argmax(K)]),color='r',linestyle=':',alpha=0.7,linewidth='2')
ax[1].set_xlabel('$log_{10}(\lambda)$')
ax[1].set_ylabel('$Curvature (\kappa)$')

plt.tight_layout()
plt.savefig('LCurve_log.png')
plt.close('all')

# Load from file
adamp_tmp = np.loadtxt('DampTested.dat')
dampnorms = np.zeros((len(adamp_tmp),3))
dampnorms[:,0] = np.loadtxt('R1NormsDamping.dat')
dampnorms[:,1] = np.loadtxt('XNormsDamping.dat')
dampnorms[:,2] = np.loadtxt('R2NormsDamping.dat')




