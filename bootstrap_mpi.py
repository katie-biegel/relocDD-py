#!/usr/bin/env python3
import os
import sys
from datetime import datetime
import shutil
import numpy as np
import glob
from mpi4py import MPI
import time

# Import needed hypoDD-py functions
from hypoDD.hypoDD import hypoDD
# Import run functions
from runOptions.runFunctions import run_input
# Import ph2dt functions
from ph2dt.ph2dt_files import ph2dt_input
# Import hypoDD functions
from hypoDD.hypoDD_files import hypoDD_input, readcts, readccs
from hypoDD.hypoDD import hypoDD
from utility.universal.raytrace import partials
from hypoDD.hypoDD_functions import dataprep, trialsrc,dtres


def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

bstart = 0

# File Locations
cwd = os.getcwd()
infol = os.path.join(cwd,'inpfiles')
outfol = os.path.join(cwd,'outputs')
datfol = os.path.join(cwd,'datfiles')
finput = os.path.join(infol,'rund.inp')

# Launch MPI
time_start = time.time()
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
status = MPI.Status()
nthreads = comm.Get_size()
tags = enum('READY', 'DONE', 'EXIT', 'START')

# Read in data on all channels
tstart = datetime.now()
# Get yesph2dt switch
[inputfol,datfol,outfol,reloctype,fileout,makedata,hypoinv,
 noiseswitch,noisediff,stdcc,stdct,nboot,nplot] = run_input(finput)

# Set number of data chunks to the number of bootstrap iterations
nch=nboot
if rank!=0:
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
    np.random.seed(bstart+rank+100)

"""
:::
Declare Folders for Outputs if Necessary
:::
"""
# Method specific outfolder
if reloctype==1:
    outfol = os.path.join(outfol,'EDD')
elif reloctype==2:
    outfol = os.path.join(outfol,'SDD')
elif reloctype==3:
    outfol = os.path.join(outfol,'DDD')
# Output file folders
tradout = os.path.join(outfol,'tradouts')
if rank==0:
    #if os.path.isdir(outfol):
    #    shutil.rmtree(outfol)
    os.makedirs(outfol,exist_ok='True')
    os.makedirs(tradout,exist_ok='True')
"""
:::
Open general log file
:::
"""
datet = str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))

"""
:::
Read in additional input files based on inputfol
:::
"""
# Set input file path names
pinput = os.path.join(inputfol,'ph2dt.inp')
hinput = os.path.join(inputfol,'hypoDD.inp')

# Read in ph2dt input file
pinputs = ph2dt_input('',pinput,dd_version=reloctype)
# Update path strings
pinputs[0] = os.path.join(datfol,pinputs[0])
pinputs[1] = os.path.join(datfol,pinputs[1])

# Read in hypoDD input files
hinputs = hypoDD_input('',hinput)
# Update data path strings
hinputs[0] = os.path.join(datfol,hinputs[0])
hinputs[1] = os.path.join(datfol,hinputs[1])
hinputs[2] = os.path.join(datfol,hinputs[2])
hinputs[3] = os.path.join(datfol,hinputs[3])
# Update hypoDD output strings 
hinputs[4] = os.path.join(tradout,hinputs[4])
hinputs[5] = os.path.join(tradout,hinputs[5])
hinputs[6] = os.path.join(tradout,hinputs[6])
hinputs[7] = os.path.join(tradout,hinputs[7])  
hinputs[8] = os.path.join(tradout,hinputs[8])
[fn_cc,fn_ct,fn_sta,fn_eve,fn_loc,fn_reloc,fn_res,fn_stares,fn_srcpar, # 8
 idata,iphase,minobs_cc,minobs_ct,amaxres_cross,amaxres_net,amaxdcc,amaxdct,maxdist, # 17
 awt_ccp,awt_ccs,awt_ctp,awt_cts,adamp,istart,maxiter,isolv,niter,aiter, # 27
 mod_nl,mod_ratio,mod_v,mod_top,iclust,ncusp,icusp] = hinputs

###
# Read in Best Fit Residual
###
datres = np.loadtxt(os.path.join(tradout,'hypoDD.res'),usecols=(6,),skiprows=1)
dattyp = np.loadtxt(os.path.join(tradout,'hypoDD.res'),usecols=(4,),skiprows=1,dtype='int')

datres_ccp = datres[dattyp==1]
datres_ccs = datres[dattyp==2]
datres_ctp = datres[dattyp==3]
datres_cts = datres[dattyp==4]

###
# Calculate traveltimes for best fit locations
###
bflocs = np.loadtxt(os.path.join(tradout,'hypoDD.reloc'),usecols=(1,2,3,))
bfcusp = np.loadtxt(os.path.join(tradout,'hypoDD.reloc'),usecols=(0,),dtype='int')
stlocs = np.loadtxt(os.path.join(datfol,'station.dat'),usecols=(1,2,))
stid = np.loadtxt(os.path.join(datfol,'station.dat'),usecols=(0,),dtype='str')
fn_srcpar = os.path.join(outfol,'bestfitlocs.srcpar')

[tmp_ttp,tmp_tts,
 tmp_xp,tmp_yp,tmp_zp] = partials(len(bflocs),bfcusp,bflocs[:,0],bflocs[:,1],bflocs[:,2],
                                  len(stlocs),stid,stlocs[:,0],stlocs[:,1],
                                  mod_nl,mod_ratio,mod_v,mod_top,fn_srcpar,return_all=False)

###
# Read in dt_dt
###
fileout = 0
maxsep_ct = amaxdct[0]
maxsep_cc = amaxdcc[0]

[ev_date,ev_time,ev_cusp,ev_lat,ev_lon,ev_dep,ev_mag,ev_herr,ev_zerr,ev_res,
 sta_lab,sta_lat,sta_lon,
 dt_sta1,dt_sta2,dt_dt,dt_qual,dt_c1,dt_c2,dt_idx,dt_ista1,dt_ista2,
 dt_ic1,dt_ic2,dt_offse,dt_offss,
 nev,nsta,ndt,ncc,nct,nccp,nccs,nctp,ncts] = dataprep(open(os.path.join(outfol,'bootmain.log'),'w'),
                                                      reloctype,fileout,fn_cc,fn_ct,fn_sta,
                                                      fn_eve,idata,iphase,ncusp,icusp,maxdist,
                                                      maxsep_ct,maxsep_cc)

# Make sure bfcusp matches evcusp
print('HypoDD.reloc: ',len(bfcusp))
print('HypoDD inp: ',len(ev_cusp))
for i in range(nev):
    if not bfcusp[i]==ev_cusp[i]:
        print('Break Boot Main: Resort events')

"""
Get initial residuals
"""
sdc0_lat = 0.
sdc0_lon = 0.
sdc0_dep = 0.
sdc0_lat = np.sum(bflocs[:,0])
sdc0_lon = np.sum(bflocs[:,1])
sdc0_dep = np.sum(bflocs[:,2])
sdc0_lat = sdc0_lat/nev
sdc0_lon = sdc0_lon/nev
sdc0_dep = sdc0_dep/nev

[nsrc,src_cusp,src_lat0,src_lon0,src_dep0,src_x0,src_y0,src_z0,src_t0,src_lat,
 src_lon,src_dep,src_x,src_y,src_z,src_t,src_xi,src_yi,src_zi,
 src_ti] = trialsrc(istart,sdc0_lat,sdc0_lon,sdc0_dep,nev,bfcusp,bflocs[:,0],bflocs[:,1],bflocs[:,2])

dt_boot,dt_res,tt = dtres(reloctype,ndt,nsrc,dt_dt,dt_idx,dt_ista1,dt_ista2,dt_ic1,dt_ic2,src_t,
                          tmp_ttp,tmp_tts)

# Change variables to remove bootstrap call function
iplot = False
plotstrin = ''
imakedata = 0
noiseswitch = 0
noisediff = 0
hypoinv = 0
stdcc = 0.
stdct = 0.
bend = 0
bsfolder = str('bootstrap_output_%i' % nboot)

# Now to cross correlations
if rank==0:  # Master
    num_workers = nthreads-1
    closed_workers = 0
    ch = 0
    if bstart!=0:
        ch = bstart
    if bend!=0:
        nboot = bend

    print('\n\nStarting bootstrapping... \n')
    """
    Set/make bootstrap output folder
    """
    if not os.path.isdir(bsfolder):
        os.mkdir(bsfolder)

    while closed_workers < num_workers: # Workers close when idle --- loop runs until all workers are closed
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        source = status.Get_source()
        tag = status.Get_tag()
        
        if tag==tags.READY:  # If worker signals it is ready
            print('Worker ',source,' ready.')
            if ch < nch:  # If there is still data to be sent
                if ch == nch-1:  # If last chunk --- different size
                    print('Sending task #',ch,' to worker ',source)
                    index_send = np.zeros(1,dtype=int) # Must reformat indexes to 1D array to be sent to worker
                    index_send[0] = ch  # First index is data chunk #
                    # Then send indexes
                    comm.send(index_send, dest=source, tag=tags.START) # Send to worker
                else:
                    print('Sending task #',ch,' to worker ',source)
                    index_send = np.zeros(1,dtype=int) # Must reformat indexes to 1D array to be sent to worker
                    index_send[0] = ch  # First index is data chunk #
                    # Then send indexes
                    comm.send(index_send, dest=source, tag=tags.START)  # Send to worker
            else:  # If there's no more work --- signal to close the worker
                print('Closing worker #',source)
                comm.send(None, dest=source, tag=tags.EXIT)
            ch += 1
        elif tag == tags.DONE: # If x-corr data is recieved from the worker --- i.e. process DONE
            data_chunk = int(data[0]) # Which data chunk received
            if data_chunk == nch-1: # If last chunk
                print('Process #',data_chunk,' from worker ',source,' completed.')
            else:
                print('Process #',data_chunk,' from worker ',source,' completed.')
        elif tag == tags.EXIT:  # If worker has been closed
           closed_workers += 1
    # Timed process ends here
    print('Time to run: ',time.time()-time_start)
elif rank>0:  # IF WORKER
    while True:
        comm.send(None, dest=0, tag=tags.READY)   # If idle send Ready
        boot_indexes = comm.recv(source=0,status=status)  # Recieve indexes for xcorr or exit tag
        tag = status.Get_tag()
        
        if tag == tags.START:  # If process to xcorr
            boot_num = boot_indexes[0]   # index 0 is the chunk #
            print('Worker #',rank,' starting cc process #',boot_num)

            """
            Run bootstrapping loop
            """
            tempstart = datetime.now()
            
            """
            Set plotting counter
            """
            plotcount=int(0)

            """
            Start bootstap loop
            """
            #print(bstart)
            #if bstart!=0:
            #    for i in range(bstart):
            #        tmp = np.random.normal(loc=0.,scale=0.005,size=1000)
            
            """
            Set boot iteration
            """
            iboot = boot_num

            # Run hypoDD
            ibootstr = os.path.join(bsfolder,str(iboot+1))
            folder1 = os.path.join(ibootstr,'tradouts')
            folder2 = os.path.join(ibootstr,'txtoutputs')
            if not os.path.isdir(ibootstr):
                os.mkdir(ibootstr)
            if not os.path.isdir(folder1):
                os.mkdir(folder1)
            if not os.path.isdir(folder2):
                os.mkdir(folder2)
            flog = open(os.path.join(ibootstr,'boot_%i.log' % boot_num),'w')  

            # Calculate noise from random choice
            bootres = np.zeros(len(dt_boot),dtype='float')
            for i in range(ndt):
                if dt_idx[i]==1:
                    bootres[i] = np.random.choice(datres_ccp/1000.,size=1,replace=True)    
                elif dt_idx[i]==2:
                    bootres[i] = np.random.choice(datres_ccs/1000.,size=1,replace=True)   
                elif dt_idx[i]==3: 
                    bootres[i] = np.random.choice(datres_ctp/1000.,size=1,replace=True) 
                elif dt_idx[i]==4:
                    bootres[i] = np.random.choice(datres_cts/1000.,size=1,replace=True) 
            bootres = bootres*(ndt/(ndt-4))

            np.savetxt(os.path.join(ibootstr,'noise_iter.txt'),bootres)

            """
            Bootstrap function guts
            """
            flog.write('Bootstrap iteration %i out of %i \n' % (iboot+1,nboot))
            print('Bootstrap iteration %i out of %i \n' % (iboot+1,nboot))
            plotcount+=1

            hypoDD(flog,hinputs,reloctype=reloctype,iboot=ibootstr,bootres=bootres,dtboot=dt_boot)

            # Once Bootstrapping Done
            print('\nTime to run bootstrap %i: %i' % (boot_num,(datetime.now()-tempstart).seconds))
            print('Sending process #',boot_num,' back to root from worker #',rank)
            send_back = np.zeros((1))
            send_back[0] = boot_num    # Index 0 is chunk #
            comm.send(send_back,dest=0,tag=tags.DONE) # Return to master

        elif tag == tags.EXIT:
            break # Break out of while loop to exit the process, exit tag sent from master

    comm.send(None, dest=0, tag=tags.EXIT)  # Exited worker

#comm.barrier()
sys.stdout.flush() # Flush system
#comm.Finalize()
sys.exit()

