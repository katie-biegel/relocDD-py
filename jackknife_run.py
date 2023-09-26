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
from hypoDD.hypoDD_files import hypoDD_input
from hypoDD.hypoDD import hypoDD

def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

# File Locations
cwd = os.getcwd()
infol = os.path.join(cwd,'inpfiles')
outfol = os.path.join(cwd,'outputs')
datfol = os.path.join(cwd,'datfiles')
finput = os.path.join(infol,'rune.inp')

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
    np.random.seed(rank)

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

#datres = np.loadtxt(os.path.join(tradout,'hypoDD.res'),usecols=(6,),skiprows=1)

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

# Change variables to remove bootstrap call function
iplot = False
plotstrin = ''
imakedata = 0
noiseswitch = 0
noisediff = 0
hypoinv = 0
stdcc = 0.
stdct = 0.
bstart = 100
bend = 200
bsfolder = str('jackknife_events')

# Read in Stations
#stas = np.loadtxt(os.path.join(tradout,'hypoDD.sta'),dtype='str',usecols=(0,))
# Set number of data chunks to the number of jackknife iterations
#nstas = len(stas)
#nch = nstas

# Read in Events
eves = np.loadtxt(os.path.join(tradout,'hypoDD.reloc'),dtype='str',usecols=(0,))
neves = len(eves)
nch = neves

# Now to cross correlations
if rank==0:  # Master
    num_workers = nthreads-1
    closed_workers = 0
    ch = bstart
    if bend!=0:
        nch = bend

    print('\n\nStarting jackknife (station testing)... \n')
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
            print('Worker #',rank,' starting jackknife process #',boot_num)

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
            ibootstr = os.path.join(cwd,bsfolder,str(eves[iboot]))
            folder1 = os.path.join(ibootstr,'tradouts')
            folder2 = os.path.join(ibootstr,'txtoutputs')
            if not os.path.isdir(ibootstr):
                os.mkdir(ibootstr)
            if not os.path.isdir(folder1):
                os.mkdir(folder1)
            if not os.path.isdir(folder2):
                os.mkdir(folder2)
            flog = open(os.path.join(ibootstr,'jackknife_%s.log' % eves[iboot]),'w')            

            """
            Bootstrap function guts
            """
            flog.write('Jackknife iteration %i out of %i \n' % (iboot+1,bend))
            print('Jackknife iteration %i out of %i \n' % (iboot+1,bend))
            plotcount+=1

            hypoDD(flog,hinputs,reloctype=reloctype,iboot=ibootstr,jke=iboot)

            # Once Bootstrapping Done
            print('\nTime to run jackknife %i: %i' % (boot_num,(datetime.now()-tempstart).seconds))
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

