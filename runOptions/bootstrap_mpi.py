#!/usr/bin/env python3
import os
import sys
from datetime import datetime
import shutil
import numpy as np
import glob
from mpi4py import MPI

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
cwd = os.get_cwd()
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

# Set number of data chunks to the number of bootstrap iterations
nch=nboot

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
log = open(os.path.join(outfol,'boot.log'),'w')
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
Run bootstrapping
:::
Bootstrapping requires fileout==1 for limited 
file/terminal input/outputs.
:::
"""
log.write('\n\nStarting bootstrapping... \n\n')
print('\n\nStarting bootstrapping... \n')

# Change variables to remove bootstrap call function
iplot = False
plotstrin = ''
imakedata = 0
noiseswitch = 0
noisediff = 0
hypoinv = 0
stdcc = 0.
stdct = 0.
bstart = 0

"""
Set/make bootstrap output folder
"""
bsfolder = str('bootstrap_output_%i' % nboot)
if not os.path.isdir(bsfolder):
    os.mkdir(bsfolder)

# Now to cross correlations
if rank==0:  # Master
    num_workers = nthreads-1
    closed_workers = 0
    ch = 0

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
                    index_send = np.zeros(2*len_ch+2,dtype=int) # Must reformat indexes to 1D array to be sent to worker
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
            print('Worker #',rank,' starting cc process #',chunk_num)

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
            print(bstart)
            if bstart!=0:
                for i in range(bstart):
                    tmp = np.random.normal(loc=0.,scale=0.005,size=1000)
            
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

            """
            Bootstrap function guts
            """
            flog.write('Bootstrap iteration %i out of %i \n' % (iboot+1,nboot))
            print('Bootstrap iteration %i out of %i \n' % (iboot+1,nboot))
            plotcount+=1

            hypoDD(flog,hinput,reloctype=reloctype,iboot=ibootstr)

            # Once Bootstrapping Done
            print('\nTime to run bootstrap %i: %i',(boot_num,(datetime.now()-tempstart).seconds))
            print('Sending process #',chunk_num,' back to root from worker #',rank)
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

