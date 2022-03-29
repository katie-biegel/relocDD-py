import numpy as np

"""
---
This script contains minor subfunctions called from the main run function.
---
This script contains the functions:
	run_input:		Reads the run.inp file
    readinfiles:    Read in relocation ouput files if needed
---
"""

def run_input(inputfile='inpfiles/run.inp'):
    """
    Read in run file input file
    :::
    Parameters:
    inputfile (str) --- File location of run input file
                        Default to 'inpfiles/run.inp'
    :::
    Returns:
    retlist (list) --- List of Run input variables
    """

    """
    Open file and readlines
    """
    inputs = open(inputfile,'r')
    cinp = inputs.readlines()
    iline = 0

    """
    Declare run.inp variables
    :::
    For each line in run.inp try to read in the
    value otherwise set the value to default values.
    :::
    """
    print('\nReading in run.inp file...')
    for var in cinp:
        """
        Skip Comment Lines
        """
        if var[0]=='*':
            continue 

        var = var.strip()

        """ 
        Folder Locations For Inputs/Outputs
        """  
        # Input File Folder - holds ph2dt.inp and hypoDD.inp
        if iline==0:
            try:
                inputfol = str(var)
                if len(inputfol)==0:
                    inputfol = str('inpfiles')  
                    print('Default input file locations to "inpfiles/".\n')
            except:
                inputfol = str('inpfiles')  
                print('Default input file locations to "inpfiles/".\n')
        # Data File Folder - holds station.dat, phase.dat, and any other data files.
        if iline==1:
            try: 
                datfol = str(var)
                if len(datfol)==0:
                    datfol = str('datfiles') 
                    print('Default data file locations to "datfiles/".\n')
            except:
                datfol = str('datfiles') 
                print('Default data file locations to "datfiles/".\n')
        # Output File Folder - folder relocation outputs are written to.
        if iline==2:
            try: 
                outfol = str(var)
                if len(outfol)==0:
                    outfol = str('outputs') 
                    print('Default output file locations to "outputs/".\n')
            except:
                outfol = str('outputs') 
                print('Default output file locations to "outputs/".\n')
        
        """
        Switch Variables
        """
        # Type of reloc switch (type of pairing used - defaults to eventpair (hypoDD))
        if iline==3:
            try:
                reloctype= int(var)  
                if reloctype!=1 and reloctype!=2 and reloctype!=3:    # Only valid values
                    reloctype=1
                    print('Default reloctype to event-pair.\n')   
            except:  
                reloctype=1
                print('Default reloctype to event-pair.\n')  
        # File outputs switch (limits disk,file,print i/o if on)   
        if iline==4:
            try:
                fileout = int(var)          
                if fileout!=0 and fileout!=1:
                   fileout=0
                   print('Default fileout to traditional (on).\n')  
            except:
                fileout=0
                print('Default fileout to traditional (on).\n')  
        # Synthetic Modelling Switch (turns on/off synthetic modeling)
        if iline==5:                
            try:
                makedata = int(var)
                if makedata!=0 and makedata!=1:
                    makedata=0
                    print('Default makedata to real data (off).\n')  
            except:
                makedata=0
                print('Default makedata to real data (off).\n')  
        # Hypoinverse Switch (only used in synthetic modeling cases)
        if iline==6:
            try:
                hypoinv = int(var)
                if hypoinv!=0 and hypoinv!=1:
                    hypoinv=0
                    print('Default hypoinv to no hypoinverse step (off).\n')
            except:         
                hypoinv=0
                print('Default hypoinv to no hypoinverse step (off).\n')

        """
        Noise Variables
        """

        # Noise Switch (only used in synthetic modeling cases)
        if iline==7:
            try:
                noiseswitch = int(var)    
                if noiseswitch!=0 and noiseswitch!=1:
                    noiseswitch=0
                    print('Default noiseswitch to no added noise (off).\n')   
            except:
                noiseswitch = 0
                print('Default noiseswitch to no added noise (off).\n') 
        # Noise difference switch 
        # ---
        # On if noise is added to traveltime measurements (i.e. data noise)
        # Off if noise is added to diff. time. measurements (i.e. theory noise)
        if iline==8:
            try:
                noisediff = int(var) 
                if noisediff!=0 and noisediff!=1:
                    noisediff=0
                    print('Default noisediff to single noise value (off).\n')   
            except:
                noisediff = 0
                print('Default noisediff to single noise value (off).\n')  
        # Std. of gaussian noise for cc measurements 
        if iline==9:
            try:
                stdcc = float(var)
            except:
                stdcc = 0.001
                print('Default stdcc to 0.001 (1 ms).\n')   
        # Std. of gaussian noise for ct measurements 
        if iline==10:
            try:
                stdct = float(var)
            except:
                stdct = 0.01
                print('Default stdct to 0.01 (10 ms).\n')  

        """
        Bootstrap Variables
        """
        
        # Number of bootstrap iterations
        if iline==11:
            try:
                nboot = int(var)
            except:
                nboot = 10
                print('Default nboot to 10 iterations.\n')  
        # Plotting step size for bootstrapping
        # ---
        # For periodic outputs from bootstrapping to check progression.
        if iline==12:
            try:
                nplot = int(var)
            except:
                nplot = 10
                print('Default nplot to 10 iterations.\n')  

        """
        Count up Line
        """
        iline += 1

    inputs.close()
    print('Run.inp read into memory.')

    """
    Return input parameters
    """
    retlist = [inputfol,datfol,outfol,
               reloctype,fileout,makedata,hypoinv,
               noiseswitch,noisediff,stdcc,stdct,
               nboot,nplot]
    return retlist


def readinfiles(noiseswitch,nclust,niter):
    """"
    Copy Over Jupyter Notebook Plotting
    -------------------
    Load Noise
    """
    tru1 = np.loadtxt('cat.loc.xy')
    for i in range(nclust):
        x = []
        if noiseswitch==1:
            x = np.loadtxt('noise.txt')

        if i==0:
            cal1=np.loadtxt('txtoutputs/cal_after_dtres_1_1.txt')
            cal4=np.loadtxt('txtoutputs/cal_after_dtres_1_%i.txt' % niter)
            dtdt1=np.loadtxt('txtoutputs/dtdt_after_dtres_1_1.txt')
            dtdt4=np.loadtxt('txtoutputs/dtdt_after_dtres_1_%i.txt' % niter)

            del1=np.loadtxt('txtoutputs/delta_source_1_1.txt')
            del4=np.loadtxt('txtoutputs/delta_source_1_%i.txt' % niter)

            abs0=np.loadtxt('txtoutputs/abs_source_1.txt')
            abs1=np.loadtxt('txtoutputs/abs_source_1_1.txt')
            abs4=np.loadtxt('txtoutputs/abs_source_1_%i.txt' % niter)

            # Load true event locations (should be all identical)
            # True locations SHOULD NOT CHANGE
            #tru1=np.loadtxt('txtoutputs/tru_source_1_8.txt')
        else:
            tmp=np.loadtxt('txtoutputs/cal_after_dtres_%i_1.txt' % i)
            cal1=np.append(cal1,tmp)
            tmp=np.loadtxt('txtoutputs/cal_after_dtres_%i_%i.txt' % (i,niter))
            cal4=np.append(cal4,tmp)
            tmp=np.loadtxt('txtoutputs/dtdt_after_dtres_%i_1.txt' % i)
            dtdt1=np.append(dtdt1,tmp)
            tmp=np.loadtxt('txtoutputs/dtdt_after_dtres_%i_%i.txt' % (i,niter))
            dtdt4=np.append(dtdt4,tmp)

            tmp=np.loadtxt('txtoutputs/delta_source_%i_1.txt' % i)
            del1=np.append(del1,tmp,axis=0)
            tmp=np.loadtxt('txtoutputs/delta_source_%i_%i.txt' % (i,niter))
            del4=np.append(del4,tmp,axis=0)

            tmp=np.loadtxt('txtoutputs/abs_source_%i.txt' % i)
            abs0=np.append(abs0,tmp,axis=0)
            tmp=np.loadtxt('txtoutputs/abs_source_%i_1.stxt' % i)
            abs1=np.append(abs1,tmp,axis=0)
            tmp=np.loadtxt('txtoutputs/abs_source_%i_%i.stxt' % (i,niter))
            abs4=np.append(abs4,tmp,axis=0)

            # Load true event locations (should be all identical)
            # True locations SHOULD NOT CHANGE
            #tmp=np.loadtxt('txtoutputs/tru_source_2_8.txt')
            #tru1=np.append(tru1,tmp)

        return [x,cal1,cal4,dtdt1,dtdt4,del1,del4,abs0,abs1,abs4,tru1]


#@profile(stream=open('mem_logs/readinnoise.mem','w+'))
def readinnoise(noiseswitch):
    x = []
    if noiseswitch==1:
        x = np.loadtxt('noise.txt')
    return x


#@profile(stream=open('mem_logs/statsOut.mem','w+'))
def statsOut(makedata,noiseswitch,reloctype,
             noisediff,stdcc,stdct,
             x,ncc,ndt,dtdt1,dtdt4,cal1,cal4,
             tru1,abs0,abs4,hypoinv):
    """
    Check Data Limitations
    """
    if makedata==1:
        print('\n\nSynthetic Data Test.')
    if noiseswitch==0:
        print('No Noise')
    elif noiseswitch==1:
        print('Noise Added.')
        if noisediff==0:
            print('Theory Noise Added (single value added to residuals)')
        elif noisediff==1:
            print('Measurement Noise Added (noise added to traveltimes, residual noise is noise difference)')
        print('Std. of Noise for CC: ',stdcc)
        print('Std. of Noise for CT: ',stdct)
        """
        Double Check if the Added Noise has the Correct STD
        ---
        If noisediff=0, should equal stdcc and stdct
        If noisediff=1 and reloctype=1 or 2, should be 2* stdcc or stdct
        If noisediff=1 and reloctype=3, should be 4* stdcc or stdct
        """
        print('Std. of Added Noise Array (CC):',np.std(x[:ncc]))
        print('Std. of Added Noise Array (CT):',np.std(x[ncc:]))
    """
    Residual Information Outputs
    """
    if reloctype==1:
        print('RMS Measured EvPair Res.: ',np.sqrt(np.mean(dtdt1*dtdt1)))
        print('RMS Cal Starting EvPair Res.: ',np.sqrt(np.mean(cal1*cal1)))
        print('RMS Cal Final EvPair Res.: ',np.sqrt(np.mean(cal4*cal4)))
    if reloctype==2:
        print('RMS Measured StPair Res.: ',np.sqrt(np.mean(dtdt1*dtdt1)))
        print('RMS Cal Starting StPair Res.: ',np.sqrt(np.mean(cal1*cal1)))
        print('RMS Cal Final StPair Res.: ',np.sqrt(np.mean(cal4*cal4)))
    if reloctype==3:
        print('RMS Measured DbPair Res.: ',np.sqrt(np.mean(dtdt1*dtdt1)))
        print('RMS Cal Starting DbPair Res.: ',np.sqrt(np.mean(cal1*cal1)))
        print('RMS Cal Final DbPair Res.: ',np.sqrt(np.mean(cal4*cal4)))
    print('RMS Starting Doub-Diff Res.: ',np.sqrt(np.mean((dtdt1-cal1)**2)))
    print('RMS Final Doub-Diff Res.: ',np.sqrt(np.mean((dtdt4-cal4)**2)))
    print('Percent Reduction in RMS Doub-Diff Res.:',
          100.*((np.sqrt(np.mean((dtdt4-cal4)**2))-np.sqrt(np.mean((dtdt1-cal1)**2)))/np.sqrt(np.mean((dtdt1-cal1)**2))))
    print('\nNo. of Data: ',ndt)
    """
    Calculated Euc. Dists.
    """
    ed_tru = np.sqrt(tru1[:,0]**2 + tru1[:,1]**2, + tru1[:,2]**2)
    ed_ini = np.sqrt(abs0[:,0]**2 + abs0[:,1]**2, + abs0[:,2]**2)
    ed_rel = np.sqrt(abs4[:,0]**2 + abs4[:,1]**2, + abs4[:,2]**2)
    nev = (tru1.shape)[0]
    """
    Absolute Location Error Outputs
    """
    #if hypoinv==1:
    print('\n\nAbsolute Location Errors:')
    cat_err = abs0[:,0:3] - tru1[:,0:3]
    cat_errx = np.sqrt(np.sum(cat_err[:,0]**2)/nev)
    cat_erry = np.sqrt(np.sum(cat_err[:,1]**2)/nev)
    cat_errz = np.sqrt(np.sum(cat_err[:,2]**2)/nev)
    cat_errt = np.sqrt(np.sum((abs0[:,3])**2)/nev)
    cat_errs = np.sqrt(np.sum((abs0[:,4])**2)/nev)
    cat_erred = np.sqrt(cat_errx**2 + cat_erry**2 + cat_errz**2)
    cat_eddists = np.sqrt(cat_err[:,0]**2 + cat_err[:,1]**2 + cat_err[:,2]**2)
    # Catalog Errors
    print('\nRMS Catalog X Error: ',cat_errx)
    print('RMS Catalog Y Error: ',cat_erry)
    print('RMS Catalog Z Error: ',cat_errz)
    if reloctype==1:
        print('RMS Catalog T Error: ',cat_errt)
    if reloctype==2:
        print('RMS Catalog Rel. St. Correc.: ',cat_errs)
    print('RMS Catalog Euc. Dist. Error: ',cat_erred)
    rel_err = abs4[:,0:3] - tru1[:,0:3]
    rel_errx = np.sqrt(np.sum(rel_err[:,0]**2)/nev)
    rel_erry = np.sqrt(np.sum(rel_err[:,1]**2)/nev)
    rel_errz = np.sqrt(np.sum(rel_err[:,2]**2)/nev)
    rel_errt = np.sqrt(np.sum((abs4[:,3])**2)/nev)
    rel_errs = np.sqrt(np.sum((abs4[:,4])**2)/nev)
    rel_erred = np.sqrt(rel_errx**2 + rel_erry**2 + rel_errz**2)
    rel_eddists = np.sqrt(rel_err[:,0]**2 + rel_err[:,1]**2 + rel_err[:,2]**2)
    # Relocation Errors
    print('\nRMS Relocation X Error: ',rel_errx)
    print('RMS Relocation Y Error: ',rel_erry)
    print('RMS Relocation Z Error: ',rel_errz)
    if reloctype==1:
        print('RMS Relocation T Error: ',rel_errt)
    if reloctype==2:
        print('RMS Relocation Rel. St. Correc.: ',rel_errs)
    print('RMS Relocation Euc. Dist. Error: ',rel_erred)
    # Euclidian Distance Errors
    print('\nMean Euc. Dist. Change from Catalog: ',np.mean(rel_eddists[:]-cat_eddists[:]))
    print('Percent Reduction in Error from Catalog (Euc. Dist. (X,Y,Z)): %2.4f %% (%f, %f, %f)' %
          ((rel_erred-cat_erred)/cat_erred*100,(rel_errx-cat_errx)/cat_errx*100,
           (rel_erry-cat_erry)/cat_erry*100,(rel_errz-cat_errz)/cat_errz*100))
    # elif hypoinv==0:
    #     print('\n\nAbsolute Location Change from Catalog:')
    #     rel_err = abs4[:,0:3] - abs0[:,0:3]
    #     rel_errx = np.sqrt(np.sum(rel_err[:,0]**2)/nev)
    #     rel_erry = np.sqrt(np.sum(rel_err[:,1]**2)/nev)
    #     rel_errz = np.sqrt(np.sum(rel_err[:,2]**2)/nev)
    #     rel_errt = np.sqrt(np.sum((abs4[:,3])**2)/nev)
    #     rel_errs = np.sqrt(np.sum((abs4[:,4])**2)/nev)
    #     rel_erred = np.sqrt(rel_errx**2 + rel_erry**2 + rel_errz**2)
    #     rel_eddists = np.sqrt(rel_err[:,0]**2 + rel_err[:,1]**2 + rel_err[:,2]**2)
    #     # Relocation Errors
    #     print('\nRMS Relocation X: ',rel_errx)
    #     print('RMS Relocation Y: ',rel_erry)
    #     print('RMS Relocation Z: ',rel_errz)
    #     if reloctype==1:
    #         print('RMS Relocation T: ',rel_errt)
    #     if reloctype==2:
    #         print('RMS Relocation Rel. St. Correc.: ',rel_errs)
    #     print('RMS Relocation Euc. Dist.: ',rel_erred)
    #     # Euclidian Distance Errors
    #     print('\nMean Euc. Dist. Change from Catalog: ',np.mean(rel_eddists[:]))

    return [ed_tru,ed_ini,ed_rel]


