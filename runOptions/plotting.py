import numpy as np
import scipy.stats as stats

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec


def resplot(plotstring,noiseswitch,reloctype,x,stdcc,stdct,
            dtdt1,dtdt4,cal1,cal4,noCC=False):
    """
    Calculate Residuals
    """
    res = dtdt4-cal4
    res_start = dtdt1-cal1
    if noCC:
        NDAT=len(res)
    else:
        NDAT=int((res.shape)[0]/2)

    """
    Main Plot
    """
    if noCC:
        fig, axs = plt.subplots(nrows=1,ncols=1,figsize=(5,5))#,sharey=True)
        nbins = 50
        #ccmin = 1.1*np.min(np.minimum(res[:NDAT],res_start[:NDAT]))
        #ccmax = 1.1*np.max(np.maximum(res[:NDAT],res_start[:NDAT]))
        #binscc = np.arange(ccmin,ccmax,(ccmax-ccmin)/nbins)
        catmin = 1.1*np.min(np.minimum(res[:NDAT],res_start[:NDAT]))
        catmax = 1.1*np.max(np.maximum(res[:NDAT],res_start[:NDAT]))
        binscat = np.arange(catmin,catmax,(catmax-catmin)/nbins)
        
        """
        Cat Residuals
        ---
        4 Parts:
            Histogram of noise added to data (black)
            Histogram of starting residuals (blue)
            Histogram of final residuals (red)
            Plot of gaussian w/ std of cc (green)
        ---
        """
        axs.set_title('catalog')
        if noiseswitch==1:
            axs.axvspan(-2*stdct,2*stdct,alpha=0.1,color='green',label='2x Std. of Noise')
            axs.axvspan(-stdct,stdct,alpha=0.2,color='green',label='Std. of Noise')
            axs.hist(x[:NDAT],binscat,color='k',density=False,histtype='step',alpha=0.75,linewidth=2,linestyle='dotted',label='Noise')
        axs.hist(res[:NDAT],binscat,color='r',density=False,histtype='step',alpha=0.75,linewidth=1,label='Final Res.')
        axs.hist(res_start[:NDAT],binscat,color='b',density=False,histtype='step',alpha=0.75,linewidth=1,label='Starting Res.')
        axs.yaxis.set_major_formatter(FormatStrFormatter('%2i'))
        axs.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        axs.tick_params(axis='both',direction='out',colors='k',grid_color='k',grid_alpha=0.5,
                           bottom=True,top=True,left=True,right=True)
        if noiseswitch==1:
            axs.set_xticks([-5*stdct,0,5*stdct])  
            #axs[1].set_xlim((-8*stdct,8*stdct))
        axs.set_xlabel('Residual (ms)')

    else:
        fig, axs = plt.subplots(nrows=1,ncols=2,figsize=(8,5),sharey=True)

        """
        CC Residuals
        ---
        4 Parts:
            Histogram of noise added to data (black)
            Histogram of starting residuals (blue)
            Histogram of final residuals (red)
            Plot of gaussian w/ std of cc (green)
        ---
        """
        nbins = 50
        ccmin = 1.1*np.min(np.minimum(res[:NDAT],res_start[:NDAT]))
        ccmax = 1.1*np.max(np.maximum(res[:NDAT],res_start[:NDAT]))
        binscc = np.arange(ccmin,ccmax,(ccmax-ccmin)/nbins)
        catmin = 1.1*np.min(np.minimum(res[NDAT:],res_start[NDAT:]))
        catmax = 1.1*np.max(np.maximum(res[NDAT:],res_start[NDAT:]))
        binscat = np.arange(catmin,catmax,(catmax-catmin)/nbins)


        axs[0].set_title('cross-corr')
        if noiseswitch==1:
            axs[0].axvspan(-2*stdcc,2*stdcc,alpha=0.1,color='green')
            axs[0].axvspan(-stdcc,stdcc,alpha=0.2,color='green')
            axs[0].hist(x[:NDAT],binscc,color='k',density=False,histtype='step',alpha=0.75,linewidth=2,linestyle='dotted')
        axs[0].hist(res[:NDAT],binscc,color='r',density=False,histtype='step',alpha=0.75,linewidth=1)
        axs[0].hist(res_start[:NDAT],binscc,color='b',density=False,histtype='step',alpha=0.75,linewidth=1)
        axs[0].yaxis.set_major_formatter(FormatStrFormatter('%2i'))
        axs[0].xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        axs[0].tick_params(axis='both',direction='out',colors='k',grid_color='k',grid_alpha=0.5,
                           bottom=True,top=True,left=True,right=True)
        if noiseswitch==1:
            axs[0].set_xticks([-5*stdcc,0,5*stdcc])  
            #axs[0].set_xlim((-8*stdcc,8*stdcc))
        axs[0].set_ylabel('Counts')
        axs[0].set_xlabel('Residual (ms)')

        """
        Cat Residuals
        ---
        4 Parts:
            Histogram of noise added to data (black)
            Histogram of starting residuals (blue)
            Histogram of final residuals (red)
            Plot of gaussian w/ std of cc (green)
        ---
        """
        axs[1].set_title('catalog')
        if noiseswitch==1:
            axs[1].axvspan(-2*stdct,2*stdct,alpha=0.1,color='green',label='2x Std. of Noise')
            axs[1].axvspan(-stdct,stdct,alpha=0.2,color='green',label='Std. of Noise')
            axs[1].hist(x[NDAT:],binscat,color='k',density=False,histtype='step',alpha=0.75,linewidth=2,linestyle='dotted',label='Noise')
        axs[1].hist(res[NDAT:],binscat,color='r',density=False,histtype='step',alpha=0.75,linewidth=1,label='Final Res.')
        axs[1].hist(res_start[NDAT:],binscat,color='b',density=False,histtype='step',alpha=0.75,linewidth=1,label='Starting Res.')
        axs[1].yaxis.set_major_formatter(FormatStrFormatter('%2i'))
        axs[1].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        axs[1].tick_params(axis='both',direction='out',colors='k',grid_color='k',grid_alpha=0.5,
                           bottom=True,top=True,left=True,right=True)
        if noiseswitch==1:
            axs[1].set_xticks([-5*stdct,0,5*stdct])  
            #axs[1].set_xlim((-8*stdct,8*stdct))
        axs[1].set_xlabel('Residual (ms)')

    """
    Scale and Add Legend
    """
    plt.locator_params(axis='y', nbins=5)
    plt.subplots_adjust(right=0.65)
    plt.legend(bbox_to_anchor=(1.25,0.5), loc="center left", borderaxespad=0)
    """
    Save Figure To File
    """
    fig.savefig('%s/fig_resnoise_hist_%i.png' % (plotstring,reloctype))


    plt.close('all')
    return None


def locplots(plotstring,reloctype,nev,
             tru1,abs0,abs1,abs4,
             ed_tru,ed_ini,ed_rel,hypoinv):

    """
    One Big Plot divided into 2 subplots 
        (a) Absolute Locations
        (b) Absolute Locations - True Locations
    ---
    Each of the 2 subplots is further divided into further
    subplots (i) the Map View and (ii) the Depth View
    ---
    Includes subfigures function which requires matplotlib v3.4.2
    """

    """
    Initiate Main Plot
    """
    fig = plt.figure(constrained_layout=True,figsize=(8,8))
    outer = fig.subfigures(nrows=2,ncols=1)

    """
    Initiate Absolute Location Plots
    """
    outer[0].suptitle('Absolute Event Locations')
    inner1 = outer[0].subplots(nrows=1,ncols=2)

    # Subplot 1: Map view absolute locations with iteration
    inner1[0].set_title('Map View')
    for i in range(nev):
        if hypoinv==1:
            inner1[0].plot(tru1[i,0],tru1[i,1],'xb',markeredgewidth=3) # Need to code in plotting all iterations
        inner1[0].plot(abs0[i,0],abs0[i,1],'or')
        #inner1[0].plot(abs2[i,0],abs2[i,1],'xr')
        #inner1[0].plot(abs3[i,0],abs3[i,1],'xg')
        inner1[0].plot(abs4[i,0],abs4[i,1],'ok')
    inner1[0].axis('equal')
    inner1[0].grid(True)
    inner1[0].set_xlabel('X (m)')
    inner1[0].set_ylabel('Y (m)')
    inner1[0].locator_params(axis='y', nbins=5)
    inner1[0].locator_params(axis='x', nbins=3)
    inner1[0].set_ylim(ymin=-1.1*abs(max(inner1[0].get_ylim(),key=abs)),
                       ymax=1.1*abs(max(inner1[0].get_ylim(),key=abs)))
    inner1[0].set_xlim(xmin=-1.1*abs(max(inner1[0].get_xlim(),key=abs)),
                       xmax=1.1*abs(max(inner1[0].get_xlim(),key=abs)))

    # Subplot 2: Depth view absolute locations with iteration
    inner1[1].set_title('Depth View')
    for i in range(nev):
        if i==0:
            if hypoinv==1:
                inner1[1].plot(tru1[i,1],tru1[i,2],'xb',markeredgewidth=3,label='True Location')
            inner1[1].plot(abs0[i,1],abs0[i,2],'or',label='Initial Catalog')
            inner1[1].plot(abs4[i,1],abs4[i,2],'ok',label='Final Relocation')
        if hypoinv==1:
            inner1[1].plot(tru1[i,1],tru1[i,2],'xb',markeredgewidth=3)
        inner1[1].plot(abs0[i,1],abs0[i,2],'or')
        #inner1[1].plot(abs2[i,1],abs2[i,2],'xr')
        #inner1[1].plot(abs3[i,1],abs3[i,2],'xg')
        inner1[1].plot(abs4[i,1],abs4[i,2],'ok')
    inner1[1].invert_yaxis()
    inner1[1].axis('equal')
    inner1[1].grid(True)
    inner1[1].set_xlabel('Y (m)')
    inner1[1].set_ylabel('Z (m)')
    inner1[1].locator_params(axis='y', nbins=5)
    inner1[1].locator_params(axis='x', nbins=3)
    inner1[1].set_ylim(ymin=-1.1*abs(max(inner1[1].get_ylim(),key=abs)),
                       ymax=1.1*abs(max(inner1[1].get_ylim(),key=abs)))
    inner1[1].set_xlim(xmin=-1.1*abs(max(inner1[1].get_xlim(),key=abs)),
                       xmax=1.1*abs(max(inner1[1].get_xlim(),key=abs)))
    inner1[1].legend(loc='upper right')

    """
    Initiate Relative to True Location Plots
    """
    outer[1].suptitle('Abs. Location Error')
    inner2 = outer[1].subplots(nrows=1,ncols=1)

    for i in range(nev):
        if hypoinv==1:
            x = [ed_tru[i],ed_tru[i]]
            y = [ed_ini[i]-ed_tru[i],ed_rel[i]-ed_tru[i]]
            inner2.plot(x,y,'k:',linewidth=0.5)
        if hypoinv==1:
            inner2.plot(ed_tru,ed_ini-ed_tru,'or')
        inner2.plot(ed_tru,ed_rel-ed_tru,'ok')
    inner2.grid(True)
    if hypoinv==1:
        inner2.set_xlabel('True Euc. Dist. From Cluster Centroid (m)')
        inner2.set_ylabel('Error from True (m)')
    else:
        inner2.set_xlabel('Catalog Euc. Dist. From Cluster Centroid (m)')
        inner2.set_ylabel('Change from Catalog (m)')

    """
    # Change to true locations map view
    inner2[0].set_title('Map View')
    for i in range(nev):
        inner2[0].plot(abs0[i,0]-tru1[i,0],abs0[i,1]-tru1[i,1],'or')
        #inner2[0].plot(abs2[i,0]-tru1[i,0],abs2[i,1]-tru1[i,1],'xr')
        #inner2[0].plot(abs3[i,0]-tru1[i,0],abs3[i,1]-tru1[i,1],'xg')
        inner2[0].plot(abs4[i,0]-tru1[i,0],abs4[i,1]-tru1[i,1],'ob')
    inner2[0].grid(True)
    inner2[0].axis('equal')
    inner2[0].set_xlabel('X (m)')
    inner2[0].set_ylabel('Y (m)')
    inner2[0].locator_params(axis='y', nbins=5)
    inner2[0].locator_params(axis='x', nbins=3)
    inner2[0].set_ylim(ymin=-1.5*abs(max(inner2[0].get_ylim(),key=abs)),
                       ymax=1.5*abs(max(inner2[0].get_ylim(),key=abs)))
    inner2[0].set_xlim(xmin=-1.5*abs(max(inner2[0].get_xlim(),key=abs)),
                       xmax=1.5*abs(max(inner2[0].get_xlim(),key=abs)))

    # Change to depth locations depth view
    inner2[1].set_title('Depth View')
    for i in range(nev):
        inner2[1].plot(abs0[i,1]-tru1[i,1],abs0[i,2]-tru1[i,2],'or')
        #inner2[1].plot(abs2[i,1]-tru1[i,1],abs2[i,2]-tru1[i,2],'xr')
        #inner2[1].plot(abs3[i,1]-tru1[i,1],abs3[i,2]-tru1[i,2],'xg')
        inner2[1].plot(abs4[i,1]-tru1[i,1],abs4[i,2]-tru1[i,2],'ob')
    inner2[1].invert_yaxis()
    inner2[1].grid(True)
    inner2[1].axis('equal')
    inner2[1].set_xlabel('Y (m)')
    inner2[1].set_ylabel('Z (m)')
    inner2[1].locator_params(axis='y', nbins=5)
    inner2[1].locator_params(axis='x', nbins=3)
    #inner2[1].set_xticks([-0.8*abs(max(inner2[1].get_xlim(),key=abs)),0,0.8*abs(max(inner2[1].get_xlim(),key=abs))])
    inner2[1].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    inner2[1].set_ylim(ymin=-1.5*abs(max(inner2[1].get_ylim(),key=abs)),
                       ymax=1.5*abs(max(inner2[1].get_ylim(),key=abs)))
    inner2[1].set_xlim(xmin=-1.5*abs(max(inner2[1].get_xlim(),key=abs)),
                       xmax=1.5*abs(max(inner2[1].get_xlim(),key=abs)))
    """

    """
    Save Figure To File
    """
    fig.savefig('%s/fig_absolute_location_errors_%i.png' % (plotstring,reloctype))

    plt.close('all')
    return None


def evpairsplot(nev,plotstring,reloctype,tru1,abs0,abs4,hypoinv):
    """
    Map Chosen Event-Pairs
    evpairs = np.loadtxt('event_pairing.txt')
    # Plot event-pair web
    fig, axs = plt.subplots(nrows=1, ncols=2)
    fig.suptitle('Event-Pair Cluster Web')
    axs[0].set_title('Map-View')
    axs[0].plot(tru1[:,0],tru1[:,1],'ok')
    for i in range(0,nev-1):
        for j in range(i,nev):
            if evpairs[i,j]==0:
                axs[0].plot([tru1[i,0],tru1[j,0]],[tru1[i,1],tru1[j,1]],'--',color='b')
                axs[1].set_title('Depth-View')
                axs[1].plot(tru1[:,1],tru1[:,2],'ok')
    for i in range(0,nev-1):
        for j in range(i,nev):
            if evpairs[i,j]==0:
                axs[1].plot([tru1[i,1],tru1[j,1]],[tru1[i,2],tru1[j,2]],'--',color='b')
                axs[1].invert_yaxis()
                plt.savefig('%s/fig_eventpairs_%i.png' % (plotstring,reloctype)) #,iboot))
    """
    """
    Relative Location Mesurements
    """
    # First find all possible event pairs
    evpairs = np.zeros((nev**2,2),dtype='int')
    nevp = 0
    for i1 in range(nev-1):
        for i2 in range(i1,nev):
            evpairs[nevp,0] = i1
            evpairs[nevp,1] = i2
            nevp += 1
    evpairs = evpairs[:nevp,:]

    # Calculate Offsets (Euc. Dist.)
    tru_offs = np.zeros((nevp,4),dtype='float')
    cat_offs = np.zeros((nevp,4),dtype='float')
    rel_offs = np.zeros((nevp,4),dtype='float')
    for ievp in range(nevp):
        # True Offsets
        offx = np.abs(tru1[evpairs[ievp,0],0]-tru1[evpairs[ievp,1],0])
        offy = np.abs(tru1[evpairs[ievp,0],1]-tru1[evpairs[ievp,1],1])
        offz = np.abs(tru1[evpairs[ievp,0],2]-tru1[evpairs[ievp,1],2])
        offs = np.sqrt(offx**2 + offy**2 + offz**2)
        tru_offs[ievp,0] = offx
        tru_offs[ievp,1] = offy
        tru_offs[ievp,2] = offz
        tru_offs[ievp,3] = offs

        # Catalog Offsets
        offx = np.abs(abs0[evpairs[ievp,0],0]-abs0[evpairs[ievp,1],0])
        offy = np.abs(abs0[evpairs[ievp,0],1]-abs0[evpairs[ievp,1],1])
        offz = np.abs(abs0[evpairs[ievp,0],2]-abs0[evpairs[ievp,1],2])
        offs = np.sqrt(offx**2 + offy**2 + offz**2)
        cat_offs[ievp,0] = offx
        cat_offs[ievp,1] = offy
        cat_offs[ievp,2] = offz
        cat_offs[ievp,3] = offs

        # Relative Offsets
        offx = np.abs(abs4[evpairs[ievp,0],0]-abs4[evpairs[ievp,1],0])
        offy = np.abs(abs4[evpairs[ievp,0],1]-abs4[evpairs[ievp,1],1])
        offz = np.abs(abs4[evpairs[ievp,0],2]-abs4[evpairs[ievp,1],2])
        offs = np.sqrt(offx**2 + offy**2 + offz**2)
        rel_offs[ievp,0] = offx
        rel_offs[ievp,1] = offy
        rel_offs[ievp,2] = offz
        rel_offs[ievp,3] = offs

    # Sort by True Location Euc. Dist From Cluster Centroid (i.e. 0,0,0)
    ind = np.argsort(tru_offs[:,3])
    tru_offs = tru_offs[ind,:]
    cat_offs = cat_offs[ind,:]
    rel_offs = rel_offs[ind,:]
    evpairs = evpairs[ind,:]

    if hypoinv==1:
        cat_everr = cat_offs[:,0:3] - tru_offs[:,0:3]
        cat_everrx = np.sqrt(np.sum(cat_everr[:,0]**2)/nevp)
        cat_everry = np.sqrt(np.sum(cat_everr[:,1]**2)/nevp)
        cat_everrz = np.sqrt(np.sum(cat_everr[:,2]**2)/nevp)
        cat_everred = np.sqrt(cat_everrx**2 + cat_everry**2 + cat_everrz**2)
        cat_eddists = np.sqrt(cat_everr[:,0]**2 + cat_everr[:,1]**2 + cat_everr[:,2]**2)

    rel_everr = rel_offs[:,0:3] - tru_offs[:,0:3]
    rel_everrx = np.sqrt(np.sum(rel_everr[:,0]**2)/nevp)
    rel_everry = np.sqrt(np.sum(rel_everr[:,1]**2)/nevp)
    rel_everrz = np.sqrt(np.sum(rel_everr[:,2]**2)/nevp)
    rel_everred = np.sqrt(rel_everrx**2 + rel_everry**2 + rel_everrz**2)
    rel_eddists = np.sqrt(rel_everr[:,0]**2 + rel_everr[:,1]**2 + rel_everr[:,2]**2)

    if hypoinv==1:
        print('\n\nRel. Location Errors:')
        print('\nRMS Catalog X Error:',cat_everrx)
        print('RMS Catalog Y Error:',cat_everry)
        print('RMS Catalog Z Error:',cat_everrz)
        print('RMS Catalog Euc. Dist. Error:',cat_everred)

        print('\nRMS Reloc X Error: ',rel_everrx)
        print('RMS Reloc Y Error: ',rel_everry)
        print('RMS Reloc Z Error: ',rel_everrz)
        print('RMS Reloc Euc. Dist. Error: ',rel_everred)
    else:
        print('\n\nRel. Location Change from Catalog:')
        print('\nRMS Reloc X: ',rel_everrx)
        print('RMS Reloc Y: ',rel_everry)
        print('RMS Reloc Z: ',rel_everrz)
        print('RMS Reloc Euc. Dist.: ',rel_everred)

    if hypoinv==1:
        print('Mean Euc. Dist. Change from Catalog: ',
              np.mean(rel_eddists[:]-cat_eddists[:]))
        print('Percent Reduction in Error from Catalog (Euc. Dist. (X,Y,Z)): %2.4f %% (%f, %f, %f)\n\n' %
              ((rel_everred-cat_everred)/cat_everred*100,(rel_everrx-cat_everrx)/cat_everrx*100,
               (rel_everry-cat_everry)/cat_everry*100,(rel_everrz-cat_everrz)/cat_everrz*100))

    #Plot First Ev Pair Figure: 
    #      Subplot 1 (Large): Euclidian Dist True (x) by Diff. from True Euc. Dist. (y)
    #                         for both the catalog and final relocation
    #      Subplots 2-4 (Small): Subplots for the X,Y,Z distances with true dist (x) and
    #                            error difference on vertical (both catalog and relocation)    

    # fig = plt.figure(constrained_layout=True,figsize=(10,10))
    # outer = fig.subfigures(nrows=2,ncols=1)

    # # Subplot 1: Eudlidian Distance
    # inner1 = outer[0].subplots(nrows=1,ncols=1)
    # for i in range(nevp):
    #     if hypoinv==1:
    #         x = [tru_offs[i,3],tru_offs[i,3]]
    #         y = [cat_offs[i,3]-tru_offs[i,3], rel_offs[i,3]-tru_offs[i,3]]
        
    #     if i==0:
    #         if hypoinv==1:
    #             inner1.plot(tru_offs[i,3],cat_offs[i,3]-tru_offs[i,3],'or',label='Catalog Error')
    #         inner1.plot(tru_offs[i,3],rel_offs[i,3]-tru_offs[i,3],'ok',label='Relocation Error')
        
    #     if hypoinv==1:
    #         inner1.plot(x,y,'k:',linewidth=0.5)
    #         inner1.plot(tru_offs[i,3],cat_offs[i,3]-tru_offs[i,3],'or')
    #     inner1.plot(tru_offs[i,3],rel_offs[i,3]-tru_offs[i,3],'ok')
    # inner1.grid(True)
    # #import pdb; pdb.set_trace()
    # if hypoinv==1:
    #     inner1.set_xlabel('True Interevent Euc. Dist. (m)')
    #     inner1.set_ylabel('Error from True (m)')
    # else:
    #     inner1.set_xlabel('Catalog Interevent Euc. Dist. (m)')
    #     inner1.set_ylabel('Change from Catalog (m)')
    # inner1.legend(loc='lower right')


    # inner = outer[1].subplots(nrows=1,ncols=3,sharey=True)
    # # Subplot 2: X Errors
    # for i in range(nevp):
    #     if hypoinv==1:
    #         inner[0].plot(tru_offs[i,0],cat_offs[i,0]-tru_offs[i,0],'or')
    #     inner[0].plot(tru_offs[i,0],rel_offs[i,0]-tru_offs[i,0],'ok')
    # inner[0].grid(True)
    # if hypoinv==1:
    #     inner[0].set_xlabel('True Interevent X Dist. (m)')
    #     inner[0].set_ylabel('Error from True (m)')
    # else:
    #     inner[0].set_xlabel('Catalog Interevent X Dist. (m)')
    #     inner[0].set_ylabel('Change from Catalog (m)')

    # # Subplot 3: Y Errors
    # for i in range(nevp):
    #     if hypoinv==1:
    #         inner[1].plot(tru_offs[i,1],cat_offs[i,1]-tru_offs[i,1],'or')
    #     inner[1].plot(tru_offs[i,1],rel_offs[i,1]-tru_offs[i,1],'ok')
    # inner[1].grid(True)
    # if hypoinv==1:
    #     inner[1].set_xlabel('True Interevent Y Dist. (m)')
    # else:
    #     inner[1].set_xlabel('Catalog Interevent Y Dist. (m)')

    # # Subplot 4: Z Errors
    # for i in range(nevp):
    #     if hypoinv==1:
    #         inner[2].plot(tru_offs[i,2],cat_offs[i,2]-tru_offs[i,2],'or')
    #     inner[2].plot(tru_offs[i,2],rel_offs[i,2]-tru_offs[i,2],'ok')
    # inner[2].grid(True)
    # if hypoinv==1:
    #     inner[2].set_xlabel('True Interevent Z Dist. (m)')
    # else:
    #     inner[2].set_xlabel('Catalog Interevent Z Dist. (m)')

    # fig.savefig('%s/fig_relative_location_errors_%i.png' % (plotstring,reloctype))

    # plt.close('all')
    return None




