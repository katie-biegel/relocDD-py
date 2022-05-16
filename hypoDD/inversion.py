import numpy as np
from datetime import datetime
import scipy.linalg as la
import scipy.sparse.linalg as sp
from scipy.sparse import csr_matrix
from scipy.sparse import csc_matrix
import sys
import matplotlib.pyplot as plt
#from memory_profiler import profile

# Other hypoDD functions
from hypoDD.hypoDD_functions import resstat
from methodTypes.eventPair.hypoDD_ev import svdg_ev,svdunpack_ev,lsqrprep_ev,lsqrunpack_ev
from methodTypes.stationPair.hypoDD_stat import svdg_st,svdunpack_st,lsqrprep_st,lsqrunpack_st
from methodTypes.doublePair.hypoDD_doub import svdg_doub,svdunpack_doub,lsqrprep_doub,lsqrunpack_doub

""""""
############################### KB Note: There are some notable biasing/uncertainty issues
####                                     with the SVD method.
####                                     We highly recommend using LSQR inversion as the 
####                                     SVD generally is shifted and overfits data.
####                                     Additionally, LSQR with bootstrapping has more 
####                                     robust uncertainty analysis and can be used for 
####                                     larger problems than SVD.
####                                     It does work however and can be used.
""""""

def svd(log,reloctype,iteri,ndt,nev,nsrc,damp,mod_ratio,idata,
        ev_cusp,src_cusp,dt_res,dt_wt,dt_ista1,dt_ista2,dt_ic1,dt_ic2,
        exav,eyav,ezav,etav,esav,dxav,dyav,dzav,dtav,dsav,
        rms_cc,rms_ct,rms_cc0,rms_ct0,rms_ccold,rms_ctold,rms_cc0old,rms_ct0old,
        tmp_xp,tmp_yp,tmp_zp,dt_idx):
    """
    Translation of lsfit_svd.f
    :::
    PARAMETERS:
    log (file obj) ---- Log file
    reloctype (int) ---- Double-difference pairing type
    iteri (int) ---- Inversion iteration index
    ndt (int) ---- No. of data
    nev (int) ---- No. of events
    nsrc (int) ---- No. of sources
    damp (float) ---- Dampig value for inversion
    mod_ratio (float) ---- VPVS ratio
    idata (int) ---- Type of data switch
    ev_cusp[nev] (int array) ---- Array holding event ids
    src_cusp[nev] (int array) ---- Array holding src location arrays
    dt_res[ndt] (float array) ---- Array holding data residuals
    dt_wt[ndt] (float array) ---- Array holding data weights
    dt_ista1[ndt] (int array) ---- Array holding station indexes
    dt_ista2[ndt] (int array) ---- Array holding station indexes
    dt_ic1[ndt] (int array) ---- Array holding ev1 indexes
    dt_ic2[ndt] (int array) ---- Array holding ev2 indexes
    exav (float) ---- Avg. x error
    eyav (float) ---- Avg. y error
    ezav (float) ---- Avg. z error
    etav (float) ---- Avg. t error
    esav (float) ---- Avg. station correction error
    dxav (float) ---- Change in x value
    dyav (float) ---- Change in y value
    dzav (float) ---- Change in z value
    dtav (float) ---- Change in t value
    dsav (float) ---- Change in station correction value
    rms_cc (float) --- RMS for cc residuals
    rms_ct (float) --- RMS for ct residuals
    rms_cc0
    rms_ct0
    rms_ccold (float) --- RMS for previous iteration cc residuals
    rms_ctold (float) --- RMS for previous iteration ct residuals
    rms_cc0old
    rms_ct0old
    tmp_xp[NSTA,NEV] (float array) --- X partial derivatives from ray tracing
    tmp_yp[NSTA,NEV] (float array) --- Y partial derivatives from ray tracing
    tmp_zp[NSTA,NEV] (float array) --- Z partial derivatives from ray tracing
    dt_idx[ndt] (float array) ---- Array holding data type indexes
    :::
    RETURNS (ALL DEFINED IN PARAMETERS)
    retlist (list) --- List of returned variables
    :::
    """

    """
    Set up full G inversion matrix    
    ######
    The G matrix holds the partial derivatives with respect to
    the hypocentral model paramters.
    """

    if reloctype==1:
        g = svdg_ev(ndt,nev,nsrc,dt_ic1,dt_ic2,dt_ista1,tmp_xp,tmp_yp,tmp_zp,
                    dt_idx,mod_ratio)
    elif reloctype==2:
        g = svdg_st(ndt,nev,nsrc,dt_ic1,dt_ista1,dt_ista2,tmp_xp,tmp_yp,tmp_zp,
                    dt_idx,mod_ratio)
    elif reloctype==3:
        g = svdg_doub(ndt,nev,nsrc,dt_ic1,dt_ic2,dt_ista1,dt_ista2,tmp_xp,tmp_yp,tmp_zp,
                      dt_idx,mod_ratio)
    
    """
    Now define inversion weights and data array
    """
    # Define arrays
    wtinv = np.zeros((ndt))
    if reloctype==2:
        wt = np.zeros((ndt+3))
        d = np.zeros((ndt+3))
    else:
        wt = np.zeros((ndt+4))
        d = np.zeros((ndt+4))

    # Set Values
    wt[0:ndt] = np.copy(dt_wt[0:ndt])   # First copy data weights
    d[0:ndt] = dt_res[0:ndt]*1000.*wt[0:ndt] # Copy and weight data array 
                                             # d array now in ms
    # Calculate the inverse weight matrix
    wtinv = np.copy(wt[0:ndt])
    wtinv = np.where(wtinv==0.,1.,1./wtinv)
    
    """
    Now weight G matrix
    ###################
    Multiply 2d matrix by 1d
    """
    g[0:ndt] = wt[0:ndt].reshape((ndt,1))*g[0:ndt,:]
    # Add four extra rows to make mean shift zero
    # This should make the design matrix non-singular
    wt[ndt:] = 1.
    if reloctype==3:
        for i in range(0,3):
            for j in range(0,nev):
                g[ndt+i,3*j-3+i] = 1.0
        nndt = ndt+3
    else:
        for i in range(0,4):
            for j in range(0,nev):
                g[ndt+i,4*j-4+i] = 1.0
        nndt = ndt+4

    """
    Scale/normalize columns
    """
    if reloctype==3:
        norm = np.zeros((3*nev))
    else:
        norm = np.zeros((4*nev))
    norm = np.sqrt(np.sum(g*g,axis=0))/nndt
    g = g/norm

    """"
    Run SVD Here
    """
    log.write('~ singular value decomposition of G %6i x %6i matrix ... \n\n ' % (nndt,nev*4))
    [u,q,v] = np.linalg.svd(g,full_matrices=False)
    v = np.transpose(v)

    """
    Check for negative singular values or
    for singular values near zero
    """
    log.write('~ backsubstitution ... \n')
    izero = 0
    qmax = np.amax(q)

    # Neg. sing.value check
    if np.amin(q)<0:
        raise Exception('Fatal Error: (svd: neg singular value)')

    # Replace small sing. values with zeros
    qmin = qmax*0.0000001
    izero = (q<qmin).sum()
    q = np.where(q<qmin,0.,q)

    # If there are sing. values close to zero note them
    if izero>0:
        log.write('>>> %3i singular values close/equal to zero \n' % izero)
        print('>>> %3i singular values close/equal to zero' % izero)

    """
    Backsubstitution 
    ################
    Get x' from Gx=d: x=v*diag(1/q)*t(u)*d
    Compute diag(1/q)*t(u)*d
    """
    if reloctype==3:
        tmp = np.zeros(3*nev)
        for j in range(0,3*nev):
            if q[j]!=0.:
                tmp[j] = np.sum(u[:,j]*d/q[j])
    else:
        tmp = np.zeros(4*nev)
        for j in range(0,4*nev):
            if q[j]!=0.:
                tmp[j] = np.sum(u[:,j]*d/q[j])
    x = np.sum(v*tmp,axis=1)

    # Rescale model vector and G
    x = x/norm
    g = g*norm

    # Unweight G matrix
    if reloctype==3:
        for i in range(0,3*nev):
            g[0:ndt,i] = wtinv*g[0:ndt,i]
    else:
        for i in range(0,4*nev):
            g[0:ndt,i] = wtinv*g[0:ndt,i]

    """
    Predict data 
    ############
    From g*x' and get residuals
    """
    dd = np.sum(g[0:ndt,:]*x,axis=1)
    dt_res[0:ndt] = dt_res[0:ndt] - dd/1000.    # Back in seconds
    
    """"
    Get covariance matrix
    #############
    cvm = v*(1/q**2)*vt
    """
    q2 = 1./(q*q)
    if reloctype==3:
        cvm = np.zeros((3*nev,3*nev),dtype='float')
        for i in range(0,3*nev):
            for j in range(0,i+1):
                s = np.sum(np.where(q!=0.,v[i,:]*v[j,:]*q2,0.))
                cvm[i,j] = s
                cvm[j,i] = s
    else:
        cvm = np.zeros((4*nev,4*nev),dtype='float')
        for i in range(0,4*nev):
            for j in range(0,i+1):
                s = np.sum(np.where(q!=0.,v[i,:]*v[j,:]*q2,0.))
                cvm[i,j] = s
                cvm[j,i] = s

    """
    Get residual statistics
    """
    resvar1 = 0.
    [rms_cc,rms_ct,rms_cc0,rms_ct0,rms_ccold,rms_ctold,
     rms_cc0old,rms_ct0old,resvar1] = resstat(log,idata,ndt,nev,dt_res,dt_wt,dt_idx,
                                              rms_cc,rms_ct,rms_cc0,rms_ct0,
                                              rms_ccold,rms_ctold,rms_cc0old,rms_ct0old,
                                              resvar1)

    """
    Errors for 95% confidence level
    """
    factor=1.96*np.sqrt(resvar1)

    if reloctype==1:
        [src_dx,src_dy,src_dz,src_dt,
         src_ex,src_ey,src_ez,src_et,src_cusp,
         exavold,eyavold,ezavold,etavold,
         dxavold,dyavold,dzavold,dtavold,
         exav,eyav,ezav,etav,
         dxav,dyav,dzav,dtav] = svdunpack_ev(log,nev,iteri,cvm,factor,norm,x,
                                             ev_cusp,exav,eyav,ezav,etav,
                                             dxav,dyav,dzav,dtav)
        src_ds = np.zeros(nev,dtype='float')
        src_es = np.zeros(nev,dtype='float')
        esav = float(0.)
        dsav = float(0.)
    elif reloctype==2:
        [src_dx,src_dy,src_dz,src_ds,
         src_ex,src_ey,src_ez,src_es,src_cusp,
         exavold,eyavold,ezavold,esavold,
         dxavold,dyavold,dzavold,dsavold,
         exav,eyav,ezav,esav,
         dxav,dyav,dzav,dsav] = svdunpack_st(log,nev,iteri,cvm,factor,norm,x,
                                             ev_cusp,exav,eyav,ezav,esav,
                                             dxav,dyav,dzav,dsav)
        src_dt = np.zeros(nev,dtype='float')
        src_et = np.zeros(nev,dtype='float')
        etav = float(0.)
        dtav = float(0.)
    elif reloctype==3:
        [src_dx,src_dy,src_dz,
         src_ex,src_ey,src_ez,src_cusp,
         exavold,eyavold,ezavold,
         dxavold,dyavold,dzavold,
         exav,eyav,ezav,
         dxav,dyav,dzav] = svdunpack_doub(log,nev,iteri,cvm,factor,norm,x,
                                          ev_cusp,exav,eyav,ezav,
                                          dxav,dyav,dzav)
        src_dt = np.zeros(nev,dtype='float')
        src_et = np.zeros(nev,dtype='float')
        etav = float(0.)
        dtav = float(0.)
        src_ds = np.zeros(nev,dtype='float')
        src_es = np.zeros(nev,dtype='float')
        esav = float(0.)
        dsav = float(0.)

    """
    Return updated locations and inversion statistics
    """
    retlist = [src_cusp,src_dx,src_dy,src_dz,src_dt,src_ds,
               src_ex,src_ey,src_ez,src_et,src_es,
               exav,eyav,ezav,etav,esav,
               dxav,dyav,dzav,dtav,dsav,
               rms_cc,rms_ct,rms_cc0,rms_ct0,
               rms_ccold,rms_ctold,rms_cc0old,rms_ct0old]
    return  retlist


#@profile(stream=open('mem_logs/lsqr_test.mem','w+'))
def lsqr(log,reloctype,iteri,ndt,nev,nsrc,damp,mod_ratio, 
         idata,ev_cusp,src_cusp,dt_res,dt_wt,dt_ista1,dt_ista2,dt_ic1,dt_ic2,
         exav,eyav,ezav,etav,esav,dxav,dyav,dzav,dtav,dsav,
         rms_cc,rms_ct,rms_cc0,rms_ct0,rms_ccold,rms_ctold,rms_cc0old,rms_ct0old,
         tmp_xp,tmp_yp,tmp_zp,dt_idx,normvar):
    """
    Translation of lsfit_lsqr.f
    :::
    PARAMETERS:
    log (file obj) ---- Log file
    reloctype (int) --- Double-difference pairing type
    iteri (int) ---- Inversion iteration index
    ndt (int) ---- No. of data
    nev (int) ---- No. of events
    nsrc (int) ---- No. of events
    damp (float) ---- Dampig value for inversion
    mod_ratio (float) ---- VPVS ratio
    idata (int) ---- Type of data switch
    ev_cusp[nev] (int array) ---- Array holding event ids
    src_cusp[nev] (int array) ---- Array holding src location arrays
    dt_res[ndt] (float array) ---- Array holding data residuals
    dt_wt[ndt] (float array) ---- Array holding data weights
    dt_ista1[ndt] (int array) ---- Array holding station indexes
    dt_ista2[ndt] (int array) ---- Array holding station indexes
    dt_ic1[ndt] (int array) ---- Array holding ev1 indexes
    dt_ic2[ndt] (int array) ---- Array holding ev2 indexes
    exav (float) ---- Avg. x error
    eyav (float) ---- Avg. y error
    ezav (float) ---- Avg. z error
    etav (float) ---- Avg. t error
    esav (float) ---- Avg. st. corr. term error
    dxav (float) ---- Change in x value
    dyav (float) ---- Change in y value
    dzav (float) ---- Change in z value
    dtav (float) ---- Change in t value
    dsav (float) ---- Change in st. corr. term
    rms_cc (float) --- RMS for cc residuals
    rms_ct (float) --- RMS for ct residuals
    rms_cc0
    rms_ct0
    rms_ccold (float) --- RMS for previous iteration cc residuals
    rms_ctold (float) --- RMS for previous iteration ct residuals
    rms_cc0old
    rms_ct0old
    tmp_xp[NSTA,NEV] (float array) --- X partial derivatives from ray tracing
    tmp_yp[NSTA,NEV] (float array) --- Y partial derivatives from ray tracing
    tmp_zp[NSTA,NEV] (float array) --- Z partial derivatives from ray tracing
    dt_idx[ndt] (float array) ---- Array holding data type indexes
    normvar (float) ----
    :::
    RETURNS (ALL DEFINED IN PARAMETERS)
    retlist (list) --- List of returned variables
    :::
    """
    # Define Common Variables
    wt = np.copy(dt_wt[0:ndt])
    wtinv = np.zeros(ndt)
    d = dt_res[0:ndt]*1000.*wt[0:ndt]   # Weighted data in ms again
    wtinv = np.where(wt==0.,1.,1./wt)
    """
    First define G matrix
    ############
    Using sparse matrix objects to conserve memory
    """
    if reloctype==1:
        [nar,nndt,row_i,col_i,rw] = lsqrprep_ev(log,ndt,nsrc,dt_idx,tmp_xp,tmp_yp,tmp_zp,
                                                dt_ista1,dt_ic1,dt_ic2,mod_ratio,wt)
    elif reloctype==2:
        [nar,nndt,row_i,col_i,rw] = lsqrprep_st(log,ndt,nsrc,dt_idx,tmp_xp,tmp_yp,tmp_zp,
                                                dt_ista1,dt_ista2,dt_ic1,mod_ratio,wt)
    elif reloctype==3:
        [nar,nndt,row_i,col_i,rw] = lsqrprep_doub(log,ndt,nsrc,dt_idx,tmp_xp,tmp_yp,tmp_zp,
                                                  dt_ista1,dt_ista2,dt_ic1,dt_ic2,
                                                  mod_ratio,wt)

    """
    Scale G matrix
    ##############
    Such that L2 norm of each column is 1
    """
    log.write('~ Scaling G Columns \n')

    """
    Least square fitting 
    ##########
    Using the algorithm of Paige and Saunders 1982*
    * Same algorithm implemented in hypoDD and in scipy sparse lsqr
    """

    # Set up input parameter first
    m = nndt
    leniw = nar
    lenrw = nar

    if reloctype==3:
        n = 3*nev
    else:
        n = 4*nev

    w1 = np.zeros(n,dtype='float')
    w2 = np.zeros(n,dtype='float')
    x = np.zeros(n,dtype='float')
    se = np.zeros(n,dtype='float')

    # Declare variables
    atol = 0.000001
    btol = 0.000001
    conlim = 100000.0
    itnlim = 100*n
    istop = 0
    anorm = 0.0
    acond = 0.0
    rnorm = 0.0
    arnorm = 0.0
    xnorm = 0.0

    datet = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    log.write('~ lsqr .... %s \n\n' % datet)

    """
    Build A (or G) inversion matrix 
    ############## 
    Pulled from lsqr.f in hypoDD v1.3
    Using the sparse column matrix formatting in scipy
    """
    if reloctype==3:
        norm = np.zeros(3*nev,dtype='float')
    else:
        norm = np.zeros(4*nev,dtype='float')
    
    for i in range(nar):
        norm[col_i[i]] += rw[i]**2
    norm = np.sqrt(norm/nndt)

    if reloctype==3:
        for i in range(3*nev):
            if norm[col_i[i]]==0.:
                import pdb; pdb.set_trace()

    for i in range(nar):
        rw[i] = rw[i]/norm[col_i[i]]
    if reloctype==3:
        a = csc_matrix((rw,(row_i,col_i)),shape=(nndt,3*nev))
    else:
        a = csc_matrix((rw,(row_i,col_i)),shape=(nndt,4*nev))

    #if reloctype==3:
    #    for i in range(3*nev):
    #        if np.abs(norm[i])>1:
    #            raise Exception('(LSQR) Normtest: bad norm value')
    #else:
    #    for i in range(4*nev):
    #        if np.abs(norm[i])>1:
    #            raise Exception('(LSQR) Normtest: bad norm value')

    #import pdb; pdb.set_trace()
    """
    Run and unpack LSQR
    """
    [x,istop,itn,r1norm,r2norm,anorm,
     acond,arnorm,xnorm,var] = sp.lsqr(a,d,damp,atol,btol,conlim,itnlim,calc_var=True)
    log.write("  istop = %1i; acond (CND)= %8.1f; anorm = %8.1f; arnorm = %8.1f; xnorm = %8.1f; itn = %i\\%i \n" % 
              (istop, acond, anorm, arnorm, xnorm, itn, itnlim))
    normvar[iteri,0] = r1norm
    normvar[iteri,1] = r2norm
    normvar[iteri,2] = acond

    """
    Sensitivities
    """
    # if iteri==0:
    #     amean = np.abs(a)
    #     amean = amean.mean(0)
    #     amean = amean.A[0]
    #     if reloctype==3:
    #         xmean = np.mean(amean[0:3*nev-2:3])
    #         ymean = np.mean(amean[1:3*nev-1:3])
    #         zmean = np.mean(amean[2:3*nev:3])
    #     else:
    #         xmean = np.mean(amean[0:4*nev-3:4])
    #         ymean = np.mean(amean[1:4*nev-2:4])
    #         zmean = np.mean(amean[2:4*nev-1:4])
    #     if reloctype==1:
    #         xmean = xmean/8
    #         ymean = ymean/8
    #         zmean = zmean/8
    #     if reloctype==2:
    #         xmean = xmean/4
    #         ymean = ymean/4
    #         zmean = zmean/4
    #     if reloctype==3:
    #         xmean = xmean/6
    #         ymean = ymean/6
    #         zmean = zmean/6
    #     print('\nParameter Sensitivities Check')
    #     print('X: %f' % (xmean))
    #     print('Y: %f' % (ymean))
    #     print('Z: %f' % (zmean))
    #     log.write('\nParameter Sensitivities\n')
    #     log.write('X: %f \n' % (xmean))
    #     log.write('Y: %f \n' % (ymean))
    #     log.write('Z: %f \n' % (zmean))

    #     # import pdb; pdb.set_trace()
    #     # admean = amean*x
    #     # xmean = np.mean(amean[0:4*nev-3:4])
    #     # ymean = np.mean(amean[1:4*nev-2:4])
    #     # zmean = np.mean(amean[2:4*nev-1:4]) 
    #     dchange = 0.001   
    #     dnew = (dt_res[0:ndt]+dchange)*1000.*wt[0:ndt]
    #     [xch,istopch,itnch,r1normch,r2normch,anormch,
    #      acondch,arnormch,xnormch,varch] = sp.lsqr(a,dnew,damp,atol,btol,conlim,itnlim,calc_var=True)
    #     xchange = np.abs(x-xch)
    #     if reloctype==3:
    #         mx = np.mean(xchange[0:3*nev-2:3])
    #         my = np.mean(xchange[1:3*nev-1:3])
    #         mz = np.mean(xchange[2:3*nev:3])
    #     else:
    #         mx = np.mean(xchange[0:4*nev-3:4])
    #         my = np.mean(xchange[1:4*nev-2:4])
    #         mz = np.mean(xchange[2:4*nev-1:4])
    #     mx = dchange/mx
    #     my = dchange/my
    #     mz = dchange/mz
    #     print('\nParameter Sensitivities Check 2')
    #     print('X: %f' % (mx))
    #     print('Y: %f' % (my))
    #     print('Z: %f' % (mz))
    #     log.write('\nParameter Sensitivities\n')
    #     log.write('X: %f \n' % (mx))
    #     log.write('Y: %f \n' % (my))
    #     log.write('Z: %f \n' % (mz))   

    """
    Calculate standard error estimates
    ########
    Calculated for each model component
    """
    se = r2norm*np.sqrt(var/nndt)
    if nsrc==1:
        nsrc = nev
    # Rescale model vector
    x = x/norm
    se = se/norm

    """
    Unweight and rescale G matrix
    """
    if reloctype==1:
        rw[0:ndt] = rw[0:ndt]*wtinv[0:ndt]*norm[col_i[0:ndt]]
        rw[ndt:2*ndt] = rw[ndt:2*ndt]*wtinv[0:ndt]*norm[col_i[0:ndt]]
        rw[2*ndt:3*ndt] = rw[2*ndt:3*ndt]*wtinv[0:ndt]*norm[col_i[0:ndt]]
        rw[3*ndt:4*ndt] = rw[3*ndt:4*ndt]*wtinv[0:ndt]*norm[col_i[0:ndt]]
        rw[4*ndt:5*ndt] = rw[4*ndt:5*ndt]*wtinv[0:ndt]*norm[col_i[0:ndt]]
        rw[5*ndt:6*ndt] = rw[5*ndt:6*ndt]*wtinv[0:ndt]*norm[col_i[0:ndt]]
        rw[6*ndt:7*ndt] = rw[6*ndt:7*ndt]*wtinv[0:ndt]*norm[col_i[0:ndt]]
        rw[7*ndt:8*ndt] = rw[7*ndt:8*ndt]*wtinv[0:ndt]*norm[col_i[0:ndt]]
    elif reloctype==2:
        rw[0:ndt] = rw[0:ndt]*wtinv[0:ndt]*norm[col_i[0:ndt]]
        rw[ndt:2*ndt] = rw[ndt:2*ndt]*wtinv[0:ndt]*norm[col_i[0:ndt]]
        rw[2*ndt:3*ndt] = rw[2*ndt:3*ndt]*wtinv[0:ndt]*norm[col_i[0:ndt]]
        rw[3*ndt:4*ndt] = rw[3*ndt:4*ndt]*wtinv[0:ndt]*norm[col_i[0:ndt]]
    elif reloctype==3:
        rw[0:ndt] = rw[0:ndt]*wtinv[0:ndt]*norm[col_i[0:ndt]]
        rw[ndt:2*ndt] = rw[ndt:2*ndt]*wtinv[0:ndt]*norm[col_i[0:ndt]]
        rw[2*ndt:3*ndt] = rw[2*ndt:3*ndt]*wtinv[0:ndt]*norm[col_i[0:ndt]]
        rw[3*ndt:4*ndt] = rw[3*ndt:4*ndt]*wtinv[0:ndt]*norm[col_i[0:ndt]]
        rw[4*ndt:5*ndt] = rw[4*ndt:5*ndt]*wtinv[0:ndt]*norm[col_i[0:ndt]]
        rw[5*ndt:6*ndt] = rw[5*ndt:6*ndt]*wtinv[0:ndt]*norm[col_i[0:ndt]]
    """
    Compute residuals from d = G*x
    """
    dt_resold = dt_res
    d = -1*np.copy(dt_res)*1000.0
    for k in range(nar):
        d[row_i[k]] = d[row_i[k]] + rw[k]*x[col_i[k]]
    dt_res = -1*d/1000.0

    # Sam value
    sam = np.zeros(3,dtype='float')
    sam[0] = r1norm
    sam[1] = xnorm
    sam[2] = r2norm

    """
    # Compute Covariance and Correlation Matrix ---- KB ADDED FOR TESTING
    """
    ##if reloctype==3:
    #    a = csc_matrix((np.abs(rw),(row_i,col_i)),shape=(nndt,3*nev))
    #else:
    #    a = csc_matrix((np.abs(rw),(row_i,col_i)),shape=(nndt,4*nev))
    #sat = a.transpose()

    # corr = np.zeros(covar.shape,dtype='float')
    # if reloctype==3:
    #     for i in range(3*nev):
    #         for j in range(3*nev):
    #             corr[i,j] = covar[i,j]/np.sqrt(covar[i,i]*covar[j,j])
    # else:
    #     for i in range(4*nev):
    #         for j in range(4*nev):
    #             corr[i,j] = covar[i,j]/np.sqrt(covar[i,i]*covar[j,j])

    # if iteri==1 or iteri==10:
    #     plt.pcolormesh(covar,cmap='seismic',vmin=-1.,vmax=1.)
    #     plt.suptitle('Model Covariance')
    #     plt.xlabel('Model Parameter Index')
    #     plt.ylabel('Model Parameter Index')
    #     plt.colorbar()
    #     plt.gca().invert_yaxis()
    #     plt.savefig('figures/model_covar_%i.png' % reloctype)
    #     plt.close('all')

    #     plt.pcolormesh(corr,cmap='seismic',vmin=-1.,vmax=1.)
    #     plt.suptitle('Model Correlation')
    #     plt.xlabel('Model Parameter Index')
    #     plt.ylabel('Model Parameter Index')
    #     plt.colorbar()
    #     plt.gca().invert_yaxis()
    #     plt.savefig('figures/model_corr_%i.png' % reloctype)
    #     plt.close('all')

    #     print('\nCorrelations Inversion Iteration %i' % iteri)
    #     print('For event 1')
    #     print('\nX to Y Correlation: ', corr[0,1])
    #     print('X to Z Correlation: ', corr[0,2])
    #     print('Y to Z Correlation: ',corr[1,2])

    """
    Get residual statistics (avrg, rms, var..)
    """
    resvar1 = 0.
    [rms_cc,rms_ct,rms_cc0,rms_ct0, rms_ccold,rms_ctold,
     rms_cc0old,rms_ct0old,resvar1] = resstat(log,reloctype,idata,ndt,nev,dt_res,dt_wt,dt_idx,
                                              rms_cc,rms_ct,rms_cc0,rms_ct0,
                                              rms_ccold,rms_ctold,rms_cc0old,rms_ct0old,
                                              resvar1)

    """
    Scale errors
    ############
    The standard error estimates returned by LSQR increase with the no. of iterations.  
    If LSQR shuts down early because of loose tolerances,
    or because the rhs-vector is special, the estimates will be too small.

    Remember that se(j) is covariance(j) / (m - n)
    where m - n = 1000000.  I've never quite understood why we
    divide by that number.
    """

    # Errors for the 95% confidence level,
    # thus multiply the standard errors by 2.7955
    factor = 2.7955

    """
    Unpack lsqr results
    """
    if reloctype==1:
        [src_dx,src_dy,src_dz,src_dt,
         src_ex,src_ey,src_ez,src_et,src_cusp,
         exavold,eyavold,ezavold,etavold,
         dxavold,dyavold,dzavold,dtavold,
         exav,eyav,ezav,etav,
         dxav,dyav,dzav,dtav] = lsqrunpack_ev(log,nev,iteri,x,se,resvar1,factor,ev_cusp,
                                              exav,eyav,ezav,etav,dxav,dyav,dzav,dtav)
        src_ds = np.zeros(nev,dtype='float')
        src_es = np.zeros(nev,dtype='float')
        esav = float(0.)
        dsav = float(0.)

    elif reloctype==2:
        [src_dx,src_dy,src_dz,src_ds,
         src_ex,src_ey,src_ez,src_es,src_cusp,
         exavold,eyavold,ezavold,esavold,
         dxavold,dyavold,dzavold,dsavold,
         exav,eyav,ezav,esav,
         dxav,dyav,dzav,dsav] = lsqrunpack_st(log,nev,iteri,x,se,resvar1,factor,ev_cusp,
                                              exav,eyav,ezav,esav,dxav,dyav,dzav,dsav)
        src_dt = np.zeros(nev,dtype='float')
        src_et = np.zeros(nev,dtype='float')
        etav = float(0.)
        dtav = float(0.)

    elif reloctype==3:
        [src_dx,src_dy,src_dz,
         src_ex,src_ey,src_ez,src_cusp,
         exavold,eyavold,ezavold,
         dxavold,dyavold,dzavold,
         exav,eyav,ezav,
         dxav,dyav,dzav] = lsqrunpack_doub(log,nev,iteri,x,se,resvar1,factor,ev_cusp,
                                           exav,eyav,ezav,dxav,dyav,dzav)
        src_dt = np.zeros(nev,dtype='float')
        src_et = np.zeros(nev,dtype='float')
        etav = float(0.)
        dtav = float(0.)
        src_ds = np.zeros(nev,dtype='float')
        src_es = np.zeros(nev,dtype='float')
        esav = float(0.)
        dsav = float(0.)


    # top = np.mean(np.abs((dt_resold[:]-dt_res[:])/dt_res[:]))
    # print('\nParameter Sensitivities Check 1')
    # print('X: %f' % (1000*top/dxav))
    # print('Y: %f' % (1000*top/dyav))
    # print('Z: %f' % (1000*top/dzav))
    # log.write('\nParameter Sensitivities\n')
    # log.write('X: %f \n' % (1000*top/dxav))
    # log.write('Y: %f \n' % (1000*top/dyav))
    # log.write('Z: %f \n' % (1000*top/dzav))
    # print('\nParameter Sensitivities Check 2')
    # print('X: %f' % (1000/acond*top/dxav))
    # print('Y: %f' % (1000/acond*top/dyav))
    # print('Z: %f' % (1000/acond*top/dzav))
    # log.write('\nParameter Sensitivities\n')
    # log.write('X: %f \n' % (1000/acond*top/dxav))
    # log.write('Y: %f \n' % (1000/acond*top/dyav))
    # log.write('Z: %f \n' % (1000/acond*top/dzav))

    # Return updated location arrays and inversion statistics
    return [src_cusp,src_dx,src_dy,src_dz,src_dt,src_ds, 
            src_ex,src_ey,src_ez,src_et,src_es,
            exav,eyav,ezav,etav,esav,dxav,dyav,dzav,dtav,dsav,
            rms_cc,rms_ct,rms_cc0,rms_ct0,
            rms_ccold,rms_ctold,rms_cc0old,rms_ct0old,acond,normvar,sam]




