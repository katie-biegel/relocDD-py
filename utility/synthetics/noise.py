# Import General Python Packages
import numpy as np

"""
This script includes the functions necessary to add synthetic
noise on a model.
---
This script contains the noise subfunctions:
    noisetofile:        Generate noise and write directly to file,
    noisearrays:        Generate noise for array passing;
And the main noise subfunction:
    generate_noise:     Generate noise for synthetic problems.
---
"""

def noisetofile(reloctype,noisediff,stdct,stdcc,fn_cc,fn_ct,phdat='phase.dat',nev=15,nsta=64):
    """
    Generate Noise and Add to Files
    """

    """
    Add single noise measurement per data (theory noise replication)
    """
    if noisediff==0:
        """
        ADD CC NOISE FIRST
        """

        """
        Open cc file
        """
        ccfile = open(fn_cc,'r')
        ccs = ccfile.readlines()
        ccfile.close()

        ccfile = open(fn_cc,'w')

        """
        Generate CC Noise Array
        """
        noise_cc = np.random.normal(0.,stdcc,(len(ccs)))

        """
        Add Noise Array to file
        """
        ncc=int(0)
        for cc in ccs:
            if cc[0]=='#':
                ccfile.write(cc)
            else:
                cc=cc.split(' ')
                cc=list(filter(None,cc))

                if reloctype==1:
                    ccfile.write('%s %f %s %s \n' % 
                                 (cc[0],float(cc[1])+noise_cc[ncc],cc[2],cc[3]))
                if reloctype==2:
                    ccfile.write('%s %s %f %s %s \n' % 
                                 (cc[0],cc[1],float(cc[2])+noise_cc[ncc],cc[3],cc[4]))
                if reloctype==3:
                    ccfile.write('%s %s %f %s %s %s \n' % 
                                 (cc[0],cc[1],float(cc[2])+noise_cc[ncc],cc[3],cc[4],cc[5]))
                ncc += 1
        ccfile.close()

        """
        ADD CT NOISE
        """

        """
        Open ct file
        """
        ctfile = open(fn_ct,'r')
        cts = ctfile.readlines()

        ctfile.close()
        ctfile = open(fn_ct,'w')

        """
        Generate CT Noise Array
        """
        noise_ct = np.random.normal(0.,stdct,(len(cts)))

        """
        Add Noise Array to file
        """
        nct=int(0)
        for ct in cts:
            if ct[0]=='#':
                ctfile.write(ct)
            else:
                ct=ct.split(' ')
                ct=list(filter(None,ct))

                """
                For the catalog data there are 2 tt measurements, but since noisediff==0
                we only want to add one noise, so noise is only added to 1 tt measurement
                """
                if reloctype==1:
                    ctfile.write('%s %f %s %s %s \n' % 
                                 (ct[0],float(ct[1])+noise_ct[nct],ct[2],ct[3],ct[4]))
                if reloctype==2:
                    ctfile.write('%s %s %f %s %s %s \n' % 
                                 (ct[0],ct[1],float(ct[2])+noise_ct[nct],ct[3],ct[4],ct[5]))
                if reloctype==3:
                    ctfile.write('%s %s %f %s %s %s %s %s \n' % 
                                 (ct[0],ct[1],float(ct[2])+noise_ct[nct],ct[3],ct[4],ct[5],ct[6],ct[7]))
                nct += 1
        ctfile.close()

        """
        Save noise to file
        """
        x = np.concatenate((noise_cc[0:ncc],noise_ct[0:nct]))
        np.savetxt('noise.txt',x)

    else:
        """
        Noise difference value
        ---
        Equivalent to noise added in ph2dt (measurement error)
        ---
        Noise added to phase measurements
        Noise then compounded (noise difference taken) for dt_dt measurements
        (For double-pair (noise x4 since four phase measurements per dt_dt measurement))
        Need to generate 1 noise value for every ray for catalog
        Same number for ccs
        ---
        """
        nscts = np.random.normal(0.,stdct,(nev,2*nsta))
        nsccs = np.random.normal(0.,stdcc,(nev,2*nsta))  

        """
        Add noise to phase.dat
        """
        phfile = open(phdat,'r')
        phases = phfile.readlines()
        phfile.close()

        iev = int(-1)
        ista = int(0)
        stas = np.empty((nsta),dtype='object')
        phfile = open(phdat,'w')
        for ph in phases:
            if ph[0]=='#':
                phfile.write(ph)
                ista=0
                iev+=1
            else:
                if iev<nev:
                    if iev==0:
                        stas[ista] = ph[0:2]
                    ph = ph.split(' ')
                    ph = list(filter(None,ph))
                    phfile.write('%s %f %s %s \n' % (ph[0],float(ph[1])+nscts[iev,ista],ph[2],ph[3]))
                    ista+=1
                elif iev>=nev:
                    ph = ph.split(' ')
                    ph = list(filter(None,ph))
                    phfile.write('%s %f %s %s \n' % (ph[0],float(ph[1])+nscts[iev-nev,nsta+ista],ph[2],ph[3]))
                    ista+=1
        phfile.close()

        """
        ADD CC NOISE FIRST
        """

        """
        Open cc file
        """
        ccfile = open(fn_cc,'r')
        ccs = ccfile.readlines()
        ccfile.close()

        ccfile = open(fn_cc,'w')

        """
        Add Noise Array to file
        """
        noise_cc = []
        ncc=int(0)
        for cc in ccs:
            if cc[0]=='#':
                ccfile.write(cc)

                cc = list(filter(None,cc.split(' ')))
                ic1 = int(cc[1])-1
                if reloctype!=2:
                    ic2 = int(cc[2])-1
            else:
                cc=cc.split(' ')
                cc=list(filter(None,cc))

                if reloctype==1:
                    sta1 = cc[0]
                    ista1 = np.where(stas==sta1)[0][0]

                    if cc[3]=='P':
                        ccfile.write('%s %f %s %s \n' % 
                                     (cc[0],float(cc[1])+nsccs[ic1,ista1]-nsccs[ic2,ista1],
                                      cc[2],cc[3]))
                        noise_cc.append(nsccs[ic1,ista1]-nsccs[ic2,ista1])
                    if cc[3]=='S':
                        ccfile.write('%s %f %s %s \n' % 
                                     (cc[0],float(cc[1])+nsccs[ic1,nsta+ista1]-nsccs[ic2,nsta+ista1],
                                      cc[2],cc[3]))
                        noise_cc.append(nsccs[ic1,nsta+ista1]-nsccs[ic2,nsta+ista1])

                if reloctype==2:
                    sta1 = cc[0]
                    sta2 = cc[1]
                    ista1 = np.where(stas==sta1)[0][0]
                    ista2 = np.where(stas==sta2)[0][0]

                    if cc[4]=='P':
                        ccfile.write('%s %s %f %s %s \n' % 
                                     (cc[0],cc[1],
                                      float(cc[2])+nsccs[ic1,ista1]-nsccs[ic1,ista2],
                                      cc[3],cc[4]))
                        noise_cc.append(nsccs[ic1,ista1]-nsccs[ic1,ista2])

                    if cc[4]=='S':
                        ccfile.write('%s %s %f %s %s \n' % 
                                     (cc[0],cc[1],
                                      float(cc[2])+nsccs[ic1,ista1+nsta]-nsccs[ic1,ista2+nsta],
                                      cc[3],cc[4]))
                        noise_cc.append(nsccs[ic1,ista1+nsta]-nsccs[ic1,ista2+nsta])

                if reloctype==3:
                    sta1 = cc[0]
                    sta2 = cc[1]
                    ista1 = np.where(stas==sta1)[0][0]
                    ista2 = np.where(stas==sta2)[0][0]

                    if cc[5]=='P':
                        ccfile.write('%s %s %f %f %s %s \n' % 
                                     (cc[0],cc[1],
                                      float(cc[2])+nsccs[ic1,ista1]-nsccs[ic2,ista1],
                                      float(cc[3])+nsccs[ic1,ista2]-nsccs[ic2,ista2],
                                      cc[4],cc[5]))
                        noise_cc.append((nsccs[ic1,ista1]-nsccs[ic2,ista1]) - 
                                        (nsccs[ic1,ista2]-nsccs[ic2,ista2]))

                    if cc[5]=='S':
                        ccfile.write('%s %s %f %f %s %s \n' % 
                                     (cc[0],cc[1],
                                      float(cc[2])+nsccs[ic1,nsta+ista1]-nsccs[ic2,nsta+ista1],
                                      float(cc[3])+nsccs[ic1,nsta+ista2]-nsccs[ic2,nsta+ista2],
                                      cc[4],cc[5]))
                        noise_cc.append((nsccs[ic1,nsta+ista1]-nsccs[ic2,nsta+ista1]) - 
                                        (nsccs[ic1,nsta+ista2]-nsccs[ic2,nsta+ista2]))

                ncc += 1
        ccfile.close()
        noise_cc = np.asarray(noise_cc)

        """
        ADD CT NOISE
        """

        """
        Open ct file
        """
        ctfile = open(fn_ct,'r')
        cts = ctfile.readlines()

        ctfile.close()
        ctfile = open(fn_ct,'w')

        """
        Add Noise Array to file
        """
        nct=int(0)
        noise_ct = []
        for ct in cts:
            if ct[0]=='#':
                ctfile.write(ct)

                ct = list(filter(None,ct.split(' ')))
                ic1 = int(ct[1])-1
                if reloctype!=2:
                    ic2 = int(ct[2])-1
            else:
                ct=ct.split(' ')
                ct=list(filter(None,ct))

                """
                For the catalog data there are 2 tt measurements, but since noisediff==0
                we only want to add one noise, so noise is only added to 1 tt measurement
                """
                if reloctype==1:
                    sta1 = ct[0]
                    ista1 = np.where(stas==sta1)[0][0]

                    if ct[4]=='P':
                        ctfile.write('%s %f %f %s %s \n' % 
                                     (ct[0],
                                      float(ct[1])+nscts[ic1,ista1],
                                      float(ct[2])+nscts[ic2,ista1],
                                      ct[3],ct[4]))
                        noise_ct.append(nscts[ic1,ista1]-nscts[ic2,ista1])

                    if ct[4]=='S':
                        ctfile.write('%s %f %f %s %s \n' % 
                                     (ct[0],
                                      float(ct[1])+nscts[ic1,nsta+ista1],
                                      float(ct[2])+nscts[ic2,nsta+ista1],
                                      ct[3],ct[4]))
                        noise_ct.append(nscts[ic1,nsta+ista1]-nscts[ic2,nsta+ista1])

                if reloctype==2:
                    sta1 = ct[0]
                    sta2 = ct[1]
                    ista1 = np.where(stas==sta1)[0][0]
                    ista2 = np.where(stas==sta2)[0][0]

                    if ct[5]=='P':
                        ctfile.write('%s %s %f %f %s %s \n' % 
                                     (ct[0],ct[1],
                                      float(ct[2])+nscts[ic1,ista1],
                                      float(ct[3])+nscts[ic1,ista2],
                                      ct[4],ct[5]))
                        noise_ct.append(nscts[ic1,ista1]-nscts[ic1,ista2])
                    if ct[5]=='S':
                        ctfile.write('%s %s %f %f %s %s \n' % 
                                     (ct[0],ct[1],
                                      float(ct[2])+nscts[ic1,nsta+ista1],
                                      float(ct[3])+nscts[ic1,nsta+ista2],
                                      ct[4],ct[5]))
                        noise_ct.append(nscts[ic1,nsta+ista1]-nscts[ic1,nsta+ista2])

                if reloctype==3:
                    sta1 = ct[0]
                    sta2 = ct[1]
                    ista1 = np.where(stas==sta1)[0][0]
                    ista2 = np.where(stas==sta2)[0][0]

                    if ct[7]=='P':
                        ctfile.write('%s %s %f %f %f %f %s %s \n' % 
                                     (ct[0],ct[1],
                                      float(ct[2])+nscts[ic1,ista1],
                                      float(ct[3])+nscts[ic1,ista2],
                                      float(ct[4])+nscts[ic2,ista1],
                                      float(ct[5])+nscts[ic2,ista2],
                                      ct[6],ct[7]))
                        noise_ct.append((nscts[ic1,ista1]-nscts[ic2,ista1]) - 
                                        (nscts[ic1,ista2]-nscts[ic2,ista2]))
                    if ct[7]=='S':
                        ctfile.write('%s %s %f %f %f %f %s %s \n' % 
                                     (ct[0],ct[1],
                                      float(ct[2])+nscts[ic1,nsta+ista1],
                                      float(ct[3])+nscts[ic1,nsta+ista2],
                                      float(ct[4])+nscts[ic2,nsta+ista1],
                                      float(ct[5])+nscts[ic2,nsta+ista2],
                                      ct[6],ct[7]))
                        noise_ct.append((nscts[ic1,nsta+ista1]-nscts[ic2,nsta+ista1]) - 
                                        (nscts[ic1,nsta+ista2]-nscts[ic2,nsta+ista2]))
                nct += 1
        ctfile.close()
        noise_ct = np.asarray(noise_ct)

        """
        Save noise to file
        """
        x = np.concatenate((noise_cc,noise_ct))
        np.savetxt('noise.txt',x)

    return None


def noisearrays(reloctype,noisediff,nsta,nev,nct,ncc,stdct,stdcc,p_time,
                dt_dt,dt_ic1,dt_ic2,dt_ista1,dt_ista2,dt_idx):
    """
    Generate Noise and Add to Data Arrays
    """

    """
    Define Noise Arrays
    """
    if noisediff==0:   
        """ 
        Single noise value
        ---
        Equivalent to noise added in hypoDD
        One noise measurment per data
        Equivalent to theory noise in relocation
        """
        noise_dt = np.random.normal(0.,stdct,(nct))
        noise_cc = np.random.normal(0.,stdcc,(ncc))

        """
        Create general Noise Array
        """
        x = np.concatenate((noise_cc,noise_dt))     # Save noise array

        """
        Update dt_dt
        """
        dt_dt += x

    else:
        """
        Noise difference value
        ---
        Equivalent to noise added in ph2dt (measurement error)
        ---
        Noise added to phase measurements
        Noise then compounded (noise difference taken) for dt_dt measurements
        (For double-pair (noise x4 since four phase measurements per dt_dt measurement))
        Need to generate 1 noise value for every ray for catalog
        Same number for ccs
        ---
        """
        nscts = np.random.normal(0.,stdct,(nev,2*nsta))
        nsccs = np.random.normal(0.,stdcc,(nev,2*nsta))           

        """
        Update p_time 
        ---
        Only if noise_diff is 1
        p_time shape is (nev,2*nsta) 
            - P saved (:,:nsta)
            - S saved (:,nsta:2*nsta)
        """
        p_time += nscts

        """
        Update dt_dt and x arrays
        """
        ndt = len(dt_dt)
        x = np.zeros(ndt,dtype='float')
        #import pdb; pdb.set_trace()
        for i in range(ndt):
            if dt_idx[i]==1:
                if reloctype==1:
                    x[i] = nsccs[dt_ic1[i],dt_ista1[i]] - nsccs[dt_ic2[i],dt_ista1[i]]
                elif reloctype==2:
                    x[i] = nsccs[dt_ic1[i],dt_ista1[i]] - nsccs[dt_ic1[i],dt_ista2[i]]
                elif reloctype==3:
                    x[i] = ((nsccs[dt_ic1[i],dt_ista1[i]] - nsccs[dt_ic2[i],dt_ista1[i]])-
                            (nsccs[dt_ic1[i],dt_ista2[i]] - nsccs[dt_ic2[i],dt_ista2[i]]))
            elif dt_idx[i]==2:
                if reloctype==1:
                    x[i] = nsccs[dt_ic1[i],nsta+dt_ista1[i]] - nsccs[dt_ic2[i],nsta+dt_ista1[i]]
                elif reloctype==2:
                    x[i] = nsccs[dt_ic1[i],nsta+dt_ista1[i]] - nsccs[dt_ic1[i],nsta+dt_ista2[i]]
                elif reloctype==3:
                    x[i] = ((nsccs[dt_ic1[i],nsta+dt_ista1[i]] - nsccs[dt_ic2[i],nsta+dt_ista1[i]])-
                            (nsccs[dt_ic1[i],nsta+dt_ista2[i]] - nsccs[dt_ic2[i],nsta+dt_ista2[i]]))
            elif dt_idx[i]==3:
                if reloctype==1:
                    x[i] = nscts[dt_ic1[i],dt_ista1[i]] - nscts[dt_ic2[i],dt_ista1[i]]
                elif reloctype==2:
                    x[i] = nscts[dt_ic1[i],dt_ista1[i]] - nscts[dt_ic1[i],dt_ista2[i]]
                elif reloctype==3:
                    x[i] = ((nscts[dt_ic1[i],dt_ista1[i]] - nscts[dt_ic2[i],dt_ista1[i]])-
                            (nscts[dt_ic1[i],dt_ista2[i]] - nscts[dt_ic2[i],dt_ista2[i]]))
            elif dt_idx[i]==4:
                if reloctype==1:
                    x[i] = nscts[dt_ic1[i],nsta+dt_ista1[i]] - nscts[dt_ic2[i],nsta+dt_ista1[i]]
                elif reloctype==2:
                    x[i] = nscts[dt_ic1[i],nsta+dt_ista1[i]] - nscts[dt_ic1[i],nsta+dt_ista2[i]]
                elif reloctype==3:
                    x[i] = ((nscts[dt_ic1[i],nsta+dt_ista1[i]] - nscts[dt_ic2[i],nsta+dt_ista1[i]])-
                            (nscts[dt_ic1[i],nsta+dt_ista2[i]] - nscts[dt_ic2[i],nsta+dt_ista2[i]]))

        """
        Uspdate dt_dt arrays
        """
        dt_dt += x

    """
    Save Noise to File
    """
    np.savetxt('noise.txt',x)

    #import pdb; pdb.set_trace()

    return x,dt_dt,p_time


def generate_noise(fileout=0,reloctype=1,noisediff=0,stdcc=0.,stdct=0.,
                   data=['dt.cc','dt.ct']):
    """
    Generate Noise for Synthetic Models
    :::
    Parameters:
    fileout (int) --- File output switch
                      0 for files on
                      1 for array passing
    reloctype (int) --- Double-difference pairing switch
    noisediff (int) --- Noise type switch
                        0 for single noise value per data (i.e. added in hypoDD)
                        1 for noise diff valuer per data (i.e. added to tt in ph2dt)
    stdcc (float) --- Std. of noise for cc measurements
    stdct (float) --- Std. of noise for ct measurements
    data (list) --- List of input data.  Dependent on fileout.
                    if fileout==0:
                        data = [fn_cc,fn_ct]
                    if fileout==1:
                        data = [nsta,nev,ncc,nct,p_time,dt_dt,dt_ic1,dt_ic2,
                                dt_ista1,dt_ista2,dt_idx]
    :::
    Returns:
    retlist (list) --- Either empty (None) if fileout==0 or 
                       list of [x,dt_dt,p_time]if fileout==1
    :::
    """

    if fileout==0:
        """
        Unpack file names
        """
        [fn_cc,fn_ct] = data

        """
        Add noise to files
        """
        noisetofile(reloctype,noisediff,stdct,stdcc,fn_cc,fn_ct)

        """
        Noise added to file
        ---
        Return Nothing
        """
        return None

    elif fileout==1:
        """
        Unpack data arrays
        """
        [nsta,nev,ncc,nct,p_time,dt_dt,dt_ic1,dt_ic2,dt_ista1,dt_ista2,dt_idx] = data

        """
        Add noise to passable arrays
        """

        noise_return = noisearrays(reloctype,noisediff,nsta,nev,nct,ncc,stdct,stdcc,p_time,
                                   dt_dt,dt_ic1,dt_ic2,dt_ista1,dt_ista2,dt_idx)

        """
        Return Noise
        """
        return noise_return

    """
    Raise exception if it makes it to here
    """
    raise Exception('Noise Error. End of Function generate_noise.')