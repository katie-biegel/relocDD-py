# Other hypoDD-py functions
from utility.synthetics.synthModel import synth_generate
from utility.synthetics.noise import generate_noise
from utility.synthetics.hypoinv import hypoinverse_wrapper

from ph2dt.ph2dt import ph2dt
from hypoDD.hypoDD import hypoDD


def limitedRun(log,pinput,hinput,reloctype,makedata,
               noiseswitch,noisediff,hypoinv,stdcc,stdct):
    """
    To run hypoDD with passing arrays and no file outputs
    """
    log.write('\n\nStarting limited run hypoDD.\n\n')
    print('Starting limited run hypoDD.')

    """
    Build synthethic model if needed Synthetic Model
    ---
    With fileout turned on, will generate data files needed for ph2dt.
    """
    if makedata==1:
        log.write('Build synthetic model.\n\n')
        print('Build synthetic model.')
        [nsta,s_lab,s_lat,s_lon,dist,az,ang,
         nev,lat,lon,depth,cuspid,dates,times,mag,herr,zerr,res,
         npha,nobs_ct,p_pha,p_sta,p_time,p_wghtr] = synth_generate(log,fileout=1)
        print('Phase.dat file generated.')

    """
    Run ph2dt
    """
    log.write('Run ph2dt.\n\n')
    print('Run ph2dt.')
    if makedata==1:
        ph2dtinput = [nsta,s_lab,s_lat,s_lon,
                      nev,lat,lon,depth,cuspid,dates,times,mag,herr,zerr,res,
                      npha,nobs_ct,p_pha,p_sta,p_time,p_wghtr]
        [dt_sta1,dt_sta2,dt_dt,dt_qual,dt_c1,dt_c2,
         dt_idx,dt_ista1,dt_ista2,dt_ic1,dt_ic2,dt_offse,dt_offss,
         ndt,nccp,nccs,nctp,ncts] = ph2dt(log,pinput,ph2dtinput,reloctype=reloctype,
                                          fileout=1,makedata=makedata)
    else:
        [dt_sta1,dt_sta2,dt_dt,dt_qual,dt_c1,dt_c2,
         dt_idx,dt_ista1,dt_ista2,dt_ic1,dt_ic2,dt_offse,dt_offss,
         ndt,nccp,nccs,nctp,ncts] = ph2dt(log,pinput,reloctype,fileout=1)
    print('Ph2dt completed.')
    ncc = nccp+nccs 
    nct = nctp+ncts

    if makedata==1:

        if noiseswitch==1:
            """
            Add noise to datafiles
            """
            log.write('Add Noise.\n\n')
            print('Add Noise.')
            [x,dt_dt,p_time] = generate_noise(1,reloctype,noisediff,stdcc,stdct,
                                              data=[nsta,nev,ncc,nct,p_time,dt_dt,
                                                    dt_ic1,dt_ic2,dt_ista1,dt_ista2,dt_idx])
            print('Add Noise completed.')
        if noiseswitch==0:
            """
            Add noise to datafiles
            """
            log.write('No Noise.\n\n')
            print('No Noise.')
            x = []

        if hypoinv==1:
            """
            Run hypoinverse catalog step if called
            """
            log.write('Updated Event Locations using HypoINVERSE-2000.\n\n')
            print('Updated Event Locations using HypoINVERSE-2000.')
            [dates,times,lat,lon,depth,mag,
             res,herr,verr,cuspid] = hypoinverse_wrapper(log,1,hinput,'synth.inp',
                                                         data=[nsta,s_lab,s_lat,s_lon,
                                                               nev,lat,lon,depth,cuspid,dates,
                                                               times,mag,herr,zerr,res,
                                                               p_time,dist,az,ang])

            print('Catalog Locations completed.')

    """
    Run hypoDD
    """
    log.write('Run relocation.\n\n')
    print('Run relocation.')

    [calstart,calend,dtdtstart,dtdtend,
     locdel,locabs,loctru] = hypoDD(log,hinput,reloctype,fileout=1,
                                    hypoDDdata=[dates,times,cuspid,lat,lon,depth,mag,herr,zerr,res,
                                                s_lab,s_lat,s_lon,
                                                dt_sta1,dt_sta2,dt_dt,dt_qual,dt_c1,dt_c2,dt_idx,
                                                dt_ista1,dt_ista2,dt_ic1,dt_ic2,dt_offse,dt_offss,
                                                nev,nsta,ndt,nccp,nccs,nctp,ncts])

    print('Relocation completed.')

    if makedata==1:
        return [calstart,calend,dtdtstart,dtdtend,locdel,locabs,loctru,x,ncc] 
    else:
        return [calstart,calend,dtdtstart,dtdtend,locdel,locabs,loctru]