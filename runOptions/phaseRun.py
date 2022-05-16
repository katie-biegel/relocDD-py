# Other hypoDD-py functions
from utility.synthetics.synthModel import synth_generate
from utility.synthetics.noise import generate_noise
from utility.synthetics.hypoinv import hypoinverse_wrapper

from ph2dt.ph2dt import ph2dt
from hypoDD.hypoDD import hypoDD


def phaseRun(log,pinput='ph2dt.inp',hinput='hypoDD.inp',reloctype=1):
    """
    Run classic syntheticRun.
    No noise or hypoinverse step. Real data files used. 
    All input/outputs files. 
    """
    log.write('\n\nStarting from phase.dat and station.dat.\n\n')
    print('Starting from phase.dat and station.dat')

    """
    Run ph2dt
    """
    log.write('Run ph2dt.\n\n')
    print('Run ph2dt.')
    hypoDDinput = ph2dt(log,pinput,reloctype=reloctype,fileout=0,returnval=True)
    print('Ph2dt completed.')

    """
    Run hypoDD
    """
    log.write('Run relocation.\n\n')
    print('Run relocation.')

    outputs = hypoDD(log,hinput,reloctype,fileout=1,
                     hypoDDdata=hypoDDinput)
    print('Relocation completed.')

    return outputs