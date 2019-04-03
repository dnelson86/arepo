""" @package examples/cosmo_box_star_formation_3d/create.py
Code creates the output list; ICs need to be present already

created by Rainer Weinberger, last modified 28.02.2019
"""


""" load libraries """
import sys    # system calls
import numpy as np    # scientific computing package
import h5py    # hdf5 format
import os    # operating system interface
from subprocess import call   # execute bash commands


## create new ics with 'music' or 'ngenic' or just 'copy' existing 
## ones to the run directory
ic_creation='copy'

""" input """
simulation_directory = str(sys.argv[1])
print("create.py " + simulation_directory)


""" initial conditions: either copy or create with code """
if ic_creation == 'copy':
    call(['cp', simulation_directory+'/L8n32/ics', simulation_directory+'/ics'])
elif ic_creation == 'music':
    status = call(['hg', 'clone', 'https://bitbucket.org/ohahn/music', simulation_directory+'/music'])
    if status != 0:
        print('CREATE: ERROR: hg clone failed!')
        sys.exit(status)
    cwd = os.getcwd()
    os.chdir(simulation_directory+'/music/')
    status = call(['make'])
    if status != 0:
        print('CREATE: ERROR: make failed!')
        sys.exit(status)
    status = call(['./MUSIC',cwd+'/param_music.txt'])
    if status != 0:
        print('CREATE: ERROR: execution failed!')
        sys.exit(status)
    os.chdir(cwd)
elif ic_creation == 'ngenic':
    status = call(['git', 'clone', 'https://gitlab.mpcdf.mpg.de/ext-c2c74fbfcdff/ngenic.git', simulation_directory+'/ngenic'])
    if status != 0:
        print('CREATE: ERROR: git clone failed!')
        sys.exit(status)
    cwd = os.getcwd()
    os.chdir(simulation_directory+'/ngenic/')
    status = call(['make'])
    if status != 0:
        print('CREATE: ERROR: make failed!')
        sys.exit(status)
    status = call(['mpiexec','-np','1','./N-GenIC',cwd+'/param_ngenic.txt'])
    if status != 0:
        print('CREATE: ERROR: execution failed!')
        sys.exit(status)
    os.chdir(cwd)
else:
    print("CREATE: ERROR: no valid option for ic creation! choose 'copy', 'music' or 'ngenic'")
    exit(1)


""" set output times """
outputTimes = np.array([0.2,0.25,0.33,0.5,0.66,1], dtype=np.float64)
ones = np.ones(outputTimes.shape, dtype=np.int)

""" copy treecool file to run directory """
if len(sys.argv) > 2:
    arepopath = sys.argv[2]
else:
    arepopath = '.'
call(['cp', arepopath + '/data/TREECOOL_ep', simulation_directory+'/TREECOOL_ep'])

""" write output list file """
data = np.array([outputTimes, ones]).T
np.savetxt(simulation_directory+"/output_list.txt",data, fmt="%g %1.f" )


""" normal exit """
sys.exit(0) 
