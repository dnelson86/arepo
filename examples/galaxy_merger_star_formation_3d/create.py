""" @package examples/galaxy_merger_star_formation_3d/create.py
Code creates the output list; ICs need to be present already
galaxy_merger_star_formation_3d uses MakeNewDisk and CombineGalaxies codes
to create the initial conditions

created by Rainer Weinberger, last modified 09.03.2019
"""


""" load libraries """
import sys    # system calls
import os  # operating system calls
import numpy as np    # scientific computing package
import h5py    # hdf5 format
from subprocess import call    # execute shell commands 


""" input """
simulation_directory = str(sys.argv[1])
print("create.py " + simulation_directory)


""" set output times """
outputTimes = np.linspace(0.0,3.0,32, dtype=np.float64)
ones = np.ones(outputTimes.shape, dtype=np.int)


""" write output list file """
data = np.array([outputTimes, ones]).T
np.savetxt(simulation_directory+"/output_list.txt",data, fmt="%g %1.f" )


""" copy treecool file to run directory """
if len(sys.argv) > 2:
    arepopath = sys.argv[2]
else:
    arepopath ='.'
call(['cp', arepopath+'/data/TREECOOL_ep', simulation_directory+'/TREECOOL_ep'])

cwd = os.getcwd()

""" create backgroundgrid ICs from SPH ICs """
## compile Arepo with ADDBACKGROUNDGRID
os.chdir(arepopath)
res = call(["make", "CONFIG="+simulation_directory+"/Config_ADDBACKGROUNDGRID.sh", \
             "BUILD_DIR="+simulation_directory+"/build_ADDBACKGROUNDGRID", \
            "EXEC="+simulation_directory+"/Arepo_ADDBACKGRUNDGRID"])
if res != 0:
    sys.exit( np.int(res) )
    
## execute Arepo with ADDBACKGROUNDGRID from run directory
os.chdir(simulation_directory)
res = call(["mpiexec", "-np", "1","./Arepo_ADDBACKGRUNDGRID", "./param_ADDBACKGROUNDGRID.txt"])
if res != 0:
    sys.exit( np.int(res) )
os.chdir(cwd)


""" normal exit """
sys.exit(0) 
