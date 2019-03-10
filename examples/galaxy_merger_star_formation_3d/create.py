""" @package examples/galaxy_merger_star_formation_3d/create.py
Code creates the output list; ICs need to be present already
galaxy_merger_star_formation_3d uses MakeNewDisk and CombineGalaxies codes
to create the initial conditions

created by Rainer Weinberger, last modified 09.03.2019
"""


""" load libraries """
import sys    # system calls
import numpy as np    # scientific computing package
import h5py    # hdf5 format
from subprocess import call    # execute shell commands 


""" input """
simulation_directory = str(sys.argv[1])
print("examples/collisionless_galaxy_3d/create.py " + simulation_directory)


""" set output times """
outputTimes = np.linspace(0.0,3.0,32, dtype=np.float64)
ones = np.ones(outputTimes.shape, dtype=np.int)


""" write output list file """
data = np.array([outputTimes, ones]).T
np.savetxt(simulation_directory+"/output_list.txt",data, fmt="%g %1.f" )


""" create backgroundgrid ICs from SPH ICs """
## compile Arepo with ADDBACKGROUNDGRID
res = call(["make", "CONFIG=./examples/galaxy_merger_star_formation_3d/Config_ADDBACKGROUNDGRID.sh", \
             "BUILD_DIR=./run/examples/galaxy_merger_star_formation_3d/build_ADDBACKGROUNDGRID", \
            "EXEC=./run/examples/galaxy_merger_star_formation_3d/Arepo_ADDBACKGRUNDGRID"])
if res != 0:
    sys.exit( np.int(res) )
    
## copy ICs to run directory
res = call(["cp", "./examples/galaxy_merger_star_formation_3d/ICs_1_1_merger_30_15_45_0_rmin10_start320_lowres.dat", \
            "./run/examples/galaxy_merger_star_formation_3d/"])
if res != 0:
    sys.exit( np.int(res) )

## execute Arepo with ADDBACKGROUNDGRID
res = call(["mpiexec", "-np", "1","./run/examples/galaxy_merger_star_formation_3d/Arepo_ADDBACKGRUNDGRID", "./examples/galaxy_merger_star_formation_3d/param_ADDBACKGROUNDGRID.txt"])
if res != 0:
    sys.exit( np.int(res) )

## clean up
res = call(["make", "CONFIG=./examples/galaxy_merger_star_formation_3d/Config_ADDBACKGROUNDGRID.sh", "clean"])
if res != 0:
    sys.exit( np.int(res) )
res = call(["rm", "./examples/galaxy_merger_star_formation_3d/param_ADDBACKGROUNDGRID.txt-usedvalues"])
if res != 0:
    sys.exit( np.int(res) )


""" normal exit """
sys.exit(0) 
