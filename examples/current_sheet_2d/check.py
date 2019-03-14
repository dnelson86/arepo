""" @package ./examples/current_sheet_2d/check.py
Code that checks results of 2d current sheet problem

created by Rainer Weinberger, last modified 13.03.2019
"""

""" load libraries """
import sys    ## system calls
import numpy as np    ## load numpy
import h5py    ## load h5py; needed to read snapshots
import matplotlib.pyplot as plt

createFigures = True

simulation_directory = str(sys.argv[1])
print("examples/current_sheet_2d/check.py: checking simulation output in directory " + simulation_directory) 

FloatType = np.float64  # double precision: np.float64, for single use np.float32

## open initial conditiions to get parameters
try:
    data = h5py.File(simulation_directory + "/IC.hdf5", "r")
except:
    print("could not open initial  conditions!")
    exit(-1)
Boxsize = FloatType(data["Header"].attrs["BoxSize"])
NumberOfCells = np.int32(data["Header"].attrs["NumPart_Total"][0]) 
CellsPerDimension = np.sqrt(NumberOfCells) ## 2d sim

## parameters for initial state
density_0 = 1.0
velocity_radial_0 = -1.0    ## radial inflow velocity
pressure_0 = 1.0e-4
gamma = 5./3.  ## note: this has to be consistent with the parameter settings for Arepo!
utherm_0 = pressure_0 / ( gamma - 1.0 ) / density_0

## maximum L1 error after one propagation; based on tests
DeltaMaxAllowed = 0.1 * (FloatType(CellsPerDimension) / 150.0)**-1

Time = []
MagneticEnergy = []

""" loop over all output files """
i_file = 0
while True:
    """ try to read in snapshot """
    directory = simulation_directory+"/output/"
    filename = "snap_%03d.hdf5" % (i_file)
    try:
        data = h5py.File(directory+filename, "r")
    except:
        break
    
    """ get simulation data """
    
    B = np.array(data["PartType0"]["MagneticField"], dtype = FloatType)
    
    B2 = B[:,0]*B[:,0] + \
        B[:,1]*B[:,1] + \
        B[:,2]*B[:,2]
        
    Time.append( FloatType(data["Header"].attrs["Time"]) )
    MagneticEnergy.append( np.sum(B2) / 8.0 / np.pi )
    
    i_file += 1
    

Time = np.array(Time)
MagneticEnergy = np.array(MagneticEnergy)

if createFigures:
    fig, ax = plt.subplots(1)
    ax.plot(Time, MagneticEnergy/MagneticEnergy[0])
    fig.savefig(simulation_directory+'/MagneticEnergy.pdf')

## if everything is ok
if np.min(MagneticEnergy/MagneticEnergy[0]) > 0.92:
    sys.exit(0)
else:
    sys.exit(1)