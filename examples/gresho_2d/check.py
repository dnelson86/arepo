""" @package ./examples/gresho_2d/check.py
Code that checks results of 2d Gresho vortex problem

created by Rainer Weinberger, last modified 21.02.2019 -- comments welcome
"""

#### load libraries
import sys    # system specific calls
import numpy as np    ## load numpy
import h5py    ## load h5py; needed to read snapshots

simulation_directory = str(sys.argv[1])
print("examples/gresho_2d/check.py: checking simulation output in directory " + simulation_directory)

FloatType = np.float64  # double precision: np.float64, for single use np.float32

## open initial conditiions to get parameters
try:
    data = h5py.File(simulation_directory + "/IC.hdf5", "r")
except:
    print("could not open initial  conditions!")
    exit(-1)
Boxsize = FloatType(data["Header"].attrs["BoxSize"])
NumberOfCells = np.int32(data["Header"].attrs["NumPart_Total"][0]) 
CellsPerDimension = np.sqrt(NumberOfCells)

## parameters
density_0 = 1.0
velocity_0 = 3.0
gamma = 5.0/3.0
gamma_minus_one = gamma - 1.0

## maximum L1 error after one propagation; empirically based
DeltaMaxAllowed = 0.03 * (FloatType(CellsPerDimension) / 40.0)**-1.4

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
    ## simulation data
    Pos = np.array(data["PartType0"]["Coordinates"], dtype = FloatType)
    Density = np.array(data["PartType0"]["Density"], dtype = FloatType)
    Mass = np.array(data["PartType0"]["Masses"], dtype = FloatType)
    Velocity = np.array(data["PartType0"]["Velocities"], dtype = FloatType)
    Velocity[:,0] -= velocity_0
    Uthermal = np.array(data["PartType0"]["InternalEnergy"], dtype = FloatType)
    Volume = Mass / Density
    
    xPosFromCenter = (Pos[:,0] - 0.5 * Boxsize)
    yPosFromCenter = (Pos[:,1] - 0.5 * Boxsize)
    Radius = np.sqrt( xPosFromCenter**2 + yPosFromCenter**2 )
    i_select, = np.where(Radius != 0.0)
    RotationVelocity = yPosFromCenter * Velocity[:,0] - xPosFromCenter * Velocity[:,1]
    RotationVelocity[i_select] /= Radius[i_select]
    RotationVelocity[Radius[:] == 0.0] = 0.0

    """ calculate analytic solution """
    ## density
    Density_ref = np.full(Density.shape, density_0, dtype=FloatType)
    
    ## different zones
    i1, = np.where( Radius < 0.2 )
    i2, = np.where( (Radius >= 0.2) & (Radius < 0.4) )
    i3, = np.where( Radius >= 0.4 )
    
    ## velocity
    RotationVelocity_ref = np.zeros(Velocity.shape[0], dtype=FloatType)
    RotationVelocity_ref[i1] = 5.0 * Radius[i1]
    RotationVelocity_ref[i2] = 2.0 - 5.0 * Radius[i2]
    RotationVelocity_ref[i3] = 0.0
    Velocity_ref = np.zeros(Velocity.shape, dtype=FloatType )
    i_select, = np.where(Radius != 0.0)
    Velocity_ref[i_select,0] = RotationVelocity_ref[i_select] * (Pos[i_select,1] - 0.5 * Boxsize) / Radius[i_select]
    Velocity_ref[i_select,1] = -RotationVelocity_ref[i_select] * (Pos[i_select,0] - 0.5 * Boxsize) / Radius[i_select]
    
    ## specific internal energy
    Pressure_ref = np.zeros(Uthermal.shape, dtype=FloatType)
    Pressure_ref[i1] = 5.0 + 12.5 * Radius[i1] * Radius[i1]
    Pressure_ref[i2] = 9.0 + 12.5 * Radius[i2] * Radius[i2] - 20.0 * Radius[i2] + 4.0 * np.log(Radius[i2] / 0.2)
    Pressure_ref[i3] = 3.0 + 4.0 * np.log(2.0)
    Uthermal_ref = Pressure_ref / density_0 / gamma_minus_one
    
    """ compare data """
    ## density
    abs_delta_dens = np.abs(Density - Density_ref)
    L1_dens = np.average(abs_delta_dens, weights = Volume)
    
    ## velocity
    abs_delta_vel = np.abs(RotationVelocity - RotationVelocity_ref)
    L1_vel = np.average(abs_delta_vel, weights = Volume)
    
    ## internal energy
    abs_delta_utherm = np.abs(Uthermal - Uthermal_ref)
    L1_utherm = np.average(abs_delta_utherm, weights = Volume)

    """ printing results """
    print("examples/Gresho_2d/check.py: L1 error of " + filename +":")
    print("\t density: %g" % L1_dens)
    print("\t velocity: %g" % L1_vel)
    print("\t specific internal energy: %g" % L1_utherm)
    print("\t tolerance: %g for %d cells per dimension" % (DeltaMaxAllowed, CellsPerDimension) )
    
    """ criteria for failing the test """
    if L1_dens > DeltaMaxAllowed or L1_vel > DeltaMaxAllowed or L1_utherm > DeltaMaxAllowed:
        sys.exit(-1)
    i_file += 1
    
print("normal exit")
""" normal exit """
sys.exit(0)
