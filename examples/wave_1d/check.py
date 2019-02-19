""" @package ./examples/wave_1d/check.py
Code that checks results of 1d wave propagation problem

created by Rainer Weinberger, last modified 19.2.2019 -- comments welcome
"""

""" load libraries """
import sys    # system specific calls
import numpy as np    # scientific computing package
import h5py    # hdf5 format

simulation_directory = str(sys.argv[1])
print("examples/wave_1d/check.py: checking simulation output in directory " + simulation_directory) 

FloatType = np.float64  # double precision: np.float64, for single use np.float32
IntType = np.int32 # integer type

""" open initial conditiions to get parameters """
try:
    data = h5py.File(simulation_directory + "/IC.hdf5", "r")
except:
    print("could not open initial  conditions!")
    exit(-1)
Boxsize = FloatType(data["Header"].attrs["BoxSize"])
NumberOfCells = np.int32(data["Header"].attrs["NumPart_Total"][0])

""" maximum L1 error after one propagation; empirically motivated """
DeltaMaxAllowed = 5e-5 * FloatType(NumberOfCells)**-2

""" initial state -- copied from create.py """
density_0 = FloatType(1.0)
velocity_0 = FloatType(0.0)
pressure_0 = FloatType(3.0) / FloatType(5.0)
gamma = FloatType(5.0) / FloatType(3.0)
gamma_minus_one = gamma - FloatType(1.0)
delta = FloatType(1e-6)    # relative density perturbation
uthermal_0 = pressure_0 / density_0 / gamma_minus_one

""" 
    loop over all output files; need to be at times when analytic
    solution equals the initial conditions
"""
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
    Velocity = np.array(data["PartType0"]["Velocities"], dtype = FloatType)
    Uthermal = np.array(data["PartType0"]["InternalEnergy"], dtype = FloatType)

    """ calculate analytic solution at new cell positions """
    Density_ref = np.full(Pos.shape[0], density_0, dtype=FloatType)
    Velocity_ref = np.zeros(Pos.shape, dtype=FloatType)
    Uthermal_ref = np.full(Pos.shape[0], uthermal_0, dtype=FloatType)
    ## perturbations
    Density_ref *=  FloatType(1.0) + delta * np.sin( FloatType(2.0) * FloatType(np.pi) * Pos[:,0] / Boxsize )
    Velocity_ref[:,0] = velocity_0
    Uthermal_ref *= (Density / density_0)**gamma_minus_one

    """ compare data """
    ## density
    abs_delta_dens = np.abs(Density - Density_ref) / Density_ref
    L1_dens = np.average(abs_delta_dens)
    
    ## velocity, here, use absolute error (velocity should be zero!)
    abs_delta_vel = np.abs(Velocity - Velocity_ref)
    L1_vel = np.average(abs_delta_vel)
    
    ## internal energy
    abs_delta_utherm = np.abs(Uthermal - Uthermal_ref) / Uthermal_ref
    L1_utherm = np.average(abs_delta_utherm)

    """ printing results """
    print("examples/wave_1d/check.py: L1 error of " + filename +":")
    print("\t density: %g" % L1_dens)
    print("\t velocity: %g" % L1_vel)
    print("\t specific internal energy: %g" % L1_utherm)
    print("\t tolerance: %g for %d cells" % (DeltaMaxAllowed, NumberOfCells) )
    
    """ criteria for failing the test """
    if L1_dens > DeltaMaxAllowed or L1_vel > DeltaMaxAllowed or L1_utherm > DeltaMaxAllowed:
        sys.exit(1)
    i_file += 1

""" normal exit """
sys.exit(0) 
