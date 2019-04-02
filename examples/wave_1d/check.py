""" @package ./examples/wave_1d/check.py
Code that checks results of 1d wave propagation problem

created by Rainer Weinberger, last modified 19.2.2019 -- comments welcome
"""

""" load libraries """
import sys    # system specific calls
import numpy as np    # scientific computing package
import h5py    # hdf5 format
import os     # file specific calls

simulation_directory = str(sys.argv[1])
print("wave_1d: checking simulation output in directory " + simulation_directory) 

FloatType = np.float64  # double precision: np.float64, for single use np.float32
IntType = np.int32 # integer type

makeplots = True
if len(sys.argv) > 2:
  if sys.argv[2] == "True":
    makeplots = True
  else:
    makeplots = False

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
    Pos = np.array(data["PartType0"]["CenterOfMass"], dtype = FloatType)
    Density = np.array(data["PartType0"]["Density"], dtype = FloatType)
    Mass = np.array(data["PartType0"]["Masses"], dtype = FloatType)
    Velocity = np.array(data["PartType0"]["Velocities"], dtype = FloatType)
    Uthermal = np.array(data["PartType0"]["InternalEnergy"], dtype = FloatType)
    Volume = Mass / Density

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
    L1_dens = np.average(abs_delta_dens, weights=Volume)
    
    ## velocity, here, use absolute error (velocity should be zero! check only x-vel, taking all components dilutes the norm!)
    abs_delta_vel = np.abs(Velocity - Velocity_ref)[:,0]
    L1_vel = np.average(abs_delta_vel, weights=Volume)
    
    ## internal energy
    abs_delta_utherm = np.abs(Uthermal - Uthermal_ref) / Uthermal_ref
    L1_utherm = np.average(abs_delta_utherm, weights=Volume)

    """ printing results """
    print("wave_1d: L1 error of " + filename +":")
    print("\t density: %g" % L1_dens)
    print("\t velocity: %g" % L1_vel)
    print("\t specific internal energy: %g" % L1_utherm)
    print("\t tolerance: %g for %d cells" % (DeltaMaxAllowed, NumberOfCells) )
    
    """ criteria for failing the test """
    if L1_dens > DeltaMaxAllowed or L1_vel > DeltaMaxAllowed or L1_utherm > DeltaMaxAllowed:
        sys.exit(1)
    
    if makeplots and i_file > 0:
      if not os.path.exists( simulation_directory+"/plots" ):
        os.mkdir( simulation_directory+"/plots" )
      
      # only import matplotlib if needed
      import matplotlib.pyplot as plt
      plt.rcParams['text.usetex'] = True
      f = plt.figure( figsize=(3.5,3.5) )
      ax = plt.axes( [0.19, 0.12, 0.75, 0.75] )
      refpos = np.linspace( 0, Boxsize, 257 )
      ax.plot( refpos, density_0 * (1. + delta * np.sin(2.0 * np.pi * refpos / Boxsize)), 'k', lw=0.7, label="Analytical solution" )
      ax.plot( Pos[:,0], Density, 'o-r', mec='r', mfc="None", label="Arepo" )
      ax.set_xlim( 0, Boxsize )
      ax.set_xlabel( "x" )
      ax.set_ylabel( "Density deviation" )
      ax.legend( loc='upper right', frameon=False, fontsize=8 )
      ax.set_title( "$\mathrm{wave\_1d:}\ \mathrm{N}=%d,\ \mathrm{L1}=%4.1e$" % (NumberOfCells,L1_dens), loc='right', size=8 )
      plt.ticklabel_format( axis='y', style='sci', scilimits=(0,0) )
      f.savefig( simulation_directory+"plots/density.pdf" )
      
      f = plt.figure( figsize=(3.5,3.5) )
      ax = plt.axes( [0.19, 0.12, 0.75, 0.75] )
      refpos = np.linspace( 0, Boxsize, 257 )
      ax.plot( refpos, np.ones(257)*velocity_0, 'k', lw=0.7, label="Analytical solution" )
      ax.plot( Pos[:,0], Velocity[:,0], 'o-r', mec='r', mfc="None", label="Arepo" )
      ax.set_xlim( 0, Boxsize )
      ax.set_xlabel( "x" )
      ax.set_ylabel( "Velocity" )
      ax.legend( loc='lower center', frameon=False, fontsize=8 )
      ax.set_title( "$\mathrm{wave\_1d:}\ \mathrm{N}=%d,\ \mathrm{L1}=%4.1e$" % (NumberOfCells,L1_vel), loc='right', size=8 )
      plt.ticklabel_format( axis='y', style='sci', scilimits=(0,0) )
      f.savefig( simulation_directory+"plots/velocity.pdf" )
    
    i_file += 1

""" normal exit """
sys.exit(0) 
