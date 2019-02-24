""" @package ./examples/noh_3d/check.py
Code that checks results of 3d Noh problem

created by Rainer Weinberger, last modified 24.02.2019
"""

""" load libraries """
import sys    ## load sys; needed for exit codes
import numpy as np    ## load numpy
import h5py    ## load h5py; needed to read snapshots
import matplotlib.pyplot as plt    ## plot stuff

createFigures = True

simulation_directory = str(sys.argv[1])
print("examples/noh_3d/check.py: checking simulation output in directory " + simulation_directory)

FloatType = np.float64  # double precision: np.float64, for single use np.float32
IntType = np.int32

## open initial conditiions to get parameters
try:
    data = h5py.File(simulation_directory + "/IC.hdf5", "r")
except:
    print("could not open initial  conditions!")
    exit(1)
Boxsize = FloatType(data["Header"].attrs["BoxSize"])
NumberOfCells = IntType(data["Header"].attrs["NumPart_Total"][0]) 
CellsPerDimension = NumberOfCells**(1./3.) ## 3d sim

## parameters for initial state
density_0 = 1.0
velocity_radial_0 = -1.0    ## radial inflow velocity
pressure_0 = 1.0e-4
gamma = 5./3.  ## note: this has to be consistent with the parameter settings for Arepo!
utherm_0 = pressure_0 / ( gamma - 1.0 ) / density_0

## maximum L1 error after one propagation; empirically based
DeltaMaxAllowed = 0.25

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
    time = FloatType(data["Header"].attrs["Time"])
    ## simulation data
    Pos = np.array(data["PartType0"]["Coordinates"], dtype = FloatType)
    Density = np.array(data["PartType0"]["Density"], dtype = FloatType)
    Mass = np.array(data["PartType0"]["Masses"], dtype = FloatType)
    Velocity = np.array(data["PartType0"]["Velocities"], dtype = FloatType)
    Uthermal = np.array(data["PartType0"]["InternalEnergy"], dtype = FloatType)
    CellVolume = Mass / Density
    
    xPosFromCenter = (Pos[:,0] - 0.5 * Boxsize)
    yPosFromCenter = (Pos[:,1] - 0.5 * Boxsize)
    zPosFromCenter = (Pos[:,2] - 0.5 * Boxsize)
    Radius = np.sqrt( xPosFromCenter**2 + yPosFromCenter**2 + zPosFromCenter**2 )
    
    vRad = Velocity[:,0] * xPosFromCenter / Radius \
      + Velocity[:,1] * yPosFromCenter / Radius \
      + Velocity[:,2] * zPosFromCenter / Radius
    
    r_shock = 1./3. * time
    
    ## exact solution, 2d
    Radius[Radius == 0] = 1e-10
    Density_ref = density_0 * (1.0 + time/Radius) * (1.0 + time/Radius) ## in 3d
    i_postshock, = np.where(Radius < r_shock)
    Density_ref[i_postshock] = 64.0    ## for 3d
    
    #### plot profiles
    if createFigures:
        fig, ax = plt.subplots(3, sharex=True, figsize=np.array([6.9,6.0]) )
        fig.subplots_adjust(left = 0.13, bottom = 0.09,right = 0.98, top = 0.98)
        
        ax[0].scatter(Radius, Density, rasterized=True, s=0.2)
        ax[1].scatter(Radius, vRad, rasterized=True, s=0.2)
        ax[2].scatter(Radius, Uthermal, rasterized=True, s=0.2)
        
        i_sorted = np.argsort(Radius)
        ax[0].plot(Radius[i_sorted], Density_ref[i_sorted], label="analytic solution")
        ax[0].set_xlim(0,0.8)
        ax[0].set_ylim(0,70)
        ax[2].set_xlabel(r"radius")
        ax[0].set_ylabel(r"density")
        ax[1].set_ylabel(r"radial velocity")
        ax[2].set_ylabel(r"spec. internal energy")
        
        fig.align_ylabels(ax[:])
        fig.savefig(directory+"/figure_%03d.pdf"%i_file)

    #### check against analytic solution
    i_compare, = np.where(Radius < 0.8) ## only check inner region; boundary has spourious effects in this testcase
    abs_delta_dens = np.abs(Density[i_compare] - Density_ref[i_compare]) / Density_ref[i_compare]
    L1_dens = np.average(abs_delta_dens, weights=CellVolume[i_compare] )
    
    L1_max = DeltaMaxAllowed * time
    print("examples/Noh_3d/check snapshot %d: L1_dens = %g, DeltaMaxAllowed = %g!" % (i_file, L1_dens, L1_max) )
    
    
    if L1_dens > L1_max:
        sys.exit(1)
    
    i_file += 1
    
## if everything is ok
sys.exit(0)
