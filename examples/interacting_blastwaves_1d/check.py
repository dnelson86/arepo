""" @package ./examples/interacting_blastwaves_1d/check.py
Code that checks results of 1d interacting blastwaves

created by Rainer Weinberger, last modified: 04.03.2019
"""

""" load libraries """
import sys    # needed for exit codes
import numpy as np    # scientific computing package
import h5py    # hdf5 format
import os      # file specific calls
import matplotlib.pyplot as plt    ## needs to be active for plotting!
plt.rcParams['text.usetex'] = True

makeplots = True
if len(sys.argv) > 2:
  if sys.argv[2] == "True":
    makeplots = True
  else:
    makeplots = False

""" check functions """
def CheckL1Error(Pos, W, gamma, i_snap):
    """
    Compare the hydrodynamical quantities of the simulation with the high 
    resolution simulation solution, calculate the L1 error and check whether
    avarge L1 error is acceptable
    
    \param[in] Pos: shape (n,) array with 1d positions
    \param[in] W: shape (n,3) array with density, velocity and pressure of each cell
    \param[in] gamma: adiabatic index of gas
    \param[in] i_snap: index of snapshot
    
    \return status: zero if everything is within tolerances, otherwise 1
    """
    
    try:
        xx, rho, v, pres = np.loadtxt(simulation_directory+"/reference_%03d.txt" % i_snap ).T
    except:
        print("Check_L1_Error: could not find: %s/reference_%03d.txt" % (simulation_directory,i_snap))
        return 1
    W_exact = np.array([rho, v, pres]).T
    
    ## match indices
    dx = xx[1] - xx[0]
    index = np.array( (Pos - 0.5*dx) / dx, dtype=np.int32 )
    residual_x = Pos - xx[index]
    
    W_predict = np.zeros([len(index), 3], dtype=np.float64)
    for k in np.arange(3):
        W_predict[:,k] = (1.0- residual_x/dx) * W_exact[index,k] + residual_x/dx * W_exact[index+1,k]
    norm = np.abs(W_predict)
    norm[:,1] = 1 ## use absolute error in velocity component
    
    ## calculate L1 norm
    delta = np.abs(W_predict - W) / norm
    L1 = np.average(delta, axis=0)
    
    ## tolarance value; found empirically, fist order convergence!
    val_max = 0.05 * 400.0 / np.float(Pos.shape[0])
    L1MaxAllowed = np.array([val_max, 4.0*val_max, val_max])
    
    if np.any(L1 > L1MaxAllowed):
        print("CheckL1Error: ERROR: L1 error too large: %g %g %g; tolerance  %g %g %g"% (L1[0], L1[1], L1[2], L1MaxAllowed[0], L1MaxAllowed[1], L1MaxAllowed[2]) )
        return 1
    else:
        print("CheckL1Error: L1 error fine: %g %g %g; tolerance  %g %g %g"% (L1[0], L1[1], L1[2], L1MaxAllowed[0], L1MaxAllowed[1], L1MaxAllowed[2]) )
        return 0


def PlotSimulationData(Pos, W, gamma, i_snap, simulation_directory):
    """
    Plot density, velocity, specific internal energy and pressure of
    Simulation and high resolution solution result on top
    
    \param[in] Pos: shape (n,) array with 1d positions
    \param[in] W: shape (n,3) array with density, velocity and pressure of each cell
    \param[in] gamma: adiabatic index of gas
    \param[in] i_snap: index of snapshot
    \param[in] simulation_directory: path to simulation
    
    \return status: zero if everything went fine, 1: could not find exact solution data
    """

    ## set  limits for plot
    min = np.zeros(4)
    max = np.zeros(4)
    min[0] = np.min( W[:,0] ) - 0.3
    max[0] = np.max( W[:,0] ) + 1.0
    min[1] = np.min( W[:,1] ) - 0.3
    max[1] = np.max( W[:,1] ) + 1.0
    min[2] = np.min( W[:,2] / W[:,0] / (gamma - 1.0) ) - 30.0
    max[2] = np.max( W[:,2] / W[:,0] / (gamma - 1.0) ) + 30.0
    min[3] = np.min( W[:,2] ) - 30.0
    max[3] = np.max( W[:,2] ) + 30.0

    try:
        ## exact data, from 20,000 cell fixed grid simulation; downsampled by factor of 20
        xx, rho, v, pres = np.loadtxt(simulation_directory+"/reference_%03d.txt" % i_snap).T
    except:
        print("PlotSimulationData: could not find: %s/reference_%03d.txt" % (simulation_directory,i_snap))
        return 1
    ## plot:
    fig, ax = plt.subplots(4, sharex=True, figsize=np.array([6.9,6.0]))
    fig.subplots_adjust(left = 0.13, bottom = 0.09,right = 0.98, top = 0.98)
    
    ax[0].plot(Pos, W[:,0], 'r+', label="Arepo cells")
    ax[1].plot(Pos, W[:,1], 'r+')
    ax[2].plot(Pos, W[:,2]/W[:,0]/(gamma - 1.0), 'r+')
    ax[3].plot(Pos, W[:,2], 'r+')
    
    ax[0].plot(Pos, W[:,0], 'b', label="Arepo")
    ax[1].plot(Pos, W[:,1], 'b')
    ax[2].plot(Pos, W[:,2]/W[:,0]/(gamma - 1.0), 'b')
    ax[3].plot(Pos, W[:,2], 'b')
    
    ax[0].plot(xx, rho, color="k", lw=0.7, label="Exact solution")
    ax[1].plot(xx, v, color="k", lw=0.7)
    ax[2].plot(xx, pres / rho / (gamma - 1.0), color="k", lw=0.7)
    ax[3].plot(xx, pres, color="k", lw=0.7)
    
    for i_plot in np.arange(4):
        ax[i_plot].set_ylim([ min[i_plot], max[i_plot] ])
        ax[i_plot].set_xlim([0.5,0.9])
      
    ## set labels
    ax[3].set_xlabel(r"pos")  
    ax[0].set_ylabel(r"density")
    ax[1].set_ylabel(r"velocity")
    ax[2].set_ylabel(r"spec. int. energy")
    ax[3].set_ylabel(r"pressure")
    fig.align_ylabels(ax[:])

    if not os.path.exists( simulation_directory+"/plots" ):
      os.mkdir( simulation_directory+"/plots" )
    fig.savefig(simulation_directory+"/plots/figure_%03d.pdf" % (i_file))
    plt.close(fig)

    return 0

simulation_directory = str(sys.argv[1])
print("wave_1d: checking simulation output in directory " + simulation_directory) 

##parameters
Dtype = np.float64  # double precision: np.float64, for single use np.float32
gamma = 1.4  ## needs to be identical to Config.sh =(5./3.) if not specified there!

""" analyze snap_001.hdf5 """
i_file = 1
ReturnFlag = 0

## read in data
directory = simulation_directory+"/output/"
filename = "snap_%03d.hdf5" % (i_file)
print("try to open "+directory+filename)
## open hdf5 file
try:
    data = h5py.File(directory+filename, "r")
except:
    pass
print("analyzing "+ filename)

## get data from snapshot
time = np.float( data["Header"].attrs["Time"] )
position = np.array(data["PartType0"]["Coordinates"], dtype = Dtype)
density = np.array(data["PartType0"]["Density"], dtype = Dtype)
vel = np.array(data["PartType0"]["Velocities"], dtype = Dtype)
internalEnergy = np.array(data["PartType0"]["InternalEnergy"], dtype = Dtype)
## convert to more useful data structure
W = np.array([density, vel[:,0], (gamma-1.0)*internalEnergy*density], dtype=Dtype).T ## shape: (n,3)

""" plot data if you want """
if makeplots:
    ReturnFlag += PlotSimulationData(position[:,0], W, gamma, i_file, simulation_directory)

""" perform checks """
ReturnFlag += CheckL1Error(position[:,0], W, gamma, i_file)
      
if ReturnFlag == 0:
    print("interacting_blastwaves_1d: success!")
else:
    print("interacting_blastwaves_1d: failed!")

sys.exit(ReturnFlag)

