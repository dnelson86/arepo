""" @package ./examples/shocktube_1d/check
Code that checks results of 1d shocktube problem

created by Rainer Weinberger, last modified: 19.02.2019
"""

""" load libraries """
import sys    # needed for exit codes
import numpy as np    # scientific computing package
import h5py    # hdf5 format
import os      # file specific calls
import matplotlib.pyplot as plt    ## needs to be active for plotting!
plt.rcParams['text.usetex'] = True

from Riemann import *    ## Riemann-solver

makeplots = True
if len(sys.argv) > 2:
  if sys.argv[2] == "True":
    makeplots = True
  else:
    makeplots = False

forceExitOnError=False  ## exits immediately when tolerance is exceeded

""" check functiions """
def CheckL1Error(Pos, W, W_L, W_R, gamma, position_0, time):
    """
    Compare the hydrodynamical quantities of the simulation with the exact 
    solution of the Riemann problem, calculate the L1 error and check whether
    avarge L1 error is acceptable
    
    \param[in] Pos: shape (n,) array with 1d positions
    \param[in] W: shape (n,3) array with density, velocity and pressure of each cell
    \param[in] W_L: shape (3,) array with left initial condition state (density, velocity and pressure)
    \param[in] W_R: shape (3,) array with right initial condition state (density, velocity and pressure)
    \param[in] gamma: adiabatic index of gas
    \param[in] position_0: initial position of discontinuity
    \param[in] time: time of snapshot data
    
    \return status: zero if everything is within tolerances, otherwise 1
    """
    xx, W_exact, PosOfCharacteristics = RiemannProblem(Pos, position_0, W_L, W_R, gamma, time)
    
    norm = np.abs(W_exact)
    norm[:,1] = 1 ## use absolute error in velocity component
    
    ## calculate L1 norm
    delta = np.abs(W_exact - W) / norm
    L1 = np.average(delta, axis=0)
    
    ## tolarance value; found empirically, fist order convergence!
    L1MaxAllowed = 2.0 / np.float(Pos.shape[0])
    
    if np.any(L1 > L1MaxAllowed):
        print("CheckL1Error: ERROR: L1 error too large: %g %g %g; tolerance %g"% (L1[0], L1[1], L1[2], L1MaxAllowed) )
        return 1
    else:
        print("CheckL1Error: L1 error fine: %g %g %g; tolerance %g"% (L1[0], L1[1], L1[2], L1MaxAllowed) )
        return 0


def CheckTotalVariation(Pos, W, W_L, W_R, gamma, position_0, time):
    """
    Compare the total variation in simulation quantities with the total 
    variation in the analytic solution of the Riemann problem
    
    \param[in] Pos: shape (n,) array with 1d positions
    \param[in] W: shape (n,3) array with density, velocity and pressure of each cell
    \param[in] W_L: shape (3,) array with left initial condition state (density, velocity and pressure)
    \param[in] W_R: shape (3,) array with right initial condition state (density, velocity and pressure)
    \param[in] gamma: adiabatic index of gas
    \param[in] position_0: initial position of discontinuity
    \param[in] time: time of snapshot data
    
    \return status: zero if everything is within tolerances, otherwise 1
    """
    xx, W_exact, PosOfCharacteristics = RiemannProblem(Pos, position_0, W_L, W_R, gamma, time)
    
    TotalVariationSim = np.zeros(3)
    TotalVariationExact = np.zeros(3)
    
    i_sorted = np.argsort(Pos) ## sorted by position
    dW = W[i_sorted[1:],:] - W[i_sorted[:-1],:] ## difference of neighbouring cells
    dW_exact = W_exact[i_sorted[1:],:] - W_exact[i_sorted[:-1],:]

    for i in np.arange(3):  
        i1_pos,  = np.where(dW[:,i] >= 0)
        i1_neg,  = np.where(dW[:,i] < 0) 
        TotalVariationSim[i] = np.sum(dW[i1_pos, i], axis=0) - np.sum(dW[i1_neg, i])
        TotalVariationExact[i] = np.sum(dW_exact[i1_pos, i], axis=0) - np.sum(dW_exact[i1_neg, i])
        
    MaxRatioAllowed = 1.01
    if np.any( TotalVariationSim / TotalVariationExact > MaxRatioAllowed):
        print("CheckTotalVariation: ERROR: TotalVariation Sim/Exact: %g %g %g, tolerance: %g" % (TotalVariationSim[0] / TotalVariationExact[0], TotalVariationSim[1] / TotalVariationExact[1], TotalVariationSim[2] / TotalVariationExact[2], MaxRatioAllowed) )
        return 1
    else:
        print("CheckTotalVariation: TotalVariation Sim/Exact fine: %g %g %g, tolerance: %g" % (TotalVariationSim[0] / TotalVariationExact[0], TotalVariationSim[1] / TotalVariationExact[1], TotalVariationSim[2] / TotalVariationExact[2], MaxRatioAllowed) )
        return 0 


def CheckWidthOfDiscontinuities(Pos, W, W_L, W_R, gamma, position_0, time):
    """
    Measure the width of the fluid discontinuities in simulation quantities 
    to assess the numerical diffusivity of the scheme
    
    \param[in] Pos: shape (n,) array with 1d positions
    \param[in] W: shape (n,3) array with density, velocity and pressure of each cell
    \param[in] W_L: shape (3,) array with left initial condition state (density, velocity and pressure)
    \param[in] W_R: shape (3,) array with right initial condition state (density, velocity and pressure)
    \param[in] gamma: adiabatic index of gas
    \param[in] position_0: initial position of discontinuity
    \param[in] time: time of snapshot data
    
    \return status: zero if everything is within tolerances, otherwise 1
    """
    xx, W_exact, PosOfCharacteristics = RiemannProblem(Pos, position_0, W_L, W_R, gamma, time)
    
    for i, pos_char in enumerate(PosOfCharacteristics):
        ReturnFlag = 0
        ## PosOfCharacteristics will give different values in index 0 and 1 (3 and 4) 
        ## if it is a rarefaction; in this case, don't check, otherwise do check
        if i == 1 or i == 4:
            continue
        if i == 0 or i == 3:
            if PosOfCharacteristics[i] != PosOfCharacteristics[i+1]:
                continue

        ## get hydrodynamic states left and right of jump
        xx = np.array( [10.0+pos_char-1e-3, 10.0+pos_char+1e-3] )
        xx, W_exact, dummy = RiemannProblem(xx, 10.0, W_L, W_R, gamma, time)

        ## get threshold values to measure how many cells are within this interval
        jump = W_exact[1,:] - W_exact[0,:]
        percentiles = np.array([W_exact[0,:] + 0.05 * jump, W_exact[0,:] + 0.95 * jump])
        percentile_05 = np.min( percentiles, axis=0 )
        percentile_95 = np.max( percentiles, axis=0 )
        i_sorted = np.argsort( np.abs(Pos-pos_char-position_0) )
        
        i_low = np.full(3, i_sorted[0], dtype=np.int32)
        i_high = np.full(3, i_sorted[0], dtype=np.int32)
        
        ## check by how many cells 5th to 95th percentile of jump are sampled
        for j in np.arange(3):
            while W[i_low[j],j] > percentile_05[j] and W[i_low[j],j] < percentile_95[j]:
                i_low[j] -= 1
            while W[i_low[j],j] > percentile_05[j] and W[i_low[j],j] < percentile_95[j]:
                i_high[j] += 1
                
        ## sufficient for exact Riemann solver
        MaxNumerOfCells = 4
                
        if(i == 2):
            print("CheckWidthOfDiscontinuities: density jump at contact discontinuity resolved by %d cells (5th to 95th precentile), tolerance: %d" \
                  % (i_high[0]-i_low[0], MaxNumerOfCells) )
        else:
            print("CheckWidthOfDiscontinuities: density, velocity and pressure jump at shock resolved by %d, %d and %d cells (5th to 95th precentile), tolerance: %d" \
                  % (i_high[0]-i_low[0], i_high[1]-i_low[1], i_high[2]-i_low[2], MaxNumerOfCells) )
        if np.any(i_high-i_low > MaxNumerOfCells):
            print("CheckWidthOfDiscontinuities: ERROR: discontinuity too wide!")
            ReturnFlag += 1
    
    return ReturnFlag

def PlotSimulationData(Pos, W, W_L, W_R, gamma, position_0, time, simulation_directory):
    """
    Plot density, velocity, specific internal energy and pressure of
    Simulation and exact solution on top
    
    \param[in] Pos: shape (n,) array with 1d positions
    \param[in] W: shape (n,3) array with density, velocity and pressure of each cell
    \param[in] W_L: shape (3,) array with left initial condition state (density, velocity and pressure)
    \param[in] W_R: shape (3,) array with right initial condition state (density, velocity and pressure)
    \param[in] gamma: adiabatic index of gas
    \param[in] position_0: initial position of discontinuity
    \param[in] time: time of snapshot data
    \param[in] simulation_directory: path to simulation
    
    \return status: zero if everything went fine
    """
    xx = np.linspace(Pos.min(),Pos.max(),1000)
    xx, W_exact, dummy = RiemannProblem(xx, position_0, W_L, W_R, gamma, time)
    
    min = np.zeros(4)
    max = np.zeros(4)
    min[0] = np.min( W_exact[:,0] ) - 0.1
    max[0] = np.max( W_exact[:,0] ) + 0.1
    min[1] = np.min( W_exact[:,1] ) - 0.1
    max[1] = np.max( W_exact[:,1] ) + 0.1
    min[2] = np.min( W_exact[:,2] / W_exact[:,0] / (gamma - 1.0) ) - 0.1
    max[2] = np.max( W_exact[:,2] / W_exact[:,0] / (gamma - 1.0) ) + 0.1
    min[3] = np.min( W_exact[:,2] ) - 0.1
    max[3] = np.max( W_exact[:,2] ) + 0.1

    ## plot:
    fig, ax = plt.subplots(4, sharex=True, figsize=np.array([6.9,6.0]) )
    fig.subplots_adjust(left = 0.13, bottom = 0.09,right = 0.98, top = 0.98)
    
    ax[0].plot(xx, W_exact[:,0], color="k", lw=0.7, label="Exact solution")
    ax[1].plot(xx, W_exact[:,1], color="k", lw=0.7)
    ax[2].plot(xx, W_exact[:,2]/W_exact[:,0]/(gamma - 1.0), color="k", lw=0.7)
    ax[3].plot(xx, W_exact[:,2], color="k", lw=0.7)
    
    ax[0].plot(Pos, W[:,0], '+r', mec='r', mfc="None", label="Arepo cells")
    ax[1].plot(Pos, W[:,1], '+r', mec='r', mfc="None")
    ax[2].plot(Pos, W[:,2]/W[:,0]/(gamma - 1.0), '+r', mec='r', mfc="None")
    ax[3].plot(Pos, W[:,2], '+r', mec='r', mfc="None")
    
    ax[0].plot(Pos, W[:,0], 'b', label="Arepo")
    ax[1].plot(Pos, W[:,1], 'b')
    ax[2].plot(Pos, W[:,2]/W[:,0]/(gamma - 1.0), 'b')
    ax[3].plot(Pos, W[:,2], 'b')
    
    ax[3].set_xlim([0,20])
    for i_plot in np.arange(4):
        ax[i_plot].set_ylim([ min[i_plot], max[i_plot] ])
    
    ax[0].legend( loc='upper right', frameon=False, fontsize=8 )
    
    ## set labels
    ax[3].set_xlabel(r"pos")  
    ax[0].set_ylabel(r"density")
    ax[1].set_ylabel(r"velocity")
    ax[2].set_ylabel(r"spec. int. energy")
    ax[3].set_ylabel(r"pressure")
    fig.align_ylabels(ax[:])

    if not os.path.exists( simulation_directory+"/plots" ):
      os.mkdir( simulation_directory+"/plots" )

    print(simulation_directory+"/plots/figure_%03d.pdf" % (i_file) )
    fig.savefig(simulation_directory+"/plots/figure_%03d.pdf" % (i_file))
    plt.close(fig)

    return 0

simulation_directory = str(sys.argv[1])
print("wave_1d: checking simulation output in directory " + simulation_directory) 

Dtype = np.float64  # double precision: np.float64, for single use np.float32
## open initial conditiions to get parameters
directory = simulation_directory+"/"
filename = "IC.hdf5" 
try:
    data = h5py.File(directory+filename, "r")
except:
    print( "Could not find file {0}".format(directory+filename) )
    sys.exit(1)
IC_position = np.array(data["PartType0"]["Coordinates"], dtype = np.float64)
IC_mass = np.array(data["PartType0"]["Masses"], dtype = np.float64)
IC_velocity = np.array(data["PartType0"]["Velocities"], dtype = np.float64)
IC_internalEnergy = np.array(data["PartType0"]["InternalEnergy"], dtype = np.float64)
NumberOfCells = np.int32(data["Header"].attrs["NumPart_Total"][0])
DeltaMaxAllowed = 2.0 / np.float(NumberOfCells)

## get parameters of Riemann problem from ICs
gamma = 1.4  ## needs to be identical to Config.sh =(5./3.) if not specified there!
dx = IC_position[1,0]-IC_position[0,0] ## assuming a uniform grid
W_L = np.array([IC_mass[0]/dx, IC_velocity[0,0], (gamma - 1.0)*IC_internalEnergy[0]*IC_mass[0]/dx])
W_R = np.array([IC_mass[-1]/dx, IC_velocity[-1,0], (gamma - 1.0)*IC_internalEnergy[-1]*IC_mass[-1]/dx])
i_0, = np.where( (IC_mass[:]/dx == W_R[0]) & (IC_velocity[:,0] == W_R[1]) & ((gamma - 1.0)*IC_internalEnergy[-1]*IC_mass[-1]/dx == W_R[2]) )
position_0 = 0.5 * (IC_position[i_0[0]-1,0]+IC_position[i_0[0],0]) ## discontinuity at interface position

""" loop over all output files """
i_file = -1
ReturnFlag = 0
while True:
    i_file += 1
    
    ## read in data
    directory = simulation_directory+"/output/"
    filename = "snap_%03d.hdf5" % (i_file)
    print("try to open "+directory+filename)
    ## open hdf5 file
    try:
        data = h5py.File(directory+filename, "r")
    except:
        break
    print("analyzing "+ filename)
    
    ## get data from snapshot
    time = np.float( data["Header"].attrs["Time"] )
    position = np.array(data["PartType0"]["Coordinates"], dtype = np.float64)
    density = np.array(data["PartType0"]["Density"], dtype = np.float64)
    vel = np.array(data["PartType0"]["Velocities"], dtype = np.float64)
    internalEnergy = np.array(data["PartType0"]["InternalEnergy"], dtype = np.float64)
    ## convert to more useful data structure
    W = np.array([density, vel[:,0], (gamma-1.0)*internalEnergy*density], dtype=np.float64).T ## shape: (n,3)
    
    
    """ plot data if you want """
    if makeplots:
        ReturnFlag += PlotSimulationData(position[:,0], W, W_L, W_R, gamma, position_0, time, simulation_directory)
        if ReturnFlag > 0 and forceExitOnError:
            print('ERROR: something went wrong in plot')
            sys.exit(ReturnFlag)
    
    """ perform checks """
    ReturnFlag += CheckL1Error(position[:,0], W, W_L, W_R, gamma, position_0, time)
    if ReturnFlag > 0 and forceExitOnError:
        print('ERROR: exceeding tolerance!')
        sys.exit(ReturnFlag)
        
    ReturnFlag += CheckTotalVariation(position[:,0], W, W_L, W_R, gamma, position_0, time)
    if ReturnFlag > 0 and forceExitOnError:
        print('ERROR: exceeding tolerance!')
        sys.exit(ReturnFlag)
    
    ReturnFlag += CheckWidthOfDiscontinuities(position[:,0], W, W_L, W_R, gamma, position_0, time)
    if ReturnFlag > 0 and forceExitOnError:
        print('ERROR: exceeding tolerance!')
        sys.exit(ReturnFlag)
        
if ReturnFlag == 0:
    print("check.py: success!")
else:
    print("check.py: failed!")

sys.exit(ReturnFlag)


