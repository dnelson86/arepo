""" @package ./examples/Noh_2d/check.py
Code that checks results of 2d Noh problem

created by Rainer Weinberger, last modified 23.02.2019
"""

""" load libraries """
import sys    # system specific calls
import numpy as np    ## load numpy
import h5py    ## load h5py; needed to read snapshots
import os      # file specific calls
import matplotlib.pyplot as plt    # needs to be active for plotting!
plt.rcParams['text.usetex'] = True

makeplots = True
if len(sys.argv) > 2:
  if sys.argv[2] == "True":
    makeplots = True
  else:
    makeplots = False

simulation_directory = str(sys.argv[1])
print("Noh_2d: checking simulation output in directory " + simulation_directory) 

FloatType = np.float64  # double precision: np.float64, for single use np.float32
IntType = np.int32

## open initial conditiions to get parameters
try:
    data = h5py.File(simulation_directory + "/IC.hdf5", "r")
except:
    print("could not open initial  conditions!")
    exit(-1)
Boxsize = FloatType(data["Header"].attrs["BoxSize"])
NumberOfCells = IntType(data["Header"].attrs["NumPart_Total"][0]) 
CellsPerDimension = np.sqrt(NumberOfCells) ## 2d sim

## parameters for initial state
density_0 = 1.0
velocity_radial_0 = -1.0    ## radial inflow velocity
pressure_0 = 1.0e-4
gamma = 5./3.  ## note: this has to be consistent with the parameter settings for Arepo!
utherm_0 = pressure_0 / ( gamma - 1.0 ) / density_0

## maximum L1 error after one propagation; empirically based for 150**2 cells
DeltaMaxAllowed = 0.05

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
    Pos = np.array(data["PartType0"]["CenterOfMass"], dtype = FloatType)
    VoronoiPos = np.array(data["PartType0"]["Coordinates"], dtype = FloatType)
    Density = np.array(data["PartType0"]["Density"], dtype = FloatType)
    Mass = np.array(data["PartType0"]["Masses"], dtype = FloatType)
    Velocity = np.array(data["PartType0"]["Velocities"], dtype = FloatType)
    Uthermal = np.array(data["PartType0"]["InternalEnergy"], dtype = FloatType)
    CellVolume = Mass / Density
    
    xPosFromCenter = (Pos[:,0] - 0.5 * Boxsize)
    yPosFromCenter = (Pos[:,1] - 0.5 * Boxsize)
    Radius = np.sqrt( xPosFromCenter**2 + yPosFromCenter**2 )
    
    vRad = Velocity[:,0] * xPosFromCenter / Radius + Velocity[:,1] * yPosFromCenter / Radius
    
    r_shock = 1./3. * time
    
    ## exact solution, 2d
    Radius[Radius == 0] = 1e-10
    Density_ref = density_0 * (1.0 + time/Radius)
    i_postshock, = np.where(Radius < r_shock)
    Density_ref[i_postshock] = 16.0    ## for 2d
    
    #### plot density profile
    if makeplots:
        fig = plt.figure( figsize=np.array([7.5,6.0]), dpi=300 )
        i_sorted = np.argsort(Radius)
        
        if r_shock == 0:
          r_shock = 0.1 * Boxsize
        
        Nplot = 512
        from scipy import spatial # needed for KDTree that we use for nearest neighbour search and Voronoi mesh
        Edges1d = np.linspace(0.5*Boxsize-1.3*r_shock, 0.5*Boxsize+1.4*r_shock, Nplot+1, endpoint=True, dtype=FloatType)
        Grid1d = 0.5 * (Edges1d[1:] + Edges1d[:-1])
        xx, yy = np.meshgrid(Grid1d, Grid1d)
        Grid2D = np.array( [xx.reshape(Nplot**2), yy.reshape(Nplot**2)] ).T      
        dist, cells = spatial.KDTree( VoronoiPos[:,:2] ).query( Grid2D, k=1 )
        
        ax = plt.axes( [0.08,0.70,0.52,0.25] )
        ax.plot(Radius[i_sorted], Density_ref[i_sorted], 'k', label="Analytic solution")
        ax.plot(Radius, Density, 'b.', rasterized=True, ms=0.2)
        ax.set_ylabel(r"Density")
        ax.set_ylim( 0, 20 )
      
        ax  = plt.axes( [0.65,0.70,0.20,0.25] )
        pc  = ax.pcolormesh( Edges1d, Edges1d, Density[cells].reshape((Nplot,Nplot)), rasterized=True, cmap=plt.get_cmap('viridis'), vmin=0, vmax=20. )
        cax = plt.axes( [0.88,0.70,0.02,0.25] )
        plt.colorbar( pc, cax=cax )
        ax.set_xlim( 0.5*Boxsize-1.3*r_shock, 0.5*Boxsize+1.3*r_shock )
        ax.set_ylim( 0.5*Boxsize-1.3*r_shock, 0.5*Boxsize+1.3*r_shock )
       
        ax = plt.axes( [0.08,0.38,0.52,0.25] )
        ax.plot(Radius, vRad, 'b.', rasterized=True, ms=0.2)
        ax.set_ylabel(r"Radial velocity")
        ax.set_ylim( -1.1, 1.1 )
        
        ax  = plt.axes( [0.65,0.38,0.20,0.25] )
        pc  = ax.pcolormesh( Edges1d, Edges1d, vRad[cells].reshape((Nplot,Nplot)), rasterized=True, cmap=plt.get_cmap('plasma'), vmin=-1.1, vmax=1.1 )
        cax = plt.axes( [0.88,0.38,0.02,0.25] )
        plt.colorbar( pc, cax=cax )
        ax.set_xlim( 0.5*Boxsize-1.3*r_shock, 0.5*Boxsize+1.3*r_shock )
        ax.set_ylim( 0.5*Boxsize-1.3*r_shock, 0.5*Boxsize+1.3*r_shock )
        
        ax = plt.axes( [0.08,0.07,0.52,0.25] )
        ax.plot(Radius, Uthermal, 'b.', rasterized=True, ms=0.2)
        ax.set_xlabel(r"Radius")
        ax.set_ylabel(r"Spec. internal energy")
        ax.set_ylim( -0.05, 0.8 )
        
        ax  = plt.axes( [0.65,0.07,0.20,0.25] )
        pc  = ax.pcolormesh( Edges1d, Edges1d, Uthermal[cells].reshape((Nplot,Nplot)), rasterized=True, cmap=plt.get_cmap('magma'), vmin=-0.05, vmax=0.8 )
        cax = plt.axes( [0.88,0.07,0.02,0.25] )
        plt.colorbar( pc, cax=cax )
        ax.set_xlim( 0.5*Boxsize-1.3*r_shock, 0.5*Boxsize+1.3*r_shock )
        ax.set_ylim( 0.5*Boxsize-1.3*r_shock, 0.5*Boxsize+1.3*r_shock )
        
        if not os.path.exists( simulation_directory+"/plots" ):
          os.mkdir( simulation_directory+"/plots" )

        print(simulation_directory+"/plots/figure_%03d.pdf" % (i_file) )
        fig.savefig(simulation_directory+"/plots/figure_%03d.pdf" % (i_file), dpi=300)
        plt.close(fig)

    #### check against analytic solution
    i_compare, = np.where(Radius < 0.8) ## only check inner region; boundary has spourious effects in this testcase
    abs_delta_dens = np.abs(Density[i_compare] - Density_ref[i_compare]) / Density_ref[i_compare]
    L1_dens = np.average(abs_delta_dens, weights=CellVolume[i_compare] )
    
    L1_max = DeltaMaxAllowed * time
    print("Noh_2d: snapshot %d: DEBUG: L1_dens = %g, DeltaMaxAllowed = %g!" % (i_file, L1_dens, L1_max) )
    
    if L1_dens > L1_max and time > 0.:
        print("Noh_2d: ERROR: snaphshot %d: L1_dens = %g, DeltaMaxAllowed = %g!" % (i_file, L1_dens, L1_max) )
        sys.exit(1)
    
    i_file += 1
    

## if everything is ok
sys.exit(0)
