""" @package ./examples/Yee_2d/check.py
Code that checks results of 2d Yee vortex problem

created by Rainer Weinberger, last modified 07.01.2018 -- comments welcome
"""

#### load libraries
import sys    ## load sys; needed for exit codes
import numpy as np    ## load numpy
import h5py    ## load h5py; needed to read snapshots
import os      # file specific calls
import matplotlib.pyplot as plt    ## needs to be active for plotting!
plt.rcParams['text.usetex'] = True

makeplots = True
if len(sys.argv) > 2:
  if sys.argv[2] == "True":
    makeplots = True
  else:
    makeplots = False

simulation_directory = str(sys.argv[1])
print("Yee_2d: checking simulation output in directory " + simulation_directory) 

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

#### parameters
Tinf = FloatType(1.0)
beta = FloatType(5.0)
gamma = FloatType(1.4)  ## note: this has to be consistent with the parameter settings for Arepo!

## maximum L1 error, empirically based
DeltaMaxAllowed = 0.005 * (FloatType(CellsPerDimension) / 50.0)**-2

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
    Time = FloatType(data["Header"].attrs["Time"])
    Pos = np.array(data["PartType0"]["CenterOfMass"], dtype = FloatType)
    VoronoiPos = np.array(data["PartType0"]["Coordinates"], dtype = FloatType)
    Density = np.array(data["PartType0"]["Density"], dtype = FloatType)
    Mass = np.array(data["PartType0"]["Masses"], dtype = FloatType)
    Velocity = np.array(data["PartType0"]["Velocities"], dtype = FloatType)
    Uthermal = np.array(data["PartType0"]["InternalEnergy"], dtype = FloatType)
    Volume = Mass / Density
    
    xPosFromCenter = (Pos[:,0] - 0.5 * Boxsize)
    yPosFromCenter = (Pos[:,1] - 0.5 * Boxsize)
    Radius = np.sqrt( xPosFromCenter**2 + yPosFromCenter**2 )
    RotationVelocity =  xPosFromCenter * Velocity[:,1] - yPosFromCenter * Velocity[:,0]
    i_select = np.where(Radius > 0)[0]
    RotationVelocity[i_select] /= Radius[i_select]

    """ calculate analytic solution """
    Temperature_ref = Tinf - ( ( gamma - 1.0 ) * beta * beta / ( 8.0 * np.pi * np.pi * gamma ) * np.exp( 1.0 - Radius*Radius ) )
    exponent = 1.0 / ( gamma - 1.0 )
    Density_ref = Temperature_ref**exponent
    Vel_ref = np.zeros([NumberOfCells, 3], dtype=FloatType)
    Vel_ref[:,0] = -0.5 * yPosFromCenter * beta / np.pi * np.exp( 0.5 * (1.0 - Radius * Radius) )
    Vel_ref[:,1] = 0.5 * xPosFromCenter * beta / np.pi * np.exp( 0.5 * (1.0 - Radius * Radius) )
    RotationVelocity_ref = xPosFromCenter * Vel_ref[:,1] - yPosFromCenter * Vel_ref[:,0]
    RotationVelocity_ref /= Radius
    Uthermal_ref = Temperature_ref / (gamma - 1.0)
    
    """ compare data """
    ## density
    abs_delta_dens = np.abs(Density - Density_ref)
    L1_dens = np.average(abs_delta_dens, weights=Volume)
    
    ## velocity, here, use absolute error (velocity can be close to zero!)
    abs_delta_vel = np.abs(RotationVelocity - RotationVelocity_ref)
    L1_vel = np.average(abs_delta_vel, weights=Volume)
    
    ## internal energy
    abs_delta_utherm = np.abs(Uthermal - Uthermal_ref)
    L1_utherm = np.average(abs_delta_utherm, weights=Volume)

    """ printing results """
    print("Yee_2d: L1 error of " + filename +":")
    print("\t density: %g" % L1_dens)
    print("\t velocity: %g" % L1_vel)
    print("\t specific internal energy: %g" % L1_utherm)
    print("\t tolerance: %g for %d cells per dimension" % (DeltaMaxAllowed, CellsPerDimension) )
    
    if makeplots:
      ## plot:
      fig = plt.figure( figsize=np.array([7.0,3.5]), dpi=300 )
      
      Nplot = 256
      from scipy import spatial # needed for KDTree that we use for nearest neighbour search
      Edges1d = np.linspace(0., Boxsize, Nplot+1, endpoint=True, dtype=FloatType)
      Grid1d = 0.5 * (Edges1d[1:] + Edges1d[:-1])
      xx, yy = np.meshgrid(Grid1d, Grid1d)
      Grid2D = np.array( [xx.reshape(Nplot**2), yy.reshape(Nplot**2)] ).T
      
      vor = spatial.Voronoi( VoronoiPos[:,:2] )
      dist, cells = spatial.KDTree( VoronoiPos[:,:2] ).query( Grid2D, k=1 )
      
      ax  = plt.axes( [0.05,0.15,0.35,0.7] )
      pc  = ax.pcolormesh( Edges1d, Edges1d, Density[cells].reshape((Nplot,Nplot)), rasterized=True )
      cax = plt.axes( [0.42,0.15,0.01,0.7] )
      plt.colorbar( pc, cax=cax )
      ax.set_title( 'Density' )
      spatial.voronoi_plot_2d( vor, ax=ax, line_colors='w', show_points=False, show_vertices=False, line_width=0.5 )
      ax.set_xlim( 0, Boxsize )
      ax.set_ylim( 0, Boxsize )
      
      error = (Density[cells]-Density_ref[cells])/Density_ref[cells]
      ax  = plt.axes( [0.53,0.15,0.35,0.7] )
      pc  = ax.pcolormesh( Edges1d, Edges1d, error.reshape((Nplot,Nplot)), rasterized=True )
      cax = plt.axes( [0.90,0.15,0.01,0.7] )
      plt.colorbar( pc, cax=cax )
      ax.set_title( 'Relative density error' )
      spatial.voronoi_plot_2d( vor, ax=ax, line_colors='w', show_points=False, show_vertices=False, line_width=0.5 )
      ax.set_xlim( 0, Boxsize )
      ax.set_ylim( 0, Boxsize )
      
      plt.text( 0.5, 0.92, "$t=%4.1f,\ L_1=%5.1e$" % (Time,L1_dens), ha='center', size=12, transform=fig.transFigure )

      if not os.path.exists( simulation_directory+"/plots" ):
        os.mkdir( simulation_directory+"/plots" )

      print(simulation_directory+"/plots/figure_%03d.pdf" % (i_file) )
      fig.savefig(simulation_directory+"/plots/figure_%03d.pdf" % (i_file), dpi=300)
      plt.close(fig)
    
    """ criteria for failing the test """
    if L1_dens > DeltaMaxAllowed or L1_vel > DeltaMaxAllowed or L1_utherm > DeltaMaxAllowed:
        sys.exit(1)
    i_file += 1
    
""" normal exit """
sys.exit(0)
    
