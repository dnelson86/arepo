""" @package ./examples/cosmo_box_star_formation_3d/check.py
Code that checks results of structure formation simulation

created by Rainer Weinberger, last modified 28.02.2019
"""

""" load libraries """
import sys    ## system calls
import numpy as np    ## load numpy
import h5py    ## load h5py; needed to read snapshots
import os      # file specific calls
import matplotlib.pyplot as plt    # needs to be active for plotting!
plt.rcParams['text.usetex'] = True

simulation_directory = str(sys.argv[1])
print("cosmo_box_star_formation_3d: checking simulation output in directory " + simulation_directory) 

""" script settings """
createReferenceSolution = False
CompareToReferenceRun = True

makeplots = True
if len(sys.argv) > 2:
  if sys.argv[2] == "True":
    makeplots = True
  else:
    makeplots = False

FloatType = np.float64  # double precision: np.float64, for single use np.float32

Boxsize = 7.5 ##Mpc/h
CellsPerDimension = 32 ## 3d sim
NumberOfCells = CellsPerDimension * CellsPerDimension * CellsPerDimension
HubbleParam = 0.6774
UnitMass = 1.0e10/HubbleParam
Volume = (Boxsize/HubbleParam) * (Boxsize/HubbleParam) * (Boxsize/HubbleParam)


""" read in star formation rate output snapshot """
directory = simulation_directory+"/output/"
filename = "sfr.txt"

data = np.loadtxt(directory+filename)

if createReferenceSolution:
    np.savetxt(simulation_directory+"/sfr_ref_L8n32_reduced.txt", data[::50,:])

scalefactor = data[:,0]
redshift_plus_one = 1.0 / data[:,0]
sfrDensity = data[:,3] / Volume
cum_mass_stars = data[:,5] * UnitMass / Volume

""" compare to reference output """
data_ref = np.loadtxt(simulation_directory+"/sfr_ref_L8n32_reduced.txt")
scalefactor_ref = data_ref[:,0]
redshift_plus_one_ref = 1.0 / data_ref[:,0]
cum_mass_stars_ref = data_ref[:,5]* UnitMass / Volume

scalefactor_to_probe = np.array([0.25,0.3333,0.5,0.66667,1])
delta_mass_stars = np.zeros(scalefactor_to_probe.shape)
avg_mass_stars = np.zeros(scalefactor_to_probe.shape)

if makeplots:
    if not os.path.exists( simulation_directory+"/plots" ):
        os.mkdir( simulation_directory+"/plots" )
    
    ## figure -- dark matter and stellar positions
    fig, ax = plt.subplots( 2, 2, figsize=np.array([6.9,6.9]), sharex=True, sharey=True )
    fig.subplots_adjust(left = 0.09, bottom = 0.09,right = 0.98, top = 0.98, hspace=0.0, wspace=0.0)
    
    for i, a in enumerate(scalefactor_to_probe):
        delta_mass_stars[i]  = np.interp(a, scalefactor, cum_mass_stars)
        delta_mass_stars[i] -= np.interp(a, scalefactor_ref, cum_mass_stars_ref)
        avg_mass_stars[i]  = 0.5 * np.interp(a, scalefactor, cum_mass_stars)
        avg_mass_stars[i] += 0.5 * np.interp(a, scalefactor_ref, cum_mass_stars_ref)
    
        filename = "snap_%03d.hdf5" % (i)
        try:
            data = h5py.File(directory+filename, "r")
        except:
            print("could not open "+directory+filename)
            sys.exit(1)
        pos = np.array(data["PartType1"]["Coordinates"], dtype = FloatType) / HubbleParam / 1000.
        posstars = np.array(data["PartType4"]["Coordinates"], dtype = FloatType) / HubbleParam / 1000.
    
        if(pos.shape[0] > 32**3):
            i_select = np.random.uniform(low=0.0, high=pos.shape[0], size=32**3).astype(np.int)
        else:
            i_select = np.arange(pos.shape[0])
        ## Try to show the different gas phases and their spatial distribution
        z = 1./ a - 1
        if i == 0:
            continue
        ax_col = np.int((i-1)%2)
        ax_row = np.int((i-1)/2)
        print('col %d, row %d'%(ax_col, ax_row))
        print(ax[ax_col][ax_row])
        ax[ax_row][ax_col].text(5.5,10.0,'redshift %.1f'%z)
        ax[ax_row][ax_col].scatter(pos[i_select, 0], pos[i_select, 1], marker='.', s=0.01, alpha=0.5, rasterized=True)
        ax[ax_row][ax_col].scatter(posstars[:, 0], posstars[:, 1], marker='*',c='r', s=0.01, alpha=1.0, rasterized=True)
        ax[ax_row][ax_col].set_xlim([0,Boxsize/HubbleParam])
        ax[ax_row][ax_col].set_ylim([0,Boxsize/HubbleParam])
    
    ax[0][0].set_ylabel('[Mpc]')
    ax[1][0].set_xlabel('[Mpc]')
    ax[1][0].set_ylabel('[Mpc]')
    ax[1][1].set_xlabel('[Mpc]')
    fig.savefig(simulation_directory+'/plots/dark_matter_and_stars_evolution.pdf', dpi=300)
    
    
    ## figure -- gas density
    fig, ax = plt.subplots( 2, 2, figsize=np.array([6.9,6.9]), sharex=True, sharey=True )
    fig.subplots_adjust(left = 0.09, bottom = 0.09,right = 0.98, top = 0.98, hspace=0.0, wspace=0.0)
    
    for i, a in enumerate(scalefactor_to_probe):
        directory = simulation_directory+"/output/"
        filename = "fof_subhalo_tab_%03d.hdf5" % (i)
        try:
            data = h5py.File(directory+filename, "r")
        except:
            print("could not open "+directory+filename)
            sys.exit(1)
    
        """ get simulation data """
        ## simulation data
        GrpPos = np.array(data["Group"]["GroupPos"], dtype=FloatType) / HubbleParam / 1000.
        
        delta_mass_stars[i]  = np.interp(a, scalefactor, cum_mass_stars)
        delta_mass_stars[i] -= np.interp(a, scalefactor_ref, cum_mass_stars_ref)
        avg_mass_stars[i]  = 0.5 * np.interp(a, scalefactor, cum_mass_stars)
        avg_mass_stars[i] += 0.5 * np.interp(a, scalefactor_ref, cum_mass_stars_ref)
    
        filename = "snap_%03d.hdf5" % (i)
        try:
            data = h5py.File(directory+filename, "r")
        except:
            print("could not open "+directory+filename)
            sys.exit(1)
        VoronoiPos = np.array(data["PartType0"]["Coordinates"], dtype = FloatType) / HubbleParam / 1000.
        Density = np.array(data["PartType0"]["Density"], dtype = FloatType)
    
        if(pos.shape[0] > 32**3):
            i_select = np.random.uniform(low=0.0, high=pos.shape[0], size=32**3).astype(np.int)
        else:
            i_select = np.arange(pos.shape[0])
        ## Try to show the different gas phases and their spatial distribution
        z = 1./ a - 1
        if i == 0:
            continue
        ax_col = np.int((i-1)%2)
        ax_row = np.int((i-1)/2)

        Nplot = 512
        from scipy import spatial # needed for KDTree that we use for nearest neighbour search and Voronoi mesh
        Edges1d = np.linspace(0., Boxsize/HubbleParam, Nplot+1, endpoint=True, dtype=FloatType)
        Grid1d = 0.5 * (Edges1d[1:] + Edges1d[:-1])
        xx, yy = np.meshgrid(Grid1d, Grid1d)        

        Grid2D = np.array( [xx.reshape(Nplot**2), yy.reshape(Nplot**2), np.ones(Nplot**2)*GrpPos[0,2]] ).T
        dist, cells = spatial.KDTree( VoronoiPos ).query( Grid2D, k=1 )
        gas_density_slice = Density[cells].reshape((Nplot,Nplot))
        
        ax[ax_row][ax_col].pcolormesh( Edges1d, Edges1d, np.log10(gas_density_slice), rasterized=True, cmap=plt.get_cmap('plasma') )        
        ax[ax_row][ax_col].text(5.5,10.0,'redshift %.1f'%z, color='w')
        ax[ax_row][ax_col].set_xlim([0,Boxsize/HubbleParam])
        ax[ax_row][ax_col].set_ylim([0,Boxsize/HubbleParam])
    
    ax[0][0].set_ylabel('[Mpc]')
    ax[1][0].set_xlabel('[Mpc]')
    ax[1][0].set_ylabel('[Mpc]')
    ax[1][1].set_xlabel('[Mpc]')
    fig.savefig(simulation_directory+'/plots/gas_density_evolution.pdf', dpi=300)

    """ figure of stellar mass density """
    fig=plt.figure()
    ax = plt.axes([0.15,0.15,0.8,0.8])
    ax.plot(redshift_plus_one, cum_mass_stars, label="new")
    ax.plot(redshift_plus_one_ref, cum_mass_stars_ref, label="ref")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim([1,11])
    ax.set_xticks([1,2,3,5,9])
    ax.set_xticklabels(["0", "1", "2", "4", "8"])
    ax.set_xlabel(r"redshift")
    ax.set_ylabel(r"Stellar mass denstiy [M$_\odot$ Mpc$^{-3}$]")
    fig.savefig(simulation_directory+"/plots/Mstar_redshift.pdf")

if CompareToReferenceRun:
    """ check if deviations within tolerance """
    abs_tolerance = 5e7 ## depends on box size...
    
    for i, a in enumerate(scalefactor_to_probe):
        delta_mass_stars[i]  = np.interp(a, scalefactor, cum_mass_stars)
        delta_mass_stars[i] -= np.interp(a, scalefactor_ref, cum_mass_stars_ref)
        avg_mass_stars[i]  = 0.5 * np.interp(a, scalefactor, cum_mass_stars)
        avg_mass_stars[i] += 0.5 * np.interp(a, scalefactor_ref, cum_mass_stars_ref)
    
        if np.max( np.abs(delta_mass_stars ) ) > abs_tolerance:
            print("Error: stellar mass deviaties from reference more than %g (Msun/Mpc^3): "% abs_tolerance )
            print("redshifts:")
            print(1./scalefactor_to_probe - 1)
            print("relative deviation:")
            print(delta_mass_stars / avg_mass_stars)
            print("absolute deviation (1e6 Msun/Mpc^3):")
            print(delta_mass_stars/1.0e6)
            sys.exit(1)

""" if everything is ok """
sys.exit(0)
