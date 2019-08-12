""" @package ./examples/cosmo_box_gravity_only_3d/check.py
Code that checks results of gravity only structure formation simulation

created by Rainer Weinberger, last modified 25.02.2019
"""

""" load libraries """
import sys    ## load sys; needed for exit codes
import numpy as np    ## load numpy
import h5py    ## load h5py; needed to read snapshots
import os      # file specific calls
import matplotlib.pyplot as plt    # needs to be active for plotting!
plt.rcParams['text.usetex'] = True

simulation_directory = str(sys.argv[1])
print("cosmo_box_gravity_only_3d: checking simulation output in directory " + simulation_directory) 

FloatType = np.float64  # double precision: np.float64, for single use np.float32

Boxsize = 50 ##Mpc/h
HubbleParam = 0.6774
UnitMass = 1.0e10
Volume = Boxsize * Boxsize * Boxsize / HubbleParam / HubbleParam / HubbleParam

Redshifts = [1, 0]
status = 0

CompareAgainstReferenceRun = True ## comparison for small L50m32 box; deactivate this when comparing against self-created ICs
CreateReferenceSolution = False
makeplots = True
if len(sys.argv) > 2:
  if sys.argv[2] == "True":
    makeplots = True
  else:
    makeplots = False

for i_file, z in enumerate(Redshifts):
    """ try to read in snapshot """
    directory = simulation_directory+"/output/"
    filename = "fof_subhalo_tab_%03d.hdf5" % (i_file)
    try:
        data = h5py.File(directory+filename, "r")
    except:
        print("could not open "+directory+filename)
        sys.exit(1)
    
    """ get simulation data """
    ## simulation data
    GrpPos = np.array(data["Group"]["GroupPos"], dtype=FloatType) / HubbleParam / 1000.
    GrpR200c = np.array(data["Group"]["Group_R_Crit200"], dtype=FloatType) / HubbleParam / 1000.
    M200c = np.array(data["Group"]["GroupMass"], dtype=FloatType) * UnitMass
    M200c = np.sort(M200c)[::-1]
    CumMassFunction = np.cumsum(np.ones(M200c.shape) ) / Volume
    
    if CreateReferenceSolution:
        np.savetxt(simulation_directory+"/Masses_L50n32_z%.1d.txt"% z, M200c)
    
    if CompareAgainstReferenceRun:
        ## comparison to reference run (sorted list of M200)
        M200c_ref = np.loadtxt(simulation_directory+"/Masses_L50n32_z%.1d.txt"% z)
        
        minLen = np.min([len(M200c), len(M200c_ref)])
        i_select = np.arange(minLen)
        
        delta = (M200c[i_select]-M200c_ref[i_select]) / M200c_ref[i_select]
        
        ## empirically based tolerances
        tolerance_average = 0.01
        tolerance_std = 0.05
        if np.abs(np.average(delta)) > tolerance_average or np.abs(np.std(delta)) > tolerance_std:
            status = 1
            print("ERROR: z=%g difference in halo masses exceeding limits!" % z)
            print("relative mass error (=delta)")
            print(delta)
            print("average delta (tolerance: %g)" % tolerance_average)
            print(np.average(delta))
            print("stddev delta (tolerance: %g)" % tolerance_std)
            print(np.std(delta))
    
    """ optional figure """
    if makeplots:
        filename = "snap_%03d.hdf5" % (i_file)
        try:
            data = h5py.File(directory+filename, "r")
        except:
            print("could not open "+directory+filename)
            sys.exit(1)
        pos = np.array(data["PartType1"]["Coordinates"], dtype=FloatType) / HubbleParam / 1000.
    
        fig = plt.figure(figsize=(6.9,6.9))
        ax = plt.axes([0.1,0.1,0.87,0.87])
        
        if(pos.shape[0] > 32**3):
            i_select = np.random.uniform(low=0.0, high=pos.shape[0], size=32**3).astype(np.int)
        else:
            i_select = np.arange(pos.shape[0])
        ax.scatter(pos[i_select, 0], pos[i_select, 1], marker='.', s=0.05, alpha=0.5, rasterized=True)
        for i in np.arange(GrpR200c.shape[0]):
            ax.add_artist(plt.Circle((GrpPos[i,0], GrpPos[i,1]), GrpR200c[i], color='k', fill=False))
        
        ax.set_xlim([0,Boxsize/HubbleParam])
        ax.set_ylim([0,Boxsize/HubbleParam])
        ax.set_xlabel('[Mpc]')
        ax.set_ylabel('[Mpc]')
        
        bx = plt.axes([0.70,0.74,0.26,0.22])
        bx.plot(M200c, CumMassFunction)
        bx.set_xscale('log')
        bx.set_yscale('log')
        bx.set_xlim([9e12,5e14])
        bx.set_ylim([9e-7,2e-4])
        bx.set_xlabel(r"$M_{200,c}\ \mathrm{[M_\odot]}$")
        bx.set_ylabel(r"$n(>M)\ \mathrm{[Mpc^{-3}]}$")
        
        if not os.path.exists( simulation_directory+"/plots" ):
          os.mkdir( simulation_directory+"/plots" )
        fig.savefig(simulation_directory+'/plots/largeScaleStructure_z%.1d.pdf'%z, dpi=300)

## if everything is ok: 0 else: 1
sys.exit(status)
