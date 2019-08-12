""" @package ./examples/cosmo_zoom_3d/check.py
Code that checks results of gravity only cosmological zoom simulation

created by Rainer Weinberger, last modified 05.03.2019
"""


""" load libraries """
import sys    ## load sys; needed for exit codes
import numpy as np    ## load numpy
import h5py    ## load h5py; needed to read snapshots
import os      # file specific calls
import matplotlib.pyplot as plt    # needs to be active for plotting!
plt.rcParams['text.usetex'] = True

createReferenceSolution = False
makeplots = True
if len(sys.argv) > 2:
  if sys.argv[2] == "True":
    makeplots = True
  else:
    makeplots = False

simulation_directory = str(sys.argv[1])
print("cosmo_zoom_3d: checking simulation output in directory " + simulation_directory) 

FloatType = np.float64  # double precision: np.float64, for single use np.float32
IntType = np.int32

Redshifts = np.array([50,4,3,2,1,0.5,0])
HubbleParam = 0.6774
UnitMass = 1.0e10

status = 0

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
    zpos = 62.
    if i_file != 0:
        SubhaloMass = np.array(data["Subhalo"]["SubhaloMass"], dtype=FloatType) * UnitMass / HubbleParam
        SubhaloPos = np.array(data["Subhalo"]["SubhaloPos"], dtype=FloatType)/ HubbleParam / 1000.
        SubhaloRad = np.array(data["Subhalo"]["SubhaloHalfmassRad"], dtype=FloatType)/ HubbleParam / 1000.
        GrpPos = np.array(data["Group"]["GroupPos"], dtype=FloatType)/ HubbleParam / 1000.
        GrpR200c = np.array(data["Group"]["Group_R_Crit200"], dtype=FloatType)/ HubbleParam / 1000.
        
        zpos = GrpPos[0,2]
        ## select subhalos within main halo only (2*R200c)
        dist  = (SubhaloPos[:,0] - GrpPos[0,0])**2
        dist += (SubhaloPos[:,1] - GrpPos[0,1])**2
        dist += (SubhaloPos[:,2] - GrpPos[0,2])**2
        dist = np.sqrt(dist)
        
        i_select, = np.where(dist < 2.0 * GrpR200c[0])
        
        n_subs = len(SubhaloMass[i_select])
        
        SubhaloMass = np.sort(SubhaloMass[i_select])[::-1]
        
        if createReferenceSolution:
            # save reference masses
            np.savetxt(simulation_directory+"/Masses_z%.1g.txt"% z, SubhaloMass)

    """ optional figure """
    if makeplots:
        filename = "snap_%03d.hdf5" % (i_file)
        try:
            data = h5py.File(directory+filename, "r")
        except:
            print("could not open "+directory+filename)
            sys.exit(1)
        pos = np.array(data["PartType1"]["Coordinates"], dtype=FloatType) / HubbleParam / 1000.
        pos2 = np.array(data["PartType2"]["Coordinates"], dtype=FloatType) / HubbleParam / 1000.
        vel2 = np.array(data["PartType2"]["Velocities"], dtype=FloatType)
        mass2 = np.array(data["PartType2"]["Masses"], dtype=FloatType) / HubbleParam * UnitMass
        id2 = np.array(data["PartType2"]["ParticleIDs"], dtype=IntType)
        
        if i_file != 0:
            cont_dist =  (pos2[:,0] - GrpPos[0,0])**2
            cont_dist += (pos2[:,1] - GrpPos[0,1])**2
            cont_dist += (pos2[:,2] - GrpPos[0,2])**2
            cont_dist = np.sqrt(cont_dist)
            print('minimum distance of contaminating (low-res) particle: %g Mpc'%np.min(cont_dist) )
            i_issue = np.where(cont_dist < GrpR200c[0])[0]
            
            if len(i_issue) > 0:
                filename = "snap_%03d.hdf5" % (0)
                data = h5py.File(directory+filename, "r")
                id2_ic = np.array(data["PartType2"]["ParticleIDs"], dtype=IntType)
                pos2_ic = np.array(data["PartType2"]["Coordinates"], dtype=FloatType)
                
                ids, i_select1, dummy = np.intersect1d(id2_ic, id2[i_issue], return_indices=True)
                print('problem at positions (code units!)')
                print(pos2_ic[i_select1,:])
                print(id2_ic[i_select1])
                exit(1)
    
        fig = plt.figure(figsize=(6.9,6.9))
        ax = plt.axes([0.1,0.1,0.87,0.87])
        
        particle_select = np.where( (pos[:,2]>zpos-2.5) & (pos[:,2]<zpos+2.5) )[0]
        particle2_select = np.where( (pos2[:,2]>zpos-2.5) & (pos2[:,2]<zpos+2.5) )[0]
        ax.scatter(pos[particle_select, 0], pos[particle_select, 1], marker='.', s=0.05, alpha=0.5, rasterized=True)
        ax.scatter(pos2[particle2_select, 0],pos2[particle2_select,1], marker='.', s=5, c='r', alpha=0.5, rasterized=True)
        

        if i_file != 0:
            ax.add_artist(plt.Circle((GrpPos[0,0], GrpPos[0,1]),GrpR200c[0],linestyle='--', color='k', fill=False))
            
            for i in np.arange(SubhaloRad.shape[0]):
                if (dist[i] < 2.0 * GrpR200c[0]):
                    ax.add_artist(plt.Circle((SubhaloPos[i,0], SubhaloPos[i,1]), SubhaloRad[i], color='k', fill=False))
        ax.set_xlim([0,100./HubbleParam])
        ax.set_ylim([0,100./HubbleParam])
        
        if i_file != 0:
            ax.set_xlim( GrpPos[0,0]-1.75, GrpPos[0,0]+1.25 )
            ax.set_ylim( GrpPos[0,1]-1.25, GrpPos[0,1]+1.75 )
        ax.set_xlabel('[Mpc]')
        ax.set_ylabel('[Mpc]')
        
        if i_file != 0:
            CumMassFunction = np.cumsum( np.ones(SubhaloMass.shape[0]))
            bx = plt.axes([0.20,0.74,0.26,0.22])
            bx.plot(SubhaloMass, CumMassFunction)
            bx.set_xscale('log')
            bx.set_yscale('log')
            bx.set_xlim([9e9,3e12])
            bx.set_ylim([0.9,150])
            bx.set_xlabel(r"$M_{h}\ \mathrm{[M_\odot]}$")
            bx.set_ylabel(r"$N(>M)$")
        
        if not os.path.exists( simulation_directory+"/plots" ):
            os.mkdir( simulation_directory+"/plots" )
        fig.savefig(simulation_directory+'/plots/haloStructure_z%.1d.pdf'%z, dpi=300)

    ## comparison to reference run (sorted list of SubhaloMass)
    if i_file > 3:
        SubhaloMass_ref = np.loadtxt(simulation_directory+"/Masses_z%.1g.txt"% z)
        minLen = np.min([len(SubhaloMass), len(SubhaloMass_ref)])
        i_select2 = np.arange(minLen)
        delta = np.array((SubhaloMass[i_select2]-SubhaloMass_ref[i_select2]) / np.sqrt(0.5 * (SubhaloMass[i_select2]+SubhaloMass_ref[i_select2]) ), dtype=np.float64)
        delta = np.abs(delta)
        tolerance = 5e4  ## corresponds to 5e11 Msun for 1e14 Msun halo and 5e10 for 2e12 Msun subhalo
        ## typical average delta values (empirically determined by varying 
        ## the number of MPI ranks) is around 3e4
        
        if np.abs(np.average(delta)) > tolerance:
            status = 1
            print("ERROR: z=%g difference in subhalo masses exceeding limits!" % z)
            print("average mass error/sqrt(mass) (=delta)")
            print(np.average(delta))
            print("tolerance: %g" % tolerance)
            print("Individual error values:")
            for element in delta:
                print(element)
            
           
 

## if everything is ok: 0 else: 1
sys.exit(status)
