""" @package examples/cosmo_zoom_3d/create.py
Code that creates ics and outuput list for a cosmological zoom 
simulaiton;

created by Rainer Weinberger, last modified 12.08.2019
"""

""" load libraries """
import sys    # system calls
import numpy as np    # scientific computing package
import h5py    # hdf5 format
import os  # operating system interface
from subprocess import call # execute bash commands

"""
  The usual procedure of producing zoom ICs is to run a uniform
  box fist, select the halo in there and trace back the particles
  to the initlal condition to define a Lagrangian region there
  which will be refined in the zoom ICs. This is done in this 
  script if the flag 'runParent' is active. If it is inactive,
  the zoom simulation is produced right away from a previously
  selected meaningful region.
  
  Note that in Arepo, the high-resolution region must not 
  include a (periodic) boundary at any time of the simulation.
"""
runParent=False


simulation_directory = str(sys.argv[1])
print("create.py " + simulation_directory)


""" output times """
outputTimes = np.array([0.0197,0.2,0.25,0.33,0.5,0.66,1], dtype=np.float64)
ones = np.ones(outputTimes.shape, dtype=np.int)
data = np.array([outputTimes, ones]).T
np.savetxt(simulation_directory+"/output_list.txt",data, fmt="%g %1.f" )


""" IC generating code MUSIC """
# check if music code is there already, otherwise get code from repository and compile IC generating code
if os.path.exists( "{0}/../../../music".format(sys.path[0]) ):
  call(['cp','-r',"{0}/../../../music".format(sys.path[0]),simulation_directory+'/music'])
else:
  print('CREATE: did not fine dir {0}/../../music, trying to clone from bitbucket.'.format(sys.path[0]))
  status = call(['hg', 'clone', 'https://bitbucket.org/ohahn/music', simulation_directory+'/music'])
  if status != 0:
      print('CREATE: ERROR: hg clone failed!')
      sys.exit(status)
cwd = os.getcwd()
os.chdir(simulation_directory+'/music/')
status = call(['make','-j'])
if status != 0:
    print('CREATE: ERROR: make failed!')
    sys.exit(status)
os.chdir(cwd)

if runParent:
    """ parent simulation """
    # execute ic generating code for parent box ICs (change working directory to do this)
    # unfortunately, many working directory changes are necessary here!
    os.chdir(simulation_directory+'/music/')
    status = call(['./MUSIC','../param_music_parent.txt'])
    if status != 0:
        print('CREATE: ERROR: IC creation failed!')
        sys.exit(status)
    os.chdir(cwd)
    
    ConfigFile = simulation_directory+'/Config_parent.sh'
    BuildDir = simulation_directory+'/build_parent'
    ExecFile = simulation_directory+'/Arepo_parent'
    
    ## compile Arepo for cosmological volume
    status = call(['make','-j', 'CONFIG='+ConfigFile, 'BUILD_DIR='+BuildDir, 'EXEC='+ExecFile])
    if status != 0:
        print('CREATE: ERROR: compilation failed!')
        sys.exit(status)
    
    ## run simulation
    os.chdir(simulation_directory)

    ExecFile = './Arepo_parent'
    ParamFile = './param_parent.txt'
    status = call(['mpiexec', '-np', '1', ExecFile, ParamFile])
    if status != 0:
        print('CREATE: ERROR: execution failed!')
        sys.exit(status)
    os.chdir(cwd)
        
    ## select halo and set limits
    sn0 = h5py.File(simulation_directory+'/output_parent/snap_000.hdf5')
    sn6 = h5py.File(simulation_directory+'/output_parent/snap_006.hdf5')
    sub6 = h5py.File(simulation_directory+'/output_parent/fof_subhalo_tab_006.hdf5')
    
    M200c = np.array(sub6['Group']['Group_M_Crit200'], dtype=np.float64)
    R200c = np.array(sub6['Group']['Group_R_Crit200'], dtype=np.float64)
    GrpPos = np.array(sub6['Group']['GroupPos'], dtype=np.float64)
    
    h = 0.6774
    massInterval = [0.5e14/1.0e10*h, 2.0e14/1.0e10*h]
    i_mass_select = np.where( (M200c>massInterval[0]) & (M200c<massInterval[1]) )[0]
    print("Number of halos after mass selection: %d"%len(i_mass_select))
    
    # ISOLATION CRITERION
    minDist = 15000.0*h
    minMassDist = 0.5*massInterval[0]
    i_perturber = np.where( (M200c>minMassDist))[0]

    i_isolated = []
    for i_candidate in i_mass_select:
        dist2 = ( GrpPos[i_candidate,0] - GrpPos[i_perturber, 0] )**2 + \
        ( GrpPos[i_candidate,1] - GrpPos[i_perturber, 1] )**2 + \
        ( GrpPos[i_candidate,2] - GrpPos[i_perturber, 2] )**2
        i_ngb = np.where(dist2 < minDist*minDist)[0]
        if len(i_ngb) < 2:   ## always should contain itself
            i_isolated.append(i_candidate)
            
    # SELECTION OF THE MOST CENTRAL GALAXY (CLOSEST TO CENTER OF BOX)
    BoxHalf = 50000.0
    
    dists = []
    for i_candidate in i_isolated:
        dist2 = ( GrpPos[i_candidate, 0] - BoxHalf)**2 + \
        ( GrpPos[i_candidate, 1] - BoxHalf)**2 + \
        ( GrpPos[i_candidate, 2] - BoxHalf)**2
        dists.append(np.sqrt(dist2))
    
    dists = np.array(dists)
    i_select = np.where( dists == np.min(dists))[0][0]
    i_halo = i_isolated[i_select]
    pos = GrpPos[i_halo, :]
    print('Number of halos after isolation criterion selection: %d'%len(i_isolated))
    print('Selected halo %d with position %g %g %g'%(i_halo, pos[0], pos[1], pos[2]))
                
    # SELECTION OF PARTICLES IN HIGH-RES REGION
    pid0 = np.array(sn0['PartType1']['ParticleIDs'], dtype=np.int64)
    pid6 = np.array(sn6['PartType1']['ParticleIDs'], dtype=np.int64)
    ppos0 = np.array(sn0['PartType1']['Coordinates'], dtype=np.float64)
    ppos6 = np.array(sn6['PartType1']['Coordinates'], dtype=np.float64)
    
    
    dist = (ppos6[:,0]-pos[0])**2 + \
    (ppos6[:,1]-pos[1])**2 + \
    (ppos6[:,2]-pos[2])**2 
    dist = np.sqrt(dist)
    i_particles = np.where(dist < 3.0*R200c[i_halo])[0]
    print('selected %d particles to trace back Lagrangian region'%(len(i_particles) ) )
    id_particle_flagged = pid6[i_particles]
    
    # SEARCH FOR THESE PARTCLES IN INITIAL CONDITION TO FIND BOUNDS
    # return_indices works only with numpy 1.15.0 and above
    ids, i_select1, dummy = np.intersect1d(pid0, id_particle_flagged, return_indices=True)
    print('number of matched ids: %d'%len(ids))
    
    #high-res limit: 
    xmin, xmax = np.min(ppos6[i_particles,0]), np.max(ppos6[i_particles,0])
    ymin, ymax = np.min(ppos6[i_particles,1]), np.max(ppos6[i_particles,1])
    zmin, zmax = np.min(ppos6[i_particles,2]), np.max(ppos6[i_particles,2])
    print('halo region at z=0')
    print(xmin, xmax, xmax-xmin, 0.5*(xmax+xmin))
    print(ymin, ymax, ymax-ymin, 0.5*(ymax+ymin))
    print(zmin, zmax, zmax-zmin, 0.5*(zmax+zmin))
    
    xmin, xmax = np.min(ppos0[i_select1,0]), np.max(ppos0[i_select1,0])
    ymin, ymax = np.min(ppos0[i_select1,1]), np.max(ppos0[i_select1,1])
    zmin, zmax = np.min(ppos0[i_select1,2]), np.max(ppos0[i_select1,2])
    print('lagrangian region at ICs')
    print(xmin, xmax, xmax-xmin, 0.5*(xmax+xmin))
    print(ymin, ymax, ymax-ymin, 0.5*(ymax+ymin))
    print(zmin, zmax, zmax-zmin, 0.5*(zmax+zmin))
    
    center = [0.5 * (xmax+xmin) / 100000.0, 0.5 * (ymax+ymin) / 100000.0, 0.5 * (zmax+zmin) / 100000.0]
    box = [( (xmax-xmin) + 1000. )/100000.,( (ymax-ymin) + 1000. )/100000.,( (zmax-zmin) + 1000.)/100000.]
    
    ## replace the high-res extent in param_music.txt
    paramFile = open(simulation_directory+'/param_music_parent.txt', 'r')
    paramOptions = []
    for line in paramFile:
        if line[:10] == 'ref_center':
            paramOptions.append('ref_center         = %.5f, %.5f, %.5f \n'%(center[0], center[1], center[2]) )
        elif line[:10] == 'ref_extent':
            paramOptions.append('ref_extent         = %.5f, %.5f, %.5f \n'%(box[0], box[1], box[2]) )
        elif line[:8] == 'levelmax':
        	paramOptions.append('levelmax           = 9 \n') ## increase resolution in zoom region
        elif line[:8] == 'filename':
            paramOptions.append('filename           = ../ics \n') ## change file name of zoom ics
        else:
            paramOptions.append(line)
    paramFile.close()
    paramFile = open(simulation_directory+'/param_music.txt', 'w')
    for line in paramOptions:
        paramFile.write(line)
    paramFile.close()


# execute ic generating code for zoom ICs (change working directory to do this)
os.chdir(simulation_directory+'/music/')
status = call(['./MUSIC', '../param_music.txt'])
if status != 0:
    print('CREATE: ERROR: execution failed!')
    sys.exit(status)
os.chdir(cwd)


""" normal exit """
sys.exit(0) 
