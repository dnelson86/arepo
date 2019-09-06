#!/bin/bash            # this line only there to enable syntax highlighting in this file

##################################################
#  Enable/Disable compile-time options as needed #
##################################################

#--------------------------------------- Basic operation mode of code; default: 3d with 6 particle types; type 0: gas >0: only gravitationally interacting
#NTYPES=6                      # number of particle types
#TWODIMS                       # 2d simulation
#ONEDIMS                       # 1d simulation
#ONEDIMS_SPHERICAL             # 1d spherically symmetric simulation

#--------------------------------------- Computational box and boundaries; default: cubic, periodic box
#LONG_X=10.0                   # stretch x extent of box by given factor
#LONG_Y=2.0                    # stretch y extent of box by given factor
#LONG_Z=10.0                   # stretch z extent of box by given factor
#REFLECTIVE_X=1 #=2            # X Boundary; 1: Reflective, 2: Inflow/Outflow; not active: periodic
#REFLECTIVE_Y=1 #=2            # Y Boundary; 1: Reflective, 2: Inflow/Outflow; not active: periodic
#REFLECTIVE_Z=1 #=2            # Z Boundary; 1: Reflective, 2: Inflow/Outflow; not active: periodic

#--------------------------------------- Hydrodynamics; default: GAMMA=5/3 ideal hydrodynamics
#NOHYDRO                       # No hydrodynamics calculation
#GAMMA=1.4                     # Adiabatic index of gas; 5/3 if not set
#ISOTHERM_EQS                  # Isothermal gas
#PASSIVE_SCALARS=3             # number of passive scalar fields advected with fluid (default: 0)
#NO_SCALAR_GRADIENTS           # disables time and spatial extrapolation for passive scalar fields (use only if you know why you're doing this)

#--------------------------------------- Magnetohydrodynamics
#MHD                           # Master switch for magnetohydrodynamics
#MHD_POWELL                    # Powell div(B) cleaning scheme for magnetohydrodynamics
#MHD_POWELL_LIMIT_TIMESTEP     # Timestep constraint due to Powell cleaning scheme
#MHD_SEEDFIELD                 # Uniform magnetic seed field of specified orientation and strength set up after reading in IC

#--------------------------------------- Riemann solver; default: exact Riemann solver
#RIEMANN_HLLC                  # HLLC approximate Riemann solver
#RIEMANN_HLLD                  # HLLD approximate Riemann solver (required to use for MHD)

#--------------------------------------- Mesh motion and regularization; default: moving mesh
#VORONOI_STATIC_MESH           # static mesh
#VORONOI_STATIC_MESH_DO_DOMAIN_DECOMPOSITION  # for VORONOI_STATIC_MESH force domain decomposition if there exist non-gas particles
#REGULARIZE_MESH_CM_DRIFT      # Mesh regularization; Move mesh generating point towards center of mass to make cells rounder.
#REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED  # Limit mesh regularization speed by local sound speed
#REGULARIZE_MESH_FACE_ANGLE    # Use maximum face angle as roundness criterion in mesh regularization

#--------------------------------------- Refinement and derefinement; default: no refinement/derefinement; criterion: target mass
#REFINEMENT_SPLIT_CELLS        # Refinement
#REFINEMENT_MERGE_CELLS        # Derefinement
#REFINEMENT_VOLUME_LIMIT       # Limit the volume of cells and the maximum volume difference between neighboring cels
#JEANS_REFINEMENT              # Refinement criterion to ensure Jeans stability of cells
#REFINEMENT_HIGH_RES_GAS       # Refinement criterion for high-resolution region in zoom simulation
#NODEREFINE_BACKGROUND_GRID    # Do not de-refine low-res gas cells in zoom simulations
#OPTIMIZE_MESH_MEMORY_FOR_REFINEMENT  # deletes the mesh structures not needed for refinement/derefinemet to lower the peak memory consumption

#--------------------------------------- non-standard phyiscs
#COOLING                       # Simple primordial cooling
#ENFORCE_JEANS_STABILITY_OF_CELLS  # this imposes an adaptive floor for the temperature
#USE_SFR                       # Star formation model, turning dense gas into collisionless partices
#SFR_KEEP_CELLS                # Do not distroy cell out of which a star has formed

#--------------------------------------- Gravity treatment; default: no gravity
#SELFGRAVITY                   # gravitational intraction between simulation particles/cells
#HIERARCHICAL_GRAVITY          # use hierarchical splitting of the time integration of the gravity
#CELL_CENTER_GRAVITY           # uses geometric centers to calculate gravity of cells, only possible with HIERARCHICAL_GRAVITY
#NO_GAS_SELFGRAVITY            # switch off gas self-gravity in tree
#GRAVITY_NOT_PERIODIC          # gravity is not treated periodically
#ALLOW_DIRECT_SUMMATION        # Performed direct summation instead of tree-based gravity if number of active particles < DIRECT_SUMMATION_THRESHOLD (= 3000 unless specified differently here)
#DIRECT_SUMMATION_THRESHOLD=1000  # Overrides maximum number of active particles for which direct summation is performed instead of tree based calculation
#EXACT_GRAVITY_FOR_PARTICLE_TYPE=4  #N-squared fashion gravity for a small number of particles of the given type
#EVALPOTENTIAL                 # computes gravitational potential

#--------------------------------------- TreePM Options; default: no Particle-Mesh
#PMGRID=512                    # Enables particle mesh; number of cells used for grid in each dimension
#ASMTH=1.25                    # This can be used to override the value assumed for the scale that defines the long-range/short-range force-split in the TreePM algorithm. The default value is 1.25, in mesh-cells.
#RCUT=6.0                      # This can be used to override the maximum radius in which the short-range tree-force is evaluated (in case the TreePM algorithm is used). The default value is 4.5, given in mesh-cells.
#PM_ZOOM_OPTIMIZED             # Particle-mesh calculation that is optimized for cosmological zoom simulations; disable for cosmological volume simulations.
#PLACEHIGHRESREGION=2          # Places a second, high-resolution PM grid; number encodes high-res particle types from which the extent of this second PM grid is calculated.
#ENLARGEREGION=1.1             # Factor by which high res region is increased with respect to max-min position of high-res particles
#GRIDBOOST=2                   # Factor by which PMGRID is increased in non-periodic (or high res) PM calculation (if not set, code uses 2)
#FFT_COLUMN_BASED              # Use column-based FFT; slightly slower, but necessary to achieve good load-balancing if number of tasks larger than PMGRID (usually the case for very large runs)

#--------------------------------------- Gravity softening
#NSOFTTYPES=4                  # Number of different softening values to which particle types can be mapped.
#MULTIPLE_NODE_SOFTENING       # If a tree node is to be used which is softened, this is done with the softenings of its different mass components
#INDIVIDUAL_GRAVITY_SOFTENING=2+4  # bitmask with particle types where the softenig type should be chosen with that of parttype 1 as a reference type
#ADAPTIVE_HYDRO_SOFTENING      # Adaptive softening of gas cells depending on their size
#NSOFTTYPES_HYDRO=64           # Overrides number of discrete softening values for gas cellls when ADAPTIVE_HYDRO_SOFTENING (default is 64)

#--------------------------------------- External gravity; default: no external potential
#EXTERNALGRAVITY               # master switch for external potential
#EXTERNALGY=0.0                # constant external gravity in y direction

#--------------------------------------- Static NFW Potential
#STATICNFW                     # static gravitational Navarro-Frenk-White (NFW) potential
#NFW_C=12                      # concentration parameter of NFW potential
#NFW_M200=100.0                # mass causing the NFW potential
#NFW_Eps=0.01                  # softening of NFW potential
#NFW_DARKFRACTION=0.87         # fraction in dark matter in NFW potential

#--------------------------------------- Static Isothermal Sphere Potential
#STATICISO                     # static gravitational isothermal sphere potential
#ISO_M200=100.0                # mass causing the isothermal sphere potential
#ISO_R200=160.0                # radius of the isothermal sphere potential
#ISO_Eps=0.1                   # softening of isothermal sphere potential
#ISO_FRACTION=0.9              # fraction in dark matter in isothermal sphere potential

#--------------------------------------- Static Hernquist Potential
#STATICHQ                      # static gravitational Hernquist potential
#HQ_M200=186.015773            # mass causing the Hernquist potential
#HQ_C=10.0                     # concentration parameter of Hernquist potential
#HQ_DARKFRACTION=0.9           # fraction in dark matter in Hernquist potential

#--------------------------------------- Time integration options
#FORCE_EQUAL_TIMESTEPS         # variable but global timestep
#TREE_BASED_TIMESTEPS          # non-local timestep criterion (take 'signal speed' into account)
#PM_TIMESTEP_BASED_ON_TYPES=2+4  # particle types that should be considered in setting the PM timestep
#NO_PMFORCE_IN_SHORT_RANGE_TIMESTEP  # PM force is not included in short-range timestep criterion
#ENLARGE_DYNAMIC_RANGE_IN_TIME # This extends the dynamic range of the integer timeline from 32 to 64 bit

#--------------------------------------- MPI
#IMPOSE_PINNING                # Enforce pinning of MPI tasks to cores if MPI does not do it
#IMPOSE_PINNING_OVERRIDE_MODE  # Override MPI pinning, if present

#--------------------------------------- Single/Double Precision
#DOUBLEPRECISION=1             # Mode of double precision: not defined: single; 1: full double precision 2: mixed, 3: mixed, fewer single precisions; unless short of memory, use 1.
#DOUBLEPRECISION_FFTW          # FFTW calculation in double precision
#OUTPUT_IN_DOUBLEPRECISION     # snapshot files will be written in double precision
#INPUT_IN_DOUBLEPRECISION      # initial conditions are in double precision
#OUTPUT_COORDINATES_IN_DOUBLEPRECISION  # will always output coordinates in double precision
#NGB_TREE_DOUBLEPRECISION      # if this is enabled, double precision is used for the neighbor node extension

#--------------------------------------- On the fly FOF groupfinder
#FOF                           # enable FoF output
#FOF_PRIMARY_LINK_TYPES=2      # 2^type for the primary dark matter type
#FOF_SECONDARY_LINK_TYPES=1+16+32  # 2^type for the types linked to nearest primaries
#FOF_SECONDARY_LINK_TARGET_TYPES=  # should normally be set to a list of all dark matter types (in zoom runs), if not set defaults to FOF_PRIMARY_LINK_TYPES
#FOF_GROUP_MIN_LEN=32          # Minimum number of particles in one group (default is 32)
#FOF_LINKLENGTH=0.16           # Linkinglength for FoF (default=0.2)
#FOF_STOREIDS                  # store IDs in group/subfind catalogue, do not order particles in snapshot files by group order

#--------------------------------------- Subfind
#SUBFIND                       # enables substructure finder
#SAVE_HSML_IN_SNAPSHOT         # stores hsml, density, and velocity dispersion values in the snapshot files
#SUBFIND_CALC_MORE             # calculates also the velocity dispersion in the local density estimate (this is automatically enabled by several other options, e.g. SAVE_HSML_IN_SNAPSHOT)
#SUBFIND_EXTENDED_PROPERTIES   # adds calculation of further quantities related to angular momentum in different components

#--------------------------------------- Things for special behaviour
#RUNNING_SAFETY_FILE           # if file './running' exists, do not start the run
#MULTIPLE_RESTARTS             # Keep restart files instead of just two copies
#EXTENDED_GHOST_SEARCH         # This extends the ghost search to the full 3x3 domain instead of the principal domain
#DOUBLE_STENCIL                # this will ensure that the boundary region of the local mesh is deep enough to have a valid double stencil for all local cells
#TETRA_INDEX_IN_FACE           # adds an index to each entry of VF[] and DC[] to one of the tetrahedra that share this edge
#NOSTOP_WHEN_BELOW_MINTIMESTEP # Simulation does not terminate when timestep drops below minimum timestep
#TIMESTEP_OUTPUT_LIMIT         # Limit timesteps to write snaps on time for output lists with huge range
#ALLOWEXTRAPARAMS              # Tolerate extra parameters that are not used
#FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES  # this can be used to load SPH ICs that contain identical particle coordinates
#RECOMPUTE_POTENTIAL_IN_SNAPSHOT  # needed for postprocess option 18 that can be used to calculate potential values for a snapshot
#ACTIVATE_MINIMUM_OPENING_ANGLE  # this does not open tree nodes under the relative opening criterion any more if their opening angle has dropped below a minimum angle
#USE_DIRECT_IO_FOR_RESTARTS    # Try to use O_DIRECT for low-level read/write operations of restart files to circumvent the linux kernel page caching
#HUGEPAGES                     # use huge pages for memory allocation, through hugetlbfs library
#DETAILEDTIMINGS               # creates individual timings entries for primary/secondary kernels to diagnose work-load balancing
#BITS_PER_DIMENSION=42         # Peano-Hilbert order

#--------------------------------------- input options
#COMBINETYPES                  # reads in the IC file types 4+5 as type 3
#LOAD_TYPES=1+2+4+16+32        # load only specific types sum(2^type)
#READ_COORDINATES_IN_DOUBLE    # read coordinates in double precision
#LONGIDS                       # Store IDs in type unsigned long long instead of unsigned int
#OFFSET_FOR_NON_CONTIGUOUS_IDS # Determines offset of IDs on startup instead of using fixed offset.
#GENERATE_GAS_IN_ICS           # Generates gas from dark matter only ICs.
#SPLIT_PARTICLE_TYPE=4+8       # Overrides splitting particle type 1 in GENERATE_GAS_IN_ICS use sum(2^type)
#SHIFT_BY_HALF_BOX             # Shift all positions by half a box size after reading in
#NTYPES_ICS=6                  # number of particle types in ICs, if not NTYPES (only works for 6, and non-HDF5 ICs!)
#READ_MASS_AS_DENSITY_IN_INPUT # Reads the mass field in the IC as density

#--------------------------------------- special input options
#IDS_OFFSET=1                  # Override offset for gas particles if created from DM
#READ_DM_AS_GAS                # reads in dark matter particles as gas cells
#TILE_ICS                      # tile ICs by TileICsFactor in each dimension

#--------------------------------------- output fields in snapshots--Default output filds are: position, velocity, ID, mass, spec. internal energy (gas), density(gas)
#OUTPUT_TASK                   # output MPI task
#OUTPUT_TIMEBIN_HYDRO          # output hydrodynamics time-bin
#OUTPUT_PRESSURE_GRADIENT      # output pressure gradient
#OUTPUT_DENSITY_GRADIENT       # output density gradient
#OUTPUT_VELOCITY_GRADIENT      # output velocity gradient
#OUTPUT_BFIELD_GRADIENT        # output magnetic field gradient
#OUTPUT_MESH_FACE_ANGLE        # output max. face angle of cells
#OUTPUT_VERTEX_VELOCITY        # output velocity of cell
#OUTPUT_VERTEX_VELOCITY_DIVERGENCE # output deivergence of cell velocity
#OUTPUT_VOLUME                 # output volume of cells; note that this can always be computat as both, density and mass of cells are by default in output
#OUTPUT_CENTER_OF_MASS         # output center of mass of cells (position is mesh-generating point)
#OUTPUT_SURFACE_AREA           # output surface area of cells
#OUTPUT_PRESSURE               # output pressure of gas
#OUTPUTPOTENTIAL               # output potential at particle position
#OUTPUTACCELERATION            # output gravitational acceleration
#OUTPUTTIMESTEP                # output timestep of particle
#OUTPUT_SOFTENINGS             # output particle softenings
#OUTPUTGRAVINTERACTIONS        # output gravitatational interactions (from the tree) of particles
#OUTPUTCOOLRATE                # outputs cooling rate, and conduction rate if enabled
#OUTPUT_DIVVEL                 # output  velocity divergence
#OUTPUT_CURLVEL                # output  velocity curl
#OUTPUT_COOLHEAT               # output actual energy loss/gain in cooling/heating routine
#OUTPUT_VORTICITY              # output vorticity of gas
#OUTPUT_CSND                   # output sound speed. This one is only used for tree-based timesteps! Calculate from hydro quantities in postprocessing if required for science applications.

#--------------------------------------- output options
#PROCESS_TIMES_OF_OUTPUTLIST   # goes through times of output list prior to starting the simulaiton to ensure that outputs are written as close to the desired time as possible (as opposed to at next possible time if this flag is not active)
#REDUCE_FLUSH                  # only flush output to log-files in predefined intervals
#OUTPUT_EVERY_STEP             # Create snapshot on every (global) synchronization point, independent of parameters choosen or output list.
#OUTPUT_CPU_CSV                # output of a cpu.csv file on top of cpu.txt
#HAVE_HDF5                     # needed when HDF5 I/O support is desired (recommended)
#HDF5_FILTERS                  # activate snapshot compression and checksum for HDF5 output
#OUTPUT_XDMF                   # writes an .xmf file for each snapshot, which can be read by visit (with the hdf5 snapshot)

#--------------------------------------- Testing and Debugging options
#DEBUG                         # enables core-dumps
#VERBOSE                       # reports readjustments of buffer sizes

#--------------------------------------- Mesh-relaxing or mesh-adding (this will not carry out a simulation)
#MESHRELAX                     # this keeps the mass constant and only regularizes the mesh
#ADDBACKGROUNDGRID=16          # Re-grid hydrodynamics quantities on a Oct-tree AMR grid. This does not perform a simulation.
