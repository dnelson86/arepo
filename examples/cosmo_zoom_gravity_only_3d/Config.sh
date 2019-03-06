#!/bin/bash        	# this line only there to enable syntax highlighting in this file


###################################################
## examples/cosmo_zoom_gravity_only_3d/Config.sh ##
###################################################


#--------------------------------------- Gravity treatment
SELFGRAVITY                              # gravitational intraction between simulation particles/cells
HIERARCHICAL_GRAVITY                     # use hierarchical splitting of the time integration of the gravity
ALLOW_DIRECT_SUMMATION                   # Performed direct summation instead of tree-based gravity if number of active particles < DIRECT_SUMMATION_THRESHOLD (= 3000 unless specified differently here)
DIRECT_SUMMATION_THRESHOLD=1024          # Overrides maximum number of active particles for which direct summation is performed instead of tree based calculation
EVALPOTENTIAL                            # computes gravitational potential


#--------------------------------------- TreePM Options
PMGRID=256                               # Enables particle mesh; number of cells used for grid in each dimension
RCUT=5.5                                 # This can be used to override the maximum radius in which the short-range tree-force is evaluated (in case the TreePM algorithm is used). The default value is 4.5, given in mesh-cells.
PLACEHIGHRESREGION=2                     # Places a second, high-resolution PM grid; number encodes high-res particle types from which the extent of this second PM grid is calculated.
ENLARGEREGION=1.1                        # Factor by which high res region is increased with respect to max-min position of high-res particles
GRIDBOOST=2                              # Factor by which PMGRID is increased in non-periodic (or high res) PM calculation (if not set, code uses 2)
PM_ZOOM_OPTIMIZED                        # Particle-mesh calculation that is optimized for cosmological zoom simulations; disable for cosmological volume simulations.


#--------------------------------------- Gravity softening
MULTIPLE_NODE_SOFTENING                  # If a tree node is to be used which is softened, this is done with the softenings of its different mass components
INDIVIDUAL_GRAVITY_SOFTENING=4+8+16+32   # bitmask with particle types where the softenig type should be chosen with that of parttype 1 as a reference type


#--------------------------------------- Time integration options
TREE_BASED_TIMESTEPS                     # non-local timestep criterion (take 'signal speed' into account)


#--------------------------------------- Single/Double Precision
DOUBLEPRECISION=1                        # Mode of double precision: not defined: single; 1: full double precision 2: mixed, 3: mixed, fewer single precisions; unless short of memory, use 1.
DOUBLEPRECISION_FFTW                     # FFTW calculation in double precision


#--------------------------------------- On the fly FOF groupfinder
FOF                                      # enable FoF output
FOF_PRIMARY_LINK_TYPES=2                 # 2^type for the primary dark matter type
FOF_SECONDARY_LINK_TYPES=1+16+32         # 2^type for the types linked to nearest primaries


#--------------------------------------- Subfind
SUBFIND                                  # enables substructure finder


#--------------------------------------- output fields in snapshots--Default output filds are: position, velocity, ID, mass, spec. internal energy (gas), density(gas)
OUTPUTPOTENTIAL                          # output potential at particle position


#--------------------------------------- output options
PROCESS_TIMES_OF_OUTPUTLIST              # goes through times of output list prior to starting the simulaiton to ensure that outputs are written as close to the desired time as possible (as opposed to at next possible time if this flag is not active)
HAVE_HDF5                                # needed when HDF5 I/O support is desired (recommended)


#--------------------------------------- Testing and Debugging options
DEBUG                                    # enables core-dumps
