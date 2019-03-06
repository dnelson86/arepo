#!/bin/bash        	# this line only there to enable syntax highlighting in this file


##########################################################
## examples/cosmo_zoom_gravity_only_3d/Config_parent.sh ##
##########################################################


#--------------------------------------- Gravity treatment
SELFGRAVITY                              # gravitational intraction between simulation particles/cells
HIERARCHICAL_GRAVITY                     # use hierarchical splitting of the time integration of the gravity
ALLOW_DIRECT_SUMMATION                   # Performed direct summation instead of tree-based gravity if number of active particles < DIRECT_SUMMATION_THRESHOLD (= 3000 unless specified differently here)
DIRECT_SUMMATION_THRESHOLD=1024          # Overrides maximum number of active particles for which direct summation is performed instead of tree based calculation


#--------------------------------------- TreePM Options
PMGRID=256                               # Enables particle mesh; number of cells used for grid in each dimension


#--------------------------------------- Gravity softening
MULTIPLE_NODE_SOFTENING                  # If a tree node is to be used which is softened, this is done with the softenings of its different mass components


#--------------------------------------- Time integration options
TREE_BASED_TIMESTEPS                     # non-local timestep criterion (take 'signal speed' into account)


#--------------------------------------- Single/Double Precision
DOUBLEPRECISION=1                        # Mode of double precision: not defined: single; 1: full double precision 2: mixed, 3: mixed, fewer single precisions; unless short of memory, use 1.
DOUBLEPRECISION_FFTW                     # FFTW calculation in double precision


#--------------------------------------- On the fly FOF groupfinder
FOF                                      # enable FoF output
FOF_PRIMARY_LINK_TYPES=2                 # 2^type for the primary dark matter type


#--------------------------------------- Subfind
SUBFIND                                  # enables substructure finder


#--------------------------------------- output options
PROCESS_TIMES_OF_OUTPUTLIST              # goes through times of output list prior to starting the simulaiton to ensure that outputs are written as close to the desired time as possible (as opposed to at next possible time if this flag is not active)
HAVE_HDF5                                # needed when HDF5 I/O support is desired (recommended)


#--------------------------------------- Testing and Debugging options
DEBUG                                    # enables core-dumps
