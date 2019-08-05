#!/bin/bash        	# this line only there to enable syntax highlighting in this file

#########################################################
#  Enable/Disable compile-time options as needed        #
#  examples/galaxy_merger_star_formation_3d/Config.sh   #
#########################################################


#--------------------------------------- Mesh motion and regularization
REGULARIZE_MESH_CM_DRIFT                 # Mesh regularization; Move mesh generating point towards center of mass to make cells rounder.
REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED  # Limit mesh regularization speed by local sound speed
REGULARIZE_MESH_FACE_ANGLE               # Use maximum face angle as roundness criterion in mesh regularization


#--------------------------------------- Refinement and derefinement
REFINEMENT_SPLIT_CELLS                   # Refinement
REFINEMENT_MERGE_CELLS                   # Derefinement
REFINEMENT_VOLUME_LIMIT                  # Limit the volume of cells and the maximum volume difference between neighboring cels
NODEREFINE_BACKGROUND_GRID               # Do not de-refine low-res gas cells in zoom simulations


#--------------------------------------- Time integration options
TREE_BASED_TIMESTEPS                     # non-local timestep criterion (take 'signal speed' into account)


#--------------------------------------- Gravity treatment
SELFGRAVITY                              # gravitational intraction between simulation particles/cells 	 
HIERARCHICAL_GRAVITY                     # use hierarchical splitting of the time integration of the gravity
CELL_CENTER_GRAVITY                      # uses geometric centers to calculate gravity of cells, only possible with HIERARCHICAL_GRAVITY
ALLOW_DIRECT_SUMMATION                   # Performed direct summation instead of tree-based gravity if number of active particles < DIRECT_SUMMATION_THRESHOLD (= 3000 unless specified differently here)
DIRECT_SUMMATION_THRESHOLD=500           # Overrides maximum number of active particles for which direct summation is performed instead of tree based calculation
GRAVITY_NOT_PERIODIC                     # gravity is not treated periodically


#--------------------------------------- Gravity softening
NSOFTTYPES=2                             # Number of different softening values to which particle types can be mapped.
MULTIPLE_NODE_SOFTENING                  # If a tree node is to be used which is softened, this is done with the softenings of its different mass components
INDIVIDUAL_GRAVITY_SOFTENING=32          # bitmask with particle types where the softenig type should be chosen with that of parttype 1 as a reference type
ADAPTIVE_HYDRO_SOFTENING                 # Adaptive softening of gas cells depending on their size


#--------------------------------------- Single/Double Precision
DOUBLEPRECISION=1                        # Mode of double precision: not defined: single; 1: full double precision 2: mixed, 3: mixed, fewer single precisions; unless short of memory, use 1.
NGB_TREE_DOUBLEPRECISION                 # if this is enabled, double precision is used for the neighbor node extension


#-------------------------------------------- Things for special behaviour
PROCESS_TIMES_OF_OUTPUTLIST              # goes through times of output list prior to starting the simulaiton to ensure that outputs are written as close to the desired time as possible (as opposed to at next possible time if this flag is not active)
OVERRIDE_PEANOGRID_WARNING               # don't stop if peanogrid is not fine enough

#--------------------------------------- Output/Input options
HAVE_HDF5                                # needed when HDF5 I/O support is desired (recommended)


#--------------------------------------- Testing and Debugging options
DEBUG                                    # enables core-dumps


#--------------------------------------- non-standard phyiscs
ENFORCE_JEANS_STABILITY_OF_CELLS         # this imposes an adaptive floor for the temperature
COOLING                                  # Simple primordial cooling
USE_SFR                                  # Star formation model, turning dense gas into collisionless partices
