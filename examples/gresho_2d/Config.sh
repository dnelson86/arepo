#!/bin/bash            # this line only there to enable syntax highlighting in this file

## examples/Gresho_2d/Config.sh
## config file for 2d Gresho vortex probelm


#--------------------------------------- Basic operation mode of code
TWODIMS                                  # 2d simulation
READ_MASS_AS_DENSITY_IN_INPUT            # Reads the mass field in the IC as density

#--------------------------------------- Mesh motion and regularization
REGULARIZE_MESH_CM_DRIFT                 # Mesh regularization; Move mesh generating point towards center of mass to make cells rounder.
REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED  # Limit mesh regularization speed by local sound speed
REGULARIZE_MESH_FACE_ANGLE               # Use maximum face angle as roundness criterion in mesh regularization

#--------------------------------------- Time integration options
TREE_BASED_TIMESTEPS                     # non-local timestep criterion (take 'signal speed' into account)

#---------------------------------------- Single/Double Precision
DOUBLEPRECISION=1                        # Mode of double precision: not defined: single; 1: full double precision 2: mixed, 3: mixed, fewer single precisions; unless short of memory, use 1.
INPUT_IN_DOUBLEPRECISION                 # initial conditions are in double precision
OUTPUT_IN_DOUBLEPRECISION                # snapshot files will be written in double precision
OUTPUT_CENTER_OF_MASS                    # output centers of cells

#--------------------------------------- Output/Input options
HAVE_HDF5                                # needed when HDF5 I/O support is desired (recommended)

#--------------------------------------- Testing and Debugging options
DEBUG                                    # enables core-dumps
