#!/bin/bash            # this line only there to enable syntax highlighting in this file

## examples/shocktube_1d/Config.sh
## config file for 1d shocktube probelm


#--------------------------------------- Basic operation mode of code
ONEDIMS                                  # 1d simulation
REFLECTIVE_X=2                           # X Boundary; 1: Reflective, 2: Inflow/Outflow; not active: periodic
GAMMA=1.4                                # Adiabatic index of gas; 5/3 if not set

#--------------------------------------- Mesh motion and regularization
REGULARIZE_MESH_CM_DRIFT                 # Mesh regularization; Move mesh generating point towards center of mass to make cells rounder.
REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED  # Limit mesh regularization speed by local sound speed
REGULARIZE_MESH_FACE_ANGLE               # Use maximum face angle as roundness criterion in mesh regularization

#--------------------------------------- Time integration options
FORCE_EQUAL_TIMESTEPS                    # variable but global timestep

#---------------------------------------- Single/Double Precision
DOUBLEPRECISION=1                        # Mode of double precision: not defined: single; 1: full double precision 2: mixed, 3: mixed, fewer single precisions; unless short of memory, use 1.
INPUT_IN_DOUBLEPRECISION                 # initial conditions are in double precision
OUTPUT_IN_DOUBLEPRECISION                # snapshot files will be written in double precision

#--------------------------------------- Output/Input options
HAVE_HDF5                                # needed when HDF5 I/O support is desired (recommended)

#--------------------------------------- Testing and Debugging options
DEBUG                                    # enables core-dumps

