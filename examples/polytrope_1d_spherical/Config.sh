#!/bin/bash            # this line only there to enable syntax highlighting in this file

## ./examples/polytrope_1d_spherical/Config.sh
## config file for 1d polytrope probelm in a spherically symmetric setup


#--------------------------------------- Basic operation mode of code
ONEDIMS                                  # 1d simulation
ONEDIMS_SPHERICAL                        # 1d spherically symmetric simulation
GAMMA=2                                  # Adiabatic index of gas


#--------------------------------------- Mesh motion and regularization
VORONOI_STATIC_MESH                      # static mesh


#--------------------------------------- Gravity treatment; default: no gravity
SELFGRAVITY                              # gravitational intraction between simulation particles/cells


#--------------------------------------- Time integration options
FORCE_EQUAL_TIMESTEPS                    # variable but global timestep


#---------------------------------------- Single/Double Precision
DOUBLEPRECISION=1                        # Mode of double precision: not defined: single; 1: full double precision 2: mixed, 3: mixed, fewer single precisions; unless short of memory, use 1.
INPUT_IN_DOUBLEPRECISION                 # initial conditions are in double precision
OUTPUT_IN_DOUBLEPRECISION                # snapshot files will be written in double precision


#--------------------------------------- Output/Input options
READ_MASS_AS_DENSITY_IN_INPUT            # Reads the mass field in the IC as density
HAVE_HDF5                                # needed when HDF5 I/O support is desired (recommended)


#--------------------------------------- Testing and Debugging options
DEBUG                                    # enables core-dumps


#--------------------------------------- output
OUTPUTACCELERATION                       # output gravitational acceleration
OUTPUT_PRESSURE_GRADIENT                 # output pressure gradient
