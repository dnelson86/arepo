#!/bin/bash        	# this line only there to enable syntax highlighting in this file

##########################################################################
#  Enable/Disable compile-time options as needed                         #
#  examples/galaxy_merger_star_formation_3d/Config_ADDBACKGROUNDGRID.sh  #
##########################################################################


#--------------------------------------- Basic operation mode of code
ADDBACKGROUNDGRID=16                     # Re-grid hydrodynamics quantities on a Oct-tree AMR grid. This does not perform a simulation.


#--------------------------------------- Mesh motion and regularization
REGULARIZE_MESH_FACE_ANGLE               # Use maximum face angle as roundness criterion in mesh regularization


#--------------------------------------- Gravity softening
NSOFTTYPES=2                             # Number of different softening values to which particle types can be mapped.


#--------------------------------------- output options
HAVE_HDF5                                # needed when HDF5 I/O support is desired (recommended)

