*************
Related codes
*************

The general code structure of AREPO is similar to the GADGET code, but with 
significant further developments and improvements. However, keeping in mind
that the code architecture originates from an Tree-PM-SPH code helps to 
understand why certain things are done in certain ways.

One of the obvious signs of this are the name of some global structures and 
functions. For example the array holding the information for gas cells only is 
called SphP. Another example would be the routines in ``src/init/density.c``, 
which used to be an SPH density calculation, but now only used for an estimate
for the typical distance to neighbouring cells. Another aspect is the 
presence of the neighbour tree, a search-tree structure on gas cells. In 
practice this means that many of the routines that depend on the state of 
neighbouring cells, including sub-grid feedback models etc., can be based on
this neighbour search. This makes it possible to port these models from earlier
versions of GADGET to AREPO.

However, it is important to realize that routines written for the GADGET code
will in general not work directly when ported to Arepo. We encourage developers
to carefully test that the implementation is still having the same effect.
Apart from the differences in implementation, it is also important to realize 
that a grid-based finite volume scheme might react different to a particular to 
a sub-grid implementation.


Differences to the development Version 
======================================

This public version of AREPO was branched off from the development version in
November 2017, substantially cut down, cleand up and the documentation completed 
by Rainer Weinberger. The main reason for doing this is to provide a code to the 
community that researchers are able to understand and use without the need of 
an extensive introduction by one of the authors. We belive that the public 
version, being better documented on a developer level, also provides a better 
starting point for new developments. 

The general idea for this public verson was to preserve the well-tested code
with a limited number of changes. However, some compile-time options have been
eliminated and to recover the same running mode **in the development version 
of AREPO**, the following compile-time flags need to be set in `Config.sh`

* ``PERIODIC`` (has become obsolete over the years)
* ``VORONOI`` (the development version also has an AMR mode)
* ``VORONOI_DYNAMIC_UPDATE`` (in practice now a default)
* ``AUTO_SWAP_ENDIAN_READIC`` (only relevant for file formats 1 and 2)
* ``CHUNKING``
* the switch ``EXACT_GRAVITY_FOR_PARTICLE_TYPE`` automatically activates ``NO_SELFGRAVITY_TYPE``, ``NO_GRAVITY_TYPE`` and ``EXACT_GRAVITY_REACTION`` (see src/main/allvars.h) in the public version. This is not the case in the development version where all these flags need to be set individually in Config.sh.

