*************
Related codes
*************

The general high-level code structure of Arepo is quite similar to the
Gadget code, but with significant changes and improvements. However,
keeping in mind that the code architecture originates from an
Tree-PM-SPH code can help to understand why certain things are done in
certain ways.

One of the obvious signs of this are the names of some global
structures and functions. For example the array holding the
information for gas cells is called SphP. Another example would be the
routines in ``src/init/density.c``, which used to be an SPH density
calculation, but is now only used for an estimate for the typical
distance to neighbouring cells. Another aspect is the presence of the
neighbour tree, a search-tree structure on gas cells. In practice this
means that many of the routines that depend on the state of
neighbouring cells, including sub-grid feedback models etc., can be
based on this neighbour search. This makes it possible to port these
models relatively easily from earlier versions of Gadget to Arepo.

However, it is important to also realize that routines written for the
Gadget code will in general not work without some porting adjustments
when inserted in to Arepo. We thus encourage developers to carefully
test that their implementations still have the intended effect when
such a porting is done.  Apart from the differences in implementation,
it is also important to realize that a grid-based finite volume scheme
can react differently to certain sub-grid implementations than a
particle-based approach.


Differences to the development Version 
======================================

This public version of Arepo was branched off from the development
version in November 2017, and then Rainer Weinberger substantially cut
it down, cleaned it up and completed the documentation. The main
reason for doing this was to provide a code to the community that
researchers are able to understand and use without the need of an
extensive introduction by one of the authors. We believe that the
public version, which is better documented at a developer level, also
provides a better starting point for new developments.

The general idea for this public verson was to preserve the
well-tested code with a limited number of changes. However, some
compile-time options have been eliminated and to recover the same
running mode **in the development version of Arepo**, the following
compile-time flags need to be set in `Config.sh`

* ``PERIODIC`` (has become obsolete over the years)
* ``VORONOI`` (the development version also has an AMR mode)
* ``VORONOI_DYNAMIC_UPDATE`` (in practice now a default)
* ``AUTO_SWAP_ENDIAN_READIC`` (only relevant for file formats 1 and 2)
* ``CHUNKING``
* the switch ``EXACT_GRAVITY_FOR_PARTICLE_TYPE`` automatically
  activates ``NO_SELFGRAVITY_TYPE``, ``NO_GRAVITY_TYPE`` and
  ``EXACT_GRAVITY_REACTION`` (see src/main/allvars.h) in the public
  version. This is not the case in the development version where all
  these flags need to be set individually in Config.sh.

