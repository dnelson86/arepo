.. Arepo documentation master file
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

**********************
Arepo documentation
**********************

Arepo is a massively parallel code for gravitational N-body systems
and magnetohydrodynamics, both on Newtonian as well as cosmological
backgrounds.  It is a flexible code that can be applied to a variety
of different types of problems, offering a number of sophisticated
simulation algorithms. A description of the numerical algorithms
employed by the code is given in the public release code paper.  For a
more in depth discussion about these algorithms, the original code
paper and subsequent publications are the best resource. This
documentation only addresses the question how to use the different
numerical algorithms.

Arepo was written by Volker Springel (vspringel@mpa-garching.mpg.de)
with further development by Rüdiger Pakmor
(rpakmor@mpa-garching.mpg.de) and contributions by many other authors
(www.arepo-code.org/people).  The public version of the code was
compiled by Rainer Weinberger (rainer.weinberger@cfa.harvard.edu).

Overview
========

The Arepo code was initially developed to combine the advantages of
finite-volume hydrodynamics with the Lagrangian convenience of
smoothed particle hydrodynamics (SPH). To this end, Arepo makes use of
an unstructured Voronoi-mesh which is, in the standard mode of
operating the code, moving with the fluid in a quasi-Lagrangian
fashion. The fluxes between cells are computed using a finite-volume
approach, and additional spatial adaptivity is provided by the
possibility to add and remove cells from the mesh according to
user-defined criteria. In addition to gas dynamics, Arepo allows for
additional collisionless particle types which interact only
gravitationally. Besides self-gravity, forces from external
gravitational potentials can also be included.

Arepo is optimized for, but not limited to, cosmological simulations
of galaxy formation, and consequently can used for simulations with
high dynamic range in space and time. This is in particular supported
by Arepo's ability to employ adaptive local timestepping for each
individual cell or particle, as well as a dynamic load and memory
balancing domain decomposition. The current version of Arepo is fully
MPI parallel, and has been tested in runs with >10,000 MPI tasks. The
exact performance is, however, highly problem- and machine-dependent.


Disclaimer
==========

It is important to note that the performance and accuracy of the code
is a sensitive function of some of the code parameters. We also stress
that Arepo comes without any warranty, and without any guarantee that
it produces correct results. If in doubt about something, reading (and
potentially improving) the source code is always the best strategy to
understand what is going on!

**Please also note the following:**

The numerical parameter values used in the examples contained in the
code distribution do not represent a specific recommendation by the
authors! In particular, we do not endorse these parameter settings in
any way as standard values, nor do we claim that they will provide
fully converged results for the example problems, or for other initial
conditions. We think that it is extremely difficult to make general
statements about what accuracy is sufficient for certain scientific
goals, especially when one desires to achieve it with the smallest
possible computational effort. For this reason we refrain from making
such recommendations. We encourage every simulator to find out for
herself/himself what integration settings are needed to achieve
sufficient accuracy for the system under study. We strongly recommend
to make convergence and resolution studies to establish the range of
validity and the uncertainty of any numerical result obtained with
Arepo.


Table of contents
=================

.. toctree::
  running.md
  simtypes.md
  config-options.md
  parameterfile.md
  snapshotformat.md
  diagnosticfiles.md
  migration.md
  development.md
