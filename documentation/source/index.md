.. Arepo documentation master file
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

**********************
Arepo documentation
**********************

Arepo is a massively parallel code for gravitational n-body systems and 
magnetohydrodynamics, both on Newtonian as well as cosmological background.
It is a flexible code that can be applied to a variety of different types of
simulations, offering a number of sophisticated simulation algorithms. A
description of the numerical algorithms employed by the code is given in the 
public release paper. 
For a more detailed discussion about these algorithms, the original code paper 
and subsequent publications are the best resource. This documentation only 
addresses the question how to use the different numerical algorithms.

Arepo was written by Volker Springel (vspringel@mpa-garching.mpg.de) with further
development by Rüdiger Pakmor (rpakmor@mpa-garching.mpg.de) and contributions by
many other authors (www.arepo-code.org/people). 
The public version of the code was compiled by Rainer Weinberger 
(rainer.weinberger@cfa.harvard.edu).

Overview
========

The Arepo code was initially developed to combine the advantages of
finite-volume hydrodynamics schemes with the Lagrangian invariance of 
smoothed particle hydrodynamics (SPH) schemes. To this end, Arepo makes use of an
unstructured Voronoi-mesh which is, in its standard setting, moving with 
the fluid in an quasi-Lagrangian fashion. The fluxes between cells are computed
using a finite-volume approach, and further spatial adaptivity is 
provided by the possibility to add and remove cells from the 
mesh according to defined criteria. In addition to gas, Arepo allows for a 
number of additional particle types which interact only gravitationally, as 
well as for forces from external gravitational potentials.

Arepo is optimized for, but not limited to, cosmological simulations of galaxy 
formation and consequently for simulations with very high dynamic ranges in 
space and time. Therefore, Arepo employs an adaptive timestepping for each 
individual cell and particle as well as a dynamic load and memory 
balancing scheme. In its current version, Arepo is fully MPI parallel, and tested 
to run with >10,000 MPI tasks. The exact performance is, however, highly problem and
machine dependent. 


Disclaimer
==========

It is important to note that the performance and accuracy of the code is a
sensitive function of some of the code parameters. We also stress that Arepo
comes without any warranty, and without any guarantee that it produces correct
results. If in doubt about something, reading (and potentially improving) the
source code is always the best strategy to understand what is going on!

**Please also note the following:**

The numerical parameter values used in the examples contained in the code
distribution do not represent a specific recommendation by the authors! In
particular, we do not endorse these parameter settings in any way as standard
values, nor do we claim that they will provide fully converged results for the
example problems, or for other initial conditions. We think that it is extremely
difficult to make general statements about what accuracy is sufficient for
certain scientific goals, especially when one desires to achieve it with the
smallest possible computational effort. For this reason we refrain from making
such recommendations. We encourage every simulator to find out for
herself/himself what integration settings are needed to achieve sufficient
accuracy for the system under study. We strongly recommend to make convergence
and resolution studies to establish the range of validity and the uncertainty of
any numerical result obtained with Arepo.


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
