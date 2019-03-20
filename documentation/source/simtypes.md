Simulation examples
**************************

AREPO is a multi-purpose code that supports a number of different types of 
simulations. Here, we present some examples that in our opinion are 
particularly relevant. This list is by no means complete, and highly biased 
to the fields of interest of the authors. 
We encourage users that have other simulation setups to make them (or small
test examples) availlable to us, in order to provide a more exhaustive list of 
examples in future code versions.


General Execution Structure of Examples
=======================================

Examples can be run using the ``test.sh`` shell script (typing ``bash test.sh``). 
The script itself, as the name suggests, is there it run a couple
of verification simulations and report whether the simulation 
results agree with specifications. In particular, if this is the case, 
the script will delete the directory in which the example run in, which 
might not be desired when getting started with the code.

The general procedure to execute an example is (replace `myexample` by the desired one):

 * create a run-directory, in which the run will be executed::
 
    mkdir -p ./run/examples/myexample
    
 * copy the entire content of the example directory of interest to the run directory::
 
    cp -r ./examples/myexample/* ./run/examples/myexample/
    
 * execute the ``create.py`` file, with the run directory as first argument::
 
    python ./run/examples/myexample/create.py \
    ./run/examples/myexample/
    
  This creates initial conditions and possible other files required by the simulation.
    
 * compile AREPO, using the provided Config.sh file::
 
    make CONFIG=./run/examples/myexample/Config.sh \
    BUILD_DIR=./run/examples/myexample/build \
    EXEC=./run/examples/myexample/Arepo
    
 * change to run directory::
 
    cd ./run/examples/myexample/
   
 * execute Arepo (here on 8 mpi ranks)::
 
    mpiexec -np 8 ./Arepo param.txt
    
 * change back to root directory::
 
    cd ../../../
    
 * execute the ``check.py`` script to do analysis on the simulation output (again specifying the relative path to the run directory)::
 
    python ./run/examples/myexample/check.py \
    ./run/examples/myexample/
    
    
Note that the only command that is crucial to be executed in a 
specific directory is the execution of Arepo ``mpiexec -np 8 ./Arepo param.txt``,
because the parameter file contains relative paths of the initial 
conditions, other required files, and the output directory. All other 
commands can be executed from an arbitrary directory, as long as
the relative path handed over is correct.


Cosmological Simulations
========================

One for the main simulation types in AREPO are simulations on an expanding 
spacetime, where the coordinates are treated as comoving coordinates.
This mode is active whenever the parameter flag ``ComovingIntegrationOn`` is set
to "1" in the parameter file. In this case the code models the expanding 
spacetime discribed by the Friedmann equations.

The halo (FoF) and subhalo finder (subfind) are analysis tools mainly for these
types of simulatons.

Cosmological volumes
--------------------

One of the standard types of cosmological simulations are 
cosmological volume simulations. These simulatons have a uniform mass 
resolution and are set up using cosmological perturbation theory, converted
to small initial position and velocity perturbations in a periodic box
with fixed comoving extent. These simulatios can be done with gravity only, 
in which case the simulation particles are all of a single type (type > 0) and 
gravity is recommended to be calculated using a TreePM
algorithm. This means that the compile time flags ``SELFGRAVITY`` and ``PMGRID`` 
are set. The example ``cosmo_box_gravity_only_3d`` is a very small 
representative of such a simulation.

In these cosmological simulations it is also possible to include gas.
This can be done, for example, by reading in gravity only initial conditions
and setting the compile-time option ``GENERATE_GAS_IN_ICS``, which  will split 
every input particle in a gas cell and a dark matter particle,
with the mass ratio following the cosmic average mass and baryon fractions,
``Omega0`` and ``OmegaBaryon``, respectively, which are defined in the parameter file. 
Including gas in the simulation, it is also useful to use mesh regularization 
options as well as refinement and derefinement. On top of this, it is also 
possible to use prescriptions for radiative cooling and star formation using 
the respective compile-time flags ``COOLING`` and ``USE_SFR``. The example 
``cosmo_box_star_formation_3d`` is using such a setup.

Cosmological zoom simulations
-----------------------------

A second class of cosmological simulations that have become very popular in 
studies of galaxy and galaxy cluster formation are so-called zoom simulations,
in which the Lagrangian region belonging to a single final object is 
sampled with significantly higher resolution than the rest of the cosmological 
volume. This way, it is possible to accurately resolve the formation of an
individual object, while still modelling the large-scale environmental 
effects. These simulations are, as the above ones, possible to run in gravity
only mode or including gas physics. A gravity only example is given as
``cosmo_zoom_gravity_only_3d`` in the examples (note however that the ICs need to be created
separately in this example, as they are slightly too large to inclue them 
into the main code repository).
Technically, such a zoom simulation needs slightly different configuration 
in the particle-mesh algorithm. In particular, it is useful to place a second
particle-mesh layer which calculates forces only in the high-resolution region.
This can be triggered by the compile-time flag ``PLACEHIGHRESREGION``. 
In addition, the flag ``PM_ZOOM_OPTIMIZED`` will change the employed 
particle-mesh algorithm to one that is specifically optimized for this type of 
simulations.

Initial conditions generation
-----------------------------

The creation of initial conditions for cosmological simulation is in general
a complicated procedure, and not included in the AREPO code. However, we included
in the examples the possibility to download and execute third-party initial condition
generating software to illustrate the general procedure how to create them
(MUSIC, Hahn&Abel 2011, MNRAS, 415, 2101 and N-GenIC, 
Springel et al. 2005, Nature, 435, 629) 
We note that these third party software packages generally 
have different library requirements than Arepo, and therefore might 
not work immediately on some machines. We also note that the
initial conditions used in the examples are choosen for illustration 
purposes and to minimize computational effort to run them. They are
generally not suited to be used directly in scientific work.


Newtonian space
===============

Apart from comoving integration, AREPO can also handle ordinary Newtonian 
spacetime by coosing the parameter option ``ComovingIntegrationOn 0``. 
While cosmological simulations usually assume periodic boundary conditions, 
simulations in Newtonian space can also have reflective or inflow/outflow
boundary conditions. 

The first example of this class of simulations is an isolated, self-gravitating
object, such as in ``isolated_galaxy_collisionless_3d``. This particular case
only contains collisionless particles, namely the dark matter and stellar 
component of a disk galaxy, and their gravitational interactions
are calculated with a tree algorithm only, assuming non-periodic forces.
One specific compile-time flag recommended in these kind of simulations is 
``RANDOMIZE_DOMAINCENTER``.  Without this flag it is possible that correlated 
force-errors cause an isolated galaxy to drift from the initial center of 
mass, which is prevented by de-correlating these errors by frequently choosing 
a new domain center. The initial conditions of this problem are created with
the GalIC code (Yurin&Springel 2014, MNRAS, 444, 62). The initial conditions
generation can be triggered in the create.py script, however is a computationally
expensive procedure. Therefore we also provide the initial conditions in the
repository.

Another popular type of simulations in galactic astrophysics are mergers of
galaxies. The example ``galaxy_merger_star_formation_3d`` shows such an example, 
which is similar to ``isolated_galaxy_collisionless_3d``, however in this 
particular case includes gas and two galaxies interacting. AREPO, as every 
grid code, requires a finite extent of the simulation box as well as non-zero 
density at every point.
This is different from Smoothed-Particle-Hydrodynamics simulations such as 
the ones performed with GADGET, where simply not placing gas particles in
the initial condition is an acceptable solution to treat low-density regions.
The initial condition for ``galaxy_merger_sfr_3d`` are, for historic reasons,
Smoothed-Particle-Hydrodynamics initial conditions. To be able to use them
in AREPO, the compile-time opiton ``ADDBACKGROUNDGRID`` was introduced.
With this mode enabled, the code does not perform a simulation, but converts 
the SPH gas particles into a hierarchical oct-tree structure of cells.
Once this is done, this output (using the recommended box size) can be used
for starting the actual simulation. Note that for this type of simulation,
we do not provide the initial-condition generating code. Generally, using
grid-code specific initial conditions is preferred over converting particle-based 
initial conditions, since the latter will always contain numerical 
noise which degrades the quality of the simulation.

Other simulations might be completely hydrodynamical, with no gravitational 
forces involved. One example of such a setup is the ``Noh_3d`` problem.
This example is more a test problem as it has an analytic solution. Thus it 
is better suited for code verification. However, also for the more complex
examples, we do check against outcomes from previous code versions in order
to ensure consistency.

Two-dimensional simulations
===========================

Two-dimensional simulations are often used as simplified examples which are
significantly cheaper to run than their three-dimensional counterparts, but 
also very useful as test problems. In particular, while the one-dimensional 
simulation code is in substantial parts detached from the rest, the two 
dimensional and three dimensional simulations largely use the same routines.
Thus these routines can be efficiently tested with 2d test problems.

One example of such a test is ``gresho_2d``, a stationary vortex problem for 
which the pressure gradient balances the centrifugal forces. The same is
true for the Yee vortex, ``yee_2d``, which has the advantage of being smooth
in the hydrodynamic properties and therefore better suited to determine the
convergence order of the code. Both these tests are sensitive to accurate 
hydrodynamical modelling and gradient estimates.

Other 2d test problems are ``noh_2d``, a converging gas flow in 2d which 
introduces a strong shock and is a challanging problem for the Riemann solver,
as well as ``current_sheet_2d``, an MHD test probing numerical reconnection 
properties of the code.


One-dimensional simulations
===========================

Most of the one-dimensional simulations are test problems for particular 
solvers or numerical methods. Note that Arepo is **NOT** MPI parallel in the 
case of 1d problems, i.e. all the following examples can only be calculated
on one MPI rank.

One of the most basic one-dimensional examples is a simple linear
wave propagation test as in the example ``wave_1d``. Such a test is well-suited
to test the convergence order of a scheme. For such applications it is 
important to ensure that both input/output and calculation are done in double 
precision by using the compile-time flags ``INPUT_IN_DOUBLEPRECISION``,
``OUTPUT_IN_DOUBLEPRECISION`` and ``DOUBLEPRECISION=1``.

Another very important basic test is a 1d shocktube problem. Such a shocktube
is calculated in the example ``shocktube_1d``, where the time evolution can 
be compared to the exact solution. This tests both the Riemann solver as well
as the gradient estimates and associated limiters. Similarly, ``mhd_shocktube_1d``
tests the time evolution of an MHD Riemann problem. Note that so far, the 
verification in the latter case is only against a high-resolution simulation.

One example of a more complex one-dimensional test problem is the example 
``interacting_blastwaves_1d``. This is checked against a high-resolution solution
calculated on a fixed grid.


One-dimensional simulations in spherical symmetry
=================================================

An important set of simulations for stellar astrophysics are one dimensional 
simulations in spherical symmetry, such as given for example in 
``polytrope_1d_spherical``. Generally, AREPO is not optimized for such kind of 
simulations and can only run these in serial. However, this option has proven 
to be useful for quick test simulations as well as for calculating radial 
profiles which are then used as initial conditions for three dimensional 
simulations. One of the key aspect for these simulation is the quality of the
boundary conditions. In this spherically symmetric mode, AREPO uses a 
reflective inner boundary condition and an inflow/outflow boundary condition
at the outer end. 


