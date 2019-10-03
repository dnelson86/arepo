Code Configuration
*************************

Basic operation mode
============================

The default running mode (without any of the flags active) is 3d with
6 particle types; type 0 is always gas; types >0 are only
gravitationally interacting.


**NTYPES=6**

The number of particle types used. Minimum: 6.

-----

**TWODIMS**

Simulation in 2d. Z coordinates and velocities are set to zero after
reading in initial conditions.

-----

**ONEDIMS**

Simulation in 1d. Y and Z coordinates and velocities are set to zero
after reading in initial conditions. For one dimensional
simulations, refinement and derefinement is not supported.
If this flag is active, the code is not MPI parallel.

-----

**ONEDIMS_SPHERICAL**

Spherically symmetric 1d simulation. Use together with ``ONEDIMS``.
The first dimension is used as the radial coordinate.

-----

Computational box
================================

The default running mode (without any of the flags active) is a cubic
box with periodic boundary conditions

**LONG_X=10.0**

These options can be used to distort the simulation cube along the
given direction with the given factor into a parallelepiped of
arbitrary aspect ratio. The box size in the given direction increases
by the factor given (e.g. if ``Boxsize`` is set to 100 and
``LONG_X=4`` is set the simulation domain extends from 0 to 400 along
X and from 0 to 100 along Y and Z.)

-----

**LONG_Y=2.0**

Stretches the y extent of the computational box by a given factor.

-----

**LONG_Z=10.0**

Stretches the z extent of the computational box by a given factor.

-----

**REFLECTIVE_X=1**

Boundary conditions in the x direction. 1: Reflective, 2:
Inflow/Outflow; not set: periodic

-----

**REFLECTIVE_Y=1**

Boundary conditions in y direction. 1: Reflective, 2: Inflow/Outflow;
not set: periodic

-----

**REFLECTIVE_Z=1**

Boundary conditions in z direction. 1: Reflective, 2: Inflow/Outflow;
not set: periodic

-----

Hydrodynamics
=============

The default mode is: ``GAMMA=5/3`` ideal hydrodynamics

**NOHYDRO**

No hydrodynamics calculation. Note that simply not including any type
0 particles has the same effect.

-----

**GAMMA=1.4**

Adiabatic index of gas. 5/3 if not set.

-----

**ISOTHERM_EQS**

Isothermal gas. Code uses an isothermal Riemann-solver.

-----

**PASSIVE_SCALARS=3**

Number of passive scalar fields advected with fluid (default: 0).

-----

**NO_SCALAR_GRADIENTS**

Disables time and spatial extrapolation for passive scalar fields. Use
only if you know why you are doing this.

-----

Magnetohydrodynamics
====================

By default, code only computes hydrodynamics. Note that for comparison
of MHD and hydrodynamical runs, it is sometimes useful to keep the MHD
settings active and to initialize the magnetic field to zero
everywhere. The equations of ideal MHD ensure that the magnetic field
stays exactly zero throughout the calculation.

**MHD**

Master switch for magnetohydrodynamics.

-----

**MHD_POWELL**

Powell div(B) cleaning scheme for magnetohydrodynamics.

-----

**MHD_POWELL_LIMIT_TIMESTEP**

Additional timestep constraint due to Powell cleaning scheme.

-----

**MHD_SEEDFIELD**

Uniform magnetic seed field of specified orientation and strength set
up after reading in IC.

-----

Riemann solver
==============

By default, an iterative, exact (hydrodynamics) Riemann solver is
used. If one of the flags below is active, this is changed. Only one
Riemann solver can be active.

**RIEMANN_HLLC**

HLLC approximate Riemann solver.

-----

**RIEMANN_HLLD**

HLLD approximate Riemann solver (required for MHD).

-----

Mesh motion 
==============================

The default mode is a moving mesh.

**VORONOI_STATIC_MESH**

Assumes the mesh to be static, i.e. to not change with time. The
vertex velocities of all mesh-generating points is set to zero and
domain decomposition is disabled.

-----

**VORONOI_STATIC_MESH_DO_DOMAIN_DECOMPOSITION**

Enables domain decomposition together with ``VORONOI_STATIC_MESH``
(which is otherwise then disabled), in case non-gas particle types
exist and the use of domain decompositions is desired. Note that on
one hand it may be advantageous in case the non-gas particles mix well
or cluster strongly, but on the other hand the mesh construction that
follows the domain decomposition is slow for a static mesh, so whether
or not using this new flag is overall advantageous depends on the
problem.

-----

**REGULARIZE_MESH_CM_DRIFT**

Mesh regularization. Move mesh generating point towards center of mass
to make cells rounder.

-----

**REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED**

Limits mesh regularization speed by local sound speed.

------

**REGULARIZE_MESH_FACE_ANGLE**

Uses maximum face angle as roundness criterion in mesh regularization.

-----

Refinement 
===========================

By default, there is no refinement and derefinement. But if enabled,
and unless set otherwise, the criterion for refinement/derefinement is
to maintain a cell target mass.

**REFINEMENT_SPLIT_CELLS**

Allows refinement.

-----

**REFINEMENT_MERGE_CELLS**

Allows derefinement.

-----

**REFINEMENT_VOLUME_LIMIT**

Limits the volume of cells and the maximum volume difference between
neighboring cells.

-----

**JEANS_REFINEMENT**

Refinement criterion to ensure resolving the Jeans length of cells.

-----

**REFINEMENT_HIGH_RES_GAS**

Limits the dynamical (de-)refinements of cells to cells which are
either already present in the ICs or are created with
``GENERATE_GAS_IN_ICS`` from type 1 particles. This adds an additional
integer quantity ``AllowRefinement`` to PartType0 in the snapshots
indicating if a gas cell is allowed to be refined and if it is, how
often this cell has already been split: if 0, no splitting allowed. If
odd (starting at 1), the cell was already present in the ICs. If even
(starting at 2), the cell was generated from a type 1 particle.  For
values of 3 or more, ``floor((AllowRefinement-1)/2.0)`` gives the
number of times the cell was split.

-----

**NODEREFINE_BACKGROUND_GRID**

The background grid will be prevented from derefining, when refinement
is used. In practice, enabling this option requires an input parameter
``MeanVolume``. Derefinement is then disallowed during the run for all
cells with ``Volume > 0.1 * MeanVolume``.

-----

**OPTIMIZE_MESH_MEMORY_FOR_REFINEMENT**

If activated some grid structures not needed for mesh refinement or
derefinement are freed before the function for refinement and
derefinement is called. The remaining mesh structures are freed after
this step as usual.

-----

Non-standard physics
====================

**COOLING**

Simple primordial cooling routine.

-----

**ENFORCE_JEANS_STABILITY_OF_CELLS**

This imposes an adaptive floor for the temperature.

-----

**USE_SFR**

Star formation model, turning dense gas into collisionless
particles. See Springel & Hernquist, (2003, MNRAS, 339, 289)

-----

**SFR_KEEP_CELLS**

Do not destroy cell out of which a star has formed.

-----

Gravity 
=================

If nothing is actived in this section, gravity is not included.

**SELFGRAVITY**

Computes gravitational interactions between simulation particles and
cells.

-----

**HIERARCHICAL_GRAVITY**

Uses hierarchical splitting of the time integration of the gravity.

-----

**CELL_CENTER_GRAVITY**

Uses geometric centers (instead of mesh-generating points) to
calculate gravity of cells, only possible with
``HIERARCHICAL_GRAVITY``.

-----

**NO_GAS_SELFGRAVITY**

Switches off gas self-gravity in tree.

-----

**GRAVITY_NOT_PERIODIC**

Gravity is not treated periodically.

-----

**ALLOW_DIRECT_SUMMATION**

Performs direct summation instead of tree-based gravity if number of
active particles < ``DIRECT_SUMMATION_THRESHOLD`` (= 3000 unless
specified differently)

-----

**DIRECT_SUMMATION_THRESHOLD=1000**

Overrides maximum number of active particles for which direct
summation is performed instead of a tree based calculation.

-----

**EXACT_GRAVITY_FOR_PARTICLE_TYPE=4**

Enables direct summation gravity calculation for the given particle
type.

-----

**EVALPOTENTIAL**

When this option is set, the code will compute the gravitational
potential energy each time a global statistics is computed. This can
be useful for testing global energy conservation.

-----


TreePM 
==============

If no option is actived here: no Particle-Mesh calculation is done.

**PMGRID=512**

Dimension of particle-mesh grid covering the domain.  This enables the
TreePM method, i.e. the long-range force is computed with a
PM-algorithm, and the short range force with the tree.  The parameter
has to be set to the size of the mesh that should be used, e.g. 256,
512, 1024 etc. The mesh dimensions need not necessarily be a power of
two, but the FFT is fastest for such a choice.  Note: If the
simulation is not in a periodic box, then a FFT method for vacuum
boundaries is employed, using a mesh with dimension twice that
specified by ``PMGRID``. Should not be used with a mesh much smaller
than 256, because the TreePM approximation is only valid if the range
of the tree calculation is small compared to the box size.

-----

**ASMTH=1.25**

This factor expressed the adopted force split scale in the TreePM
approach in units of the grid cell size. Setting this value overrides
the default value of 1.25, in mesh-cells, which defines the
long-range/short-range force split.

-----

**RCUT=6.0**

This determines the maximum radius, in units of the force split scale,
out to which the tree calculation in TreePM mode considers tree
nodes. If a tree node is more distant, the corresponding branch is
discarded. The default value is
4.5, given in mesh-cells.

-----

**PM_ZOOM_OPTIMIZED**

This option enables a different communication algorithm in the PM
calculations which works well independent of the data layout, in
particular it can cope well with highly clustered particle
distributions that occupy only a small subset of the total simulated
volume. However, this method is a bit slower than the default approach
(used when the option is disabled), which is best matched for
homogeneously sampled periodic boxes.

-----

**PLACEHIGHRESREGION=2**

If this option is set (will only work together with ``PMGRID``), then
the long range force is computed in two stages: One Fourier-grid is
used to cover the whole simulation volume, allowing the computation of
the large-scale force.  A second Fourier mesh is placed on the region
occupied by "high-resolution" particles, allowing the computation of
an intermediate-scale force. Finally, the force on very small scales
is computed by the tree. This procedure can be useful for
"zoom-simulations", where the majority of particles (the high-res
particles) are occupying only a small fraction of the volume. To
activate this option, the parameter needs to be set to an integer that
encodes the particle type(s) that make up the high-res particles in
the form of a bit mask. For example, if types 0, 1, and 4 are the
high-res particles, then the parameter should be set to
``PLACEHIGHRESREGION=1+2+16``, i.e. to the sum 2^0 + 2^1 + 2^4. The
spatial region covered by the high-res grid is determined
automatically from the initial conditions. The region is recalculated
if one of the selected particles is falling outside of the
high-resolution region. Note: If a periodic box is used, the high-res
zone is not allowed to intersect the box boundaries.

-----

**ENLARGEREGION=1.1**

This is only relevant when ``PLACEHIGHRESREGION`` is activated. The
size of the high resolution box will be automatically determined as
the minimum size required to contain the selected particle type(s), in
a "shrink-wrap" fashion.  This region is expanded on the fly, if
needed (see above). However, in order to prevent a situation where
this size needs to be enlarged frequently, such as when the particle
set is (slowly) expanding, the minimum size is multiplied by the
factor ``ENLARGEREGION`` (if defined). Then even if the set is
expanding, this will only rarely trigger a recalculation of the high
resolution mesh geometry, which is in general also associated with a
change of the force split scale.

-----

**GRIDBOOST=2**

Normally, if ``PLACEHIGHRESREGION`` is enabled, the code will try to
employ an effective grid size for the high-resolution patch that is
equivalent to ``PMGRID``. Because zero-padding has to be used for the
high-res inset, this gives a total mesh twice as large, which
corresponds to ``GRIDBOOST=2``. This value can here be modified by
hand, to e.g. 1, 3, 4, etc., to decrease or increase the size of
the high-res PM grid relative to that covering the full box. The total
mesh size used for the high-resolution FFTs is given by
``GRIDBOOST*PMGRID``.


-----

**FFT_COLUMN_BASED**

When this is enabled, the FFT calculations are not parallelized in
terms of a slab-decomposition but rather through a column based
approach. This scales to larger number of MPI ranks but is slower in
absolute terms as twice as many transpose operations need to be
performed. It is hence only worthwhile to use this option for a very
large number of MPI ranks that exceeds the 1D mesh dimension.

-----

Gravity softening
=================

In the default configuration, the code uses a small table of possible
gravitational softening lengths, which are specified in the parameter
file through the ``SofteningComovingTypeX`` and
``SofteningMaxPhysTypeX`` options, where X is an integer that gives
the "softening type". Each particle type is mapped to one of these
softening types through the ``SofteningTypeOfPartTypeY`` parameters,
where ``Y`` gives the particle type.  The number of particle types and
the number of softening types do not necessarily have to be
equal. Several particle types can be mapped to the same softening if
desired.

**NSOFTTYPES=4**

This can be changed to modify the number of available softening
types. These must be explicitly input as SofteningComovingTypeX
parameters, and so the value of ``NSOFTTYPES`` must match the number
of these entries in the parameter file.

-----

**MULTIPLE_NODE_SOFTENING**

If the tree walk wants to use a 'softened node' (i.e. where the
maximum gravitational softening of some particles in the node is
larger than the node distance and larger than the target particle's
softening), the node is opened by default (because there could be mass
components with a still smaller softening hidden in the node). This
can cause a substantial performance penalty in some cases. By setting
this option, this can be avoided. The code will then be allowed to use
softened nodes, but it does that by evaluating the node-particle
interaction for each mass component with different softening type
separately (but by neglecting possible shifts in their centers of
masses).  This also requires that each tree node computes and stores a
vector with these different masses. It is therefore desirable to not
make the table of softening types excessively large. This option can
be combined with adaptive hydro softening. In this case, particle type
0 needs to be mapped to softening type 0 in the parameter file, and no
other particle type may be mapped to softening type 0 (the code will
issue an error message if one doesn't obey to this).

-----

**INDIVIDUAL_GRAVITY_SOFTENING=2+4**

The code can also be asked to set the softening types of some of the
particle types automatically based on particle mass. The particle
types to which this is applied are set by this compile time option
through a bitmask encoding the types. The code by default assumes that
the softening of particle type 1 should be the reference. To this end,
the code determines the average mass of type 1 particles, and the
types selected through this option then compute a desired softening
length by scaling the type-1 softening with the cube root of the mass
ratio. Then, the softening type that is closest to this desired
softening is assigned to the particle (*choosing only from those
softening values explicitly input as a SofteningComovingTypeX
parameter*). This option is primarily useful for zoom simulations,
where one may for example lump all boundary dark matter particles
together into type 2 or 3, but yet provide a set of softening types
over which they are automatically distributed according to their
mass. If both ``ADAPTIVE_HYDRO_SOFTENING`` and
``MULTIPLE_NODE_SOFTENING`` are set, the softening types considered
for assignment exclude softening type 0. Note: particles that accrete
matter (black holes or sinks) get their softening updated if needed.

-----

**ADAPTIVE_HYDRO_SOFTENING**

When this is enabled, the gravitational softening lengths of hydro
cells are varied along with their radius. To this end, the radius of a
cell is multiplied by the parameter ``GasSoftFactor``. Then, the
closest softening from a logarithmically spaced table of possible
softenings is adopted for the cell. The minimum softening in the table
is specified by the parameter ``MinimumComovingHydroSoftening``, and
the larger ones are spaced a factor ``AdaptiveHydroSofteningSpacing``
apart. The resulting minimum and maximum softening values are reported
in the stdout log file.

-----

**NSOFTTYPES_HYDRO=64**

This is only relevant if ``ADAPTIVE_HYDRO_SOFTENING`` is enabled and
can be set to override the default value of 64 for the length of the
logarithmically spaced softening table. The sum of ``NSOFTTYPES`` and
``NSOFTTYPES_HYDRO`` may not exceed 254 (this is checked).

-----

External gravity
================

By default, there is no external gravitational potential.

**EXTERNALGRAVITY**

Master switch for external potential.

-----

**EXTERNALGY=0.0**

Constant external gravity in the y-direction

-----

NFW Potential
--------------------

**STATICNFW**

Static gravitational Navarro-Frenk-White (NFW) potential.

-----

**NFW_C=12**

Concentration parameter of NFW potential.

-----

**NFW_M200=100.0**

Mass causing the NFW potential.

-----

**NFW_Eps=0.01**

Softening of NFW potential.

-----

**NFW_DARKFRACTION=0.87**

Fraction of dark matter in NFW potential. The potential will be
reduced by this factor (with the idea being that the complement is
respresented by gas mass included explicitly in the simulation).

-----

Isothermal Sphere
----------------------------------

**STATICISO**

Static gravitational isothermal sphere potential.

-----

**ISO_M200=100.0**

Mass causing the isothermal sphere potential.

-----

**ISO_R200=160.0**

Radius of the isothermal sphere potential.

-----

**ISO_Eps=0.1**

Softening of isothermal sphere potential.

-----

**ISO_FRACTION=0.9**

Fraction in dark matter in isothermal sphere potential. Potential will
be reduced by this factor.

-----

Hernquist Potential
--------------------------

**STATICHQ**

Static gravitational Hernquist potential.

-----

**HQ_M200=186.015773**

Mass causing the Hernquist potential.

-----

**HQ_C=10.0**

Concentration parameter of Hernquist potential.

-----

**HQ_DARKFRACTION=0.9**

Fraction in dark matter in Hernquist potential. Potential will be
reduced by this factor.

-----

Time integration
========================

**FORCE_EQUAL_TIMESTEPS**

Variable but global timestep. Here the tightest timestep criterion
evaluated for any of the particles determines the timestep of all
particles.

-----

**TREE_BASED_TIMESTEPS**

Non-local timestep criterion (which takes the 'signal speed' of
hydrodynamical waves arriving from any point into account).

-----

**PM_TIMESTEP_BASED_ON_TYPES=2+4**

Particle types that should be considered in setting the PM timestep.

-----

**NO_PMFORCE_IN_SHORT_RANGE_TIMESTEP**

PM force is not included in short-range timestep criterion.

-----

**ENLARGE_DYNAMIC_RANGE_IN_TIME**

This extends the dynamic range of the integer timeline from 32 to 64
bits.

-----

Message Passing Interface
========================================

**IMPOSE_PINNING**

Enforce pinning of MPI tasks to cores if MPI does not do it.

-----

**IMPOSE_PINNING_OVERRIDE_MODE**

Override MPI pinning, if present.

-----

Single/Double Precision 
=======================

**DOUBLEPRECISION=1**

Mode of numerical precision: not set: single; 1: full double precision
2: mixed, 3: mixed, fewer single precisions; unless extremely short of
memory, we recommend to always use 1.

-----

**DOUBLEPRECISION_FFTW**

FFTW calculation in double precision.

-----

**OUTPUT_IN_DOUBLEPRECISION**

Snapshot files will be written in double precision.

-----

**INPUT_IN_DOUBLEPRECISION**

Initial conditions are in double precision.

-----

**OUTPUT_COORDINATES_IN_DOUBLEPRECISION**

Will always output coordinates in double precision.

-----

**NGB_TREE_DOUBLEPRECISION**

If this is enabled, double precision is aslo used for storing the
spatial neighbor node extension (the precision requirements for this
are less demanding than for other quantities).

-----

Groupfinder
========================================

**FOF**

Master switch to enable the friends-of-friends group finder in the
code. This will then usually be applied automatically before snapshot
files are written (unless disabled selectively for certain output
dumps).

-----

**FOF_PRIMARY_LINK_TYPES=2**

This option selects the particle types that are processed by the
friends-of-friends linking algorithm. A default linking length of 0.2
is assumed for this particle type unless specified otherwise.  The
specified value corresponds to Sum(2^type) for the primary dark matter
type(s).

-----

**FOF_SECONDARY_LINK_TYPES=1+16+32**

With this option, FOF groups can be augmented by particles/cells of
other particle types that they "enclose". To this end, for each
particle among the types selected by the bit mask specified with
``FOF_SECONDARY_LINK_TYPES``, the nearest among
``FOF_PRIMARY_LINK_TYPES`` is found and then the particle is attached
to whatever group this particle is in. The specified values
corresponds to sum(2^type) for the types linked to nearest primaries.

-----

**FOF_SECONDARY_LINK_TARGET_TYPES= 2**

An option to make the secondary linking work better in zoom runs
(after the FOF groups have been found, the tree is newly constructed
for all the secondary link targets). This should normally be set to
all dark matter particle types. If not set, it defaults to
``FOF_PRIMARY_LINK_TYPES``, which reproduces the old behavior.

-----

**FOF_GROUP_MIN_LEN=32**

Minimum number of particles (primary+secondary) in one group (default
is 32).

-----

**FOF_LINKLENGTH=0.16**

Linking length for FoF in units of the mean inter-particle separation
(default=0.2).

-----

**FOF_STOREIDS**

Normally, the snapshots produced with a FOF group catalogue are stored
in group order, such that the particle set making up a group can be
inferred as a contiguous block of particles in the snapshot file,
making it redundant to separately store the IDs of the particles
forming a group in the group catalogue. By activating this option, one
can nevertheless enforce the creation of the corresponding lists of
IDs as part of the group catalogue output.

-----

Subfind 
===========================

**SUBFIND**

When enabled, this automatically applies the Subfind subtructure
finder to all FOF groups after they have been found. Also, the
snapshot files are brought into subhalo order within each group.

-----

**SAVE_HSML_IN_SNAPSHOT**

When activated, this will store the smoothing kernel lengths used for
estimating the total matter density around every point, and the
corresponding densities, in the snapshot files associated with a run
of Subfind.

-----

**SUBFIND_CALC_MORE**

Additional calculations are carried out in the Subfind algorithm,
which are not always needed.
(i) The velocity dispersion in the kernel volume used for estimating
the local density.
(ii) The DM density around every particle is stored in the snapshot if
this is set together with ``SAVE_HSML_IN_SNAPSHOT``.

-----

**SUBFIND_EXTENDED_PROPERTIES**

Additional calculations are carried out in the Subfind algorithm,
which are not always needed and may be expensive.
(i) Further quantities related to the angular momentum in different
components.
(ii) The kinetic, thermal and potential binding energies for spherical
overdensity halos.


-----

Special behavior 
============================

**RUNNING_SAFETY_FILE**

If file './running' exists, do not start the run. Can be used to
prevent that a simulation is executed twice at the same time.

-----

**MULTIPLE_RESTARTS**

Keep several restart files instead of just last two copies.

-----

**EXTENDED_GHOST_SEARCH**

This extends the ghost search to the full 3x3 domain instead of the
principal domain. This can be needed for a successful mesh
construction if the box is sampled only with a couple of cells per
dimension.

-----

**DOUBLE_STENCIL**

This will ensure that the boundary region of the local mesh is deep
enough to have a valid double stencil for all local cells. This is not
needed for the default algorithms but can be useful for code
extensions.

-----

**TETRA_INDEX_IN_FACE**

Adds an extra index to each entry of VF[] and DC[] to one of the
tetrahedra that share this edge. This may be useful for code
extensions.

-----

**NOSTOP_WHEN_BELOW_MINTIMESTEP**

Simulation does not terminate when timestep drops below the specified
minimum timestep size, instead it continues with this timestep
floor. 

-----

**TIMESTEP_OUTPUT_LIMIT**

Limits timesteps such that the requested output times are honored even
if their spacing is finer than the smallest timestep the code makes,
i.e. the code uses the output spacing as an additional timestep
criterion.

-----

**ALLOWEXTRAPARAMS**

Tolerate extra parameters in the parameter file that are not
used. Normally, the code aborts with a complaint if such parameters
are encountered.

-----

**FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES**

This can be used to load SPH ICs that contain particles at identical
coordinates.

-----

**RECOMPUTE_POTENTIAL_IN_SNAPSHOT**

Needed for post-processing option 18 that can be used to calculate potential 
values for a snapshot.

-----

**ACTIVATE_MINIMUM_OPENING_ANGLE**

This does not open tree nodes under the relative opening criterion any
more if their opening angle has dropped below a minimum angle.

-----

**USE_DIRECT_IO_FOR_RESTARTS**

Try to use O_DIRECT for low-level read/write operations of restart
files to circumvent linux kernel page caching.

-----

**HUGEPAGES**

Use huge pages for memory allocation, through hugetlbfs library. Only
possible if the machine supports this.

-----

**DETAILEDTIMINGS**

Creates individual timing entries for primary/secondary kernels to
help in diagnosing work-load balancing.

-----

**BITS_PER_DIMENSION=42**

Bits per dimension used for computing Peano-Hilbert keys. (default: 42)

-----


Input options 
=============

**COMBINETYPES**

Reads in the IC file types 4+5 as type 3.

-----

**LOAD_TYPES=1+2+4+16+32**

Load only specific types sum(2^type).

-----

**READ_COORDINATES_IN_DOUBLE**

Read coordinates in double precision.

-----

**LONGIDS**

If this is set, the code stores particle-IDs as 64-bit long
integers. This is only really needed if you want to go beyond ~2
billion particles.

-----

**OFFSET_FOR_NON_CONTIGUOUS_IDS**

Determines offset of IDs on startup instead of using fixed offset.

-----

**GENERATE_GAS_IN_ICS**

Generates gas from dark matter only ICs (using particle type 1 by default).

-----

**SPLIT_PARTICLE_TYPE=4+8**

Overrides splitting particle type 1 in ``GENERATE_GAS_IN_ICS`` use sum(2^type).

-----

**SHIFT_BY_HALF_BOX**

Shift all positions by half a box size after reading in.

-----

**NTYPES_ICS=6**

Number of particle types in ICs, if not ``NTYPES``.

-----

**READ_MASS_AS_DENSITY_IN_INPUT**

Reads the mass field in the IC as density.

-----

Special input options 
=====================

**IDS_OFFSET=1**

Override offset for gas cell IDs if created from dark matter
particles.

-----

**READ_DM_AS_GAS**

Reads in dark matter particles as gas cells.

-----

**TILE_ICS**

Tile ICs by TileICsFactor (specified as parameter) in each dimension.

-----

Output fields 
==========================

Default output fields are: ``position``, ``velocity``, ``ID``, ``mass``, 
``specific internal energy`` (gas only), ``density`` (gas only)

**OUTPUT_TASK**

Output of MPI rank on which a certainl cell or particle resides.

-----

**OUTPUT_TIMEBIN_HYDRO**

Output of hydrodynamical time-bin.

-----

**OUTPUT_PRESSURE_GRADIENT**

Output of pressure gradient.

-----

**OUTPUT_DENSITY_GRADIENT**

Output of density gradient.

-----

**OUTPUT_VELOCITY_GRADIENT**

Output of velocity gradient.

-----

**OUTPUT_BFIELD_GRADIENT**

Output of magnetic field gradient.

-----

**OUTPUT_MESH_FACE_ANGLE**

Output of maximum face angle of cells.

-----

**OUTPUT_VERTEX_VELOCITY**

Output of velocity of mesh-generating points.

-----

**OUTPUT_VOLUME**

Output of volume of cells; note that this can always be computed from density 
and mass of cells, which are included by default in the output.

-----

**OUTPUT_CENTER_OF_MASS**

Output of center of mass of cells (``Pos`` is the position of the
mesh-generating point).

-----

**OUTPUT_SURFACE_AREA**

Output of surface area of cells as well as the number of faces.

-----

**OUTPUT_PRESSURE**

Output of pressure of gas.

-----

**OUTPUTPOTENTIAL**

This will force the code to compute gravitational potential values for
all particles and cells each time a snapshot file is generated. These
values are then included in the snapshot files. Note that the
computation of the values of the potential costs additional time.

-----

**OUTPUTACCELERATION**

Output of gravitational acceleration.

-----

**OUTPUTTIMESTEP**

Output of timestep of particle.

-----

**OUTPUT_SOFTENINGS**

Output of particle softenings.

-----

**OUTPUTGRAVINTERACTIONS**

Output of gravitational interaction count (from the gravitational
tree) of particles, which is used in the work-load balancing
algorithm.

-----

**OUTPUTCOOLRATE**

Output of cooling rate.

-----

**OUTPUT_DIVVEL**

Output of velocity divergence.

-----

**OUTPUT_CURLVEL**

Output of velocity curl.

-----

**OUTPUT_COOLHEAT**

Output of actual energy loss/gain in cooling/heating routine.

-----

**OUTPUT_VORTICITY**

Output of vorticity of gas.

-----

**OUTPUT_CSND**

Output of sound speed. This field is only used for tree-based timesteps. 
Calculate from hydro quantities in post-processing if required for science 
applications.

-----

Output options 
==============

**PROCESS_TIMES_OF_OUTPUTLIST**

Goes through times of output list prior to starting the simulation to
ensure that outputs are written as close to the desired time as
possible (i.e. also up to half a timestep size before it, as opposed
to always after at the next possible time if this flag is not active).

-----

**REDUCE_FLUSH**

If enabled, files and stdout are only flushed after a certain time
defined in the parameter file (standard behavior: everything is
flushed whenever something is written to it).

-----

**OUTPUT_EVERY_STEP**

Create snapshot on every (global) synchronization point, independent of 
parameters chosen or output list. 

-----

**OUTPUT_CPU_CSV**

Output of a cpu.csv file on top of cpu.txt.

-----

**HAVE_HDF5**

If this is set, the code will be compiled with support for input and
output in the HDF5 format. You need to have the HDF5 libraries and
headers installed on your computer for this option to work. The HDF5
format can then be selected as format "3" in Arepo's parameterfile.

-----

**HDF5_FILTERS**

Activate snapshot compression and checksum for HDF5 output.

-----

**OUTPUT_XDMF**

Writes an ``.xmf`` file for each snapshot, which can be read by the
visualization toolkit Visit (with the hdf5 snapshot). Note: so far
only working if the snapshot is stored in one file.

-----

Testing and Debugging
=============================

**DEBUG**

Enables core-dumps.

-----

**VERBOSE**

Reports readjustments of buffer sizes.

-----

Re-gridding
============================

These options are auxiliary modes to prepare/convert/relax initial
conditions and will not carry out a simulation.

**MESHRELAX**

This keeps the mass constant and only regularizes the mesh.

-----

**ADDBACKGROUNDGRID=16**

Re-grid hydrodynamical quantities on an oct-tree AMR grid. This does
not perform a simulation. This "converts" an SPH initial condition
into a (moving) mesh initial condition.
