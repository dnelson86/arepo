*******************
Getting started
*******************


Library requirements
====================

Arepo requires the following libraries for compilation:

Always needed:

1. **mpi** The `Message Passing Interface`, version 2.0 or higher.  In
   practice Arepo is well-tested on the open source libraries MPICH
   <https://www.mpich.org> and OpenMPI <http://www.open-mpi.org>.  For
   some machines, vendor supplied versions may exist and can also be
   used.  Arepo relies exclusively on MPI for its parallelization.
2. **gsl** `The GNU Scientific Library`. This open-source library is available 
   under <http://www.gnu.org/software/gsl>.
   
Needed in some cases:

3. **fftw3** `Fastest Fourier Transform in the West`
   (<http://www.fftw.org>) This free software package for
   fast-fourier-transformations is used in the particle-mesh
   algorithm. Consequently, it is only needed if ``PMGRID`` is set in
   ``Config.sh``.
4. **hdf5** `Hierarchical Data Format`
   <https://support.hdfgroup.org/HDF5/doc/H5.intro.html> is a data
   format specification that is used for the input/output format 3 in
   Arepo, which is the de-facto standard I/O format (even though the
   old Gadget formats are supported as well). The library is needed if
   the ``HAVE_HDF5`` flag is set in `Config.sh`.
5. **hwloc** : The `hardware locality library` is useful for allowing
   the code to probe the processor topology it is running on and
   enable a pinning of MPI threads to individual cores. This library
   is optional and only required if the ``IMPOSE_PINNING`` flag is
   set.


Arepo needs a working C compiler supporting at least the ``c11`` standard, as 
well as GNU-Make and Python, which are used in the build process.

Building the code 
=================

Systype variable
----------------

The first thing after obtaining a local copy of the code is to set a
``SYSTYPE`` variable, either by typing::

    export SYSYTPE="MySystype"

(assuming you run the bash-shell) or by copying the
``Template-Makefile.systype`` file::

    cp Template-Makefile.systype Makefile.systype

and uncommenting the intended ``SYSTYPE`` variable in
``Makefile.systype``.  The ``SYSTYPE`` is a variable defining a
specific machine. Each machine has its specific compiler settings and
library paths defined in the Makefile. New ``SYSTYPE`` keywords can be
added in the Makefile when needed.

Config.sh -- Compile time flags
-------------------------------

The next step is to set the compile-time configuration of the
code. This is done best by copying the ``Template-Config.sh`` file::

    cp Template-Config.sh Config.sh

Note that the compile-time configuration file does not need to be
called ``Config.sh``, however, if it is not, this needs to be
specified when calling the Makefile.  For example, for the code
examples this file is located in a different location, and they thus
make use of this possibility.

The ``Config.sh`` file enables/disables possible modules of Arepo,
such as, e.g. ``SELFGRAVITY``, which switches on gravitational
interactions between cells/particles. All flags not out-commented with
a '#' will be processed and cast into a ``#define MyModule`` in an
automatically produced header file ``arepoconfig.h`` that is placed
into the code build directory.

Compiling the code
------------------

Arepo can be compiled with::

    make

which will first execute a check on the compile time flags (see below)
followed by the actual build. It is also possible to skip the check by
typing::

    make build

The compilation can be accelerated on multi-core machines by compiling
files in parallel with an option like ``-j 8``, which will use 8 threads
for the compilation.

Compile-time checks of flags
----------------------------

The normal ``make`` command will first execute a check using the
python script ``check.py``. This ensures that all flags used in the
code via ``#if defined(MyFlag)`` or ``#ifdef MyFlag`` occur either in
``Template-Config.sh`` or ``defines_extra``. The main reason for doing
this is to prevent typos in the code such as ``#if defined(MyFlg)``,
which are generally very hard to find because the code might compile
normally, but the code followed by ``#if defined(MyFlg)`` just does
not appear in the executable.

New flags that can be activated should be included to
``Template-Config.sh``, while internally defined flags (such as header
file guards) are listed in ``defines_extra`` in order to exempt them
from the check.

Starting the code
=================

param.txt -- Runtime settings
-----------------------------

After compiling Arepo, the executable ``Arepo`` should be present in
the main directory. To start a simulation, Arepo needs be handed over
the path to a file containing the runtime parameters, usually called
``param.txt``.  The settings that need to be set generally depend on
the compile-time flags used. All possible parameters are listed under
:ref:`Parameterfile`, including the compile-time flags they belong
to. Starting Arepo with a parameter file that has too few or obsolete
parameter settings will result in an output telling which parameters
are required or still needed, respectively.

Starting the simulation
-----------------------

A simulation is started with a command like::

    mpirun -np 32 Arepo param.txt

which will start Arepo on 32 MPI ranks using the ``param.txt`` file
for the run-time parameter settings. Note that on some systems, mpirun
may instead be called mpiexec, or be another launch-command such as srun,
with a slightly different syntax. We refer to the user guide of individual
machines for such peculiarities.

Interrupting a run
==================

Arepo supports running a simulations in several consecutive
installments, and resuming an already started simulation from
so-called restart files (a mechanism sometimes called checkpointing),
which are essentially memory dumps written to the output
directory. The code is written in such a way that the results are
completely independent of the number of restarts performed, which is a
crucial aspect for debugging the code if necessary, as this ensures
deterministic reproducability.

A (regular) code termination always results in the writing of such
restart files, and is normally triggered when more than 85% of the
specified maximum runtime has been reached, or when the final
simulation time has been reached.  It is important, in particular in
larger simulations, that the code is able to complete the output of
the restart file before, e.g. a job time-limit is reached. Therefore
the specified maximum runtime should never exceed the runtime limits
of the machine.

It is also possible to trigger a (regular) code interruption by
introducing a file called ``stop`` in the output directory.::

    touch stop

Furthermore, Arepo usually produces restart-files in regular, user
specified time intervals to ensure runs can be resumed without a major
loss of computation time even if a run was terminated abruptly by the
operating system.

Restarting a run
================

A simulation can be restarted in two ways in Arepo, which is indicated with the
additional flags 1 or 2 in the execution call. The flag 1 is the restarting
from a restart file and generally recommended if possible. This, however, 
has the restriction that a run needs to be resumed with **exactly as many MPI 
ranks** as it was originally started from. The restart-command then looks 
as follows::

    mpirun -np 32 ./Arepo param.txt 1

Alternatively, a flag value of 2 signals that one wants to resume a
simulation from a regular snapshot output. This is possible but
discouraged because here the simulation results will not be binary
identical to a run without this restart. In particular in runs where
some fields (especially positions) are written in single precision,
this can also have effects on the overall quality of the simulation. A
command to restart Arepo from a snapshot looks like the following::

    mpirun -np 32 ./Arepo param.txt 2 7

Here Arepo would read in the snapshot number 7 and resume the
simulation from its corresponding time. This second method is only
recommened if the regular start from restart files is not possible for
some reason.

Starting post-processing tasks
==============================

There are a number of post-processing options in the public version of 
Arepo. The main ones are 

* to run the halo and/or subhalo finder in post-processing (flag 3)
* to convert a snapshot to another format (flag 6),
* to write the Voronoi-mesh data (flag 14)
* to calculate and output the gradients of hydrodynamical quantities (flag 17)
* to calculate and output the gravitational potential (flag 18)

Note that the non-continuity of the numbers originates from the goal
to maintain consistency with the developer version, which has more
options. Other flags are not supported in the current public version
of Arepo.
