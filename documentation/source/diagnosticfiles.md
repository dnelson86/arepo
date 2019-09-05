*****************
Diagnostic output
*****************

Arepo will not only output the simulation snapshot and reduced data
via the halo-finder files, but also a number of (mostly ascii)
diagnostic log- files which contain important information about the
code performance and runtime behavior.

In practice, to quickly check the performance of large production
runs, it is useful to check the ``timebins.txt`` and ``cpu.txt``
files. The former will give information how many simulation elements
are evolved with which timestep sizes, i.e. characteristics of the
system simulated, the latter provides information about the
computational time spent in each part of the code, which can be
influenced to some degree by the values of the code parameters.

For ongoing simulations, these can be checked via

.. code-block :: python

    tail -n 30 timebins.txt
    tail -n 200 cpu.txt


stdout 
======

The standard output contains general information about the simulation
status and many of the main routines will print general information in
it.  The output itself is mainly relevant for reconstructing what the
simulation did, which is needed e.g. for debugging purposes.

balance.txt 
===========

Output of fractional cpu time used in each individual step, optimized
to be machine readable (while cpu.txt is more human readable).

Symbol key

.. code-block :: python

    total                = '-' / '-'
    treegrav             = 'a' / ')'
    treebuild            = 'b' / '('
    insert               = 'c' / '*'
    branches             = 'd' / '&'
    toplevel             = 'e' / '^'
    treecostm            = 'f' / '%'
    treewalk             = 'g' / '$'
    treewalk1            = 'h' / '#'
    treewalk2            = 'i' / '@'
    treebalsndrcv        = 'j' / '!'
    treeback             = 'm' / '7'
    treedirect           = 'r' / '2'
    pm_grav              = 's' / '1'
    ngbtreebuild         = 't' / 'Z'
    ngbtreevelupdate     = 'u' / 'Y'
    voronoi              = 'v' / 'X'
    insert               = 'w' / 'W'
    findpoints           = 'x' / 'V'
    cellcheck            = 'y' / 'U'
    geometry             = 'z' / 'T'
    exchange             = 'A' / 'S'
    dynamic              = 'B' / 'R'
    hydro                = 'C' / 'Q'
    gradients            = 'D' / 'P'
    fluxes               = 'F' / 'N'
    fluxcomm             = 'H' / 'L'
    updates              = 'J' / 'j'
    vertex vel           = 'K' / 'I'
    mhd                  = '4' / 'p'
    domain               = 'U' / 'y'
    peano                = 'V' / 'x'
    drift/kicks          = 'W' / 'w'
    timeline             = 'X' / 'v'
    treetimesteps        = 'Y' / 'u'
    i/o                  = 'Z' / 't'
    logs                 = '1' / 's'
    sfrcool              = '2' / 'r'
    fof                  = '#' / 'h'
    subfind              = '$' / 'g'
    refine               = '%' / 'f'
    mesh_derefine        = '^' / 'e'
    images               = '&' / 'd'
    initializ.           = '*' / 'c'
    restart              = '(' / 'b'
    misc                 = ')' / 'a'

example:

.. code-block :: python

    Step=  13147  sec=     0.212 Nsync-grv=      5181 Nsync-hyd=       
    342  dhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhttwwwwwxxzBBBBBBBBBBBBBBB
    BBBBBBBBBBBBBBBBBBBBBBDFFFUUUUUUUUUUUUUUUUUUVVVVWXY2222^^
    Step=  13148  sec=     0.001 Nsync-grv=         1 Nsync-hyd=       
      1  rvvwwwwwwwwwwwwwwwxxxxxxxxxxxxxyyzzzzzABBBDDFFFFFFFFHVVVVWWWYYY
    Y1111111111111111111111111111112222222222222^))))))))))))


cpu.txt 
=======

For each sync-point, such a block is written. This file reports
measurements of the different timers built into Arepo. Each
computationally expensive operation has a different timer attached to
it, thus allowing to closely monitor where the computational time is
spent.  Some of the timers (e.g. treegrav) have sub-timers for
individual operations.  This is denoted by the indentation hierarchy
in the first column.  The fraction of time spent in different code
parts, as well as the absolute amount, is highly problem
dependent. The timers make it possible to identify inefficient parts
of the overall algorithm and concentrate on the most time-consuming
parts of the code. There is also the option ``OUTPUT_CPU_CSV`` which
returns thes same data as a more easily machine-readable ``cpu.csv``
file.

The different columns are: name; wallclock time (in s) this step;
percentage this step; wallclock time (in s) cumulative; percentage up
to this step. A typical block of cpu.txt looks the following (here a
gravity-only, tree-only run):

.. code-block :: python

    Step 131, Time: 0.197266, CPUs: 1, MultiDomains: 8, HighestActiveTim
    eBin: 20
                              diff               cumulative
    total                     1.03  100.0%      99.22  100.0%
      treegrav                1.00   97.6%      96.12        96.9%
        treebuild             0.01    0.7%       0.76         0.8%
          insert              0.00    0.5%       0.54         0.5%
          branches            0.00    0.2%       0.21         0.2%
          toplevel            0.00    0.0%       0.00         0.0%
        treecostm             0.00    0.3%       0.30         0.3%
        treewalk              0.99   96.5%      94.96        95.7%
          treewalk1           0.99   96.5%      94.96        95.7%
          treewalk2           0.00    0.0%       0.00         0.0%
        treebalsndrcv         0.00    0.0%       0.00         0.0%
        treeback              0.00    0.0%       0.00         0.0%
        treedirect            0.00    0.0%       0.01         0.0%
      ngbtreebuild            0.00    0.0%       0.00         0.0%
      ngbtreevelupdate        0.00    0.0%       0.01         0.0%
      voronoi                 0.00    0.0%       0.00         0.0%
        insert                0.00    0.0%       0.00         0.0%
        findpoints            0.00    0.0%       0.00         0.0%
        cellcheck             0.00    0.0%       0.00         0.0%
        geometry              0.00    0.0%       0.00         0.0%
        exchange              0.00    0.0%       0.00         0.0%
        dynamic               0.00    0.0%       0.00         0.0%
      hydro                   0.00    0.0%       0.01         0.0%
        gradients             0.00    0.0%       0.00         0.0%
        fluxes                0.00    0.0%       0.00         0.0%
        fluxcomm              0.00    0.0%       0.00         0.0%
        updates               0.00    0.0%       0.00         0.0%
        vertex vel            0.00    0.0%       0.00         0.0%
        mhd                   0.00    0.0%       0.00         0.0%
      domain                  0.01    1.4%       1.65         1.7%
      peano                   0.01    0.5%       0.59         0.6%
      drift/kicks             0.00    0.4%       0.74         0.7%
      timeline                0.00    0.0%       0.03         0.0%
      treetimesteps           0.00    0.0%       0.00         0.0%
      i/o                     0.00    0.0%       0.03         0.0%
      logs                    0.00    0.0%       0.03         0.0%
      sfrcool                 0.00    0.0%       0.00         0.0%
      refine                  0.00    0.0%       0.00         0.0%
      mesh_derefine           0.00    0.0%       0.00         0.0%
      images                  0.00    0.0%       0.00         0.0%
      initializ.              0.00    0.0%       0.00         0.0%
      restart                 0.00    0.0%       0.00         0.0%
      misc                    0.00    0.0%       0.02         0.0%


domain.txt
==========

The load-balancing (cpu work and memory) both for gravity and hydro
calculations is reported for each timebin individually. Reported every
sync-point.  Ideally balanced runs have a value 1, the higher the
value, the more imbalanced the simulation.

.. code-block :: python

    DOMAIN BALANCE, Sync-Point 13314, Time: 0.997486
    Timebins:       Gravity       Hydro  cumulative      grav-balance   
        hydro-balance
     |bin=18          47268       26570       64920  m  1.000 | 1.000   
        1.000 | 1.000
     |bin=16          12523        2399       17652  m  1.000 | 1.000   
        1.000 | 1.000
    >|bin=14           5129         331        5129  m  1.000 | 1.000  *
        1.000 | 1.000
    --------------------------------------------------------------------
    -----------------
    BALANCE,  LOAD:   1.000       1.000       1.000  WORK:      1.000   
                1.000
    --------------------------------------------------------------------
    -----------------

energy.txt
==========

In specified intervals (in simulation time, specified by the parameter
`TimeBetStatistics`) the total energy and its components are computed
and written into the file `energy.txt`. This file also contains the
cumulative energy that had to be injected into the system to ensure
positivity in thermal energy.  All output is in code units. Note: this
only works with up to 6 particle types.  The columns are

.. code-block :: python

    1. simulation time/ scalefactor 
    2. total thermal energy
    3. total potential energy
    4. total kinetic energy
    5. internal energy particle type 0
    6. potential energy particle type 0
    7. kinetic energy particle type 0
    8. internal energy particle type 1
    9. potential energy particle type 1
    10. kinetic energy particle type 1
    11. internal energy particle type 2
    12. potential energy particle type 2
    13. kinetic energy particle type 2
    14. internal energy particle type 3
    15. potential energy particle type 3
    16. kinetic energy particle type 3
    17. internal energy particle type 4
    18. potential energy particle type 4
    19. kinetic energy particle type 4
    20. internal energy particle type 5
    21. potential energy particle type 5
    22. kinetic energy particle type 5
    23. total mass in particle type 0
    24. total mass in particle type 1
    25. total mass in particle type 2
    26. total mass in particle type 3
    27. total mass in particle type 4
    28. total mass in particle type 5
    29. total injected energy due to positivity enforcement of thermal e
    nergy

Two example lines:

.. code-block :: python

    0.96967 3.29069e+06 0 4.27406e+07 3.29069e+06 0 1.65766e+06 0 0 3.93
    02e+07 0 0 0 0 0 0 0 0 1.78097e+06 0 0 0 503.489 3047.89 0 0 65.5756
     0 7.71477
    0.978903 3.28666e+06 0 4.18631e+07 3.28666e+06 0 1.67505e+06 0 0 3.8
    4203e+07 0 0 0 0 0 0 0 0 1.76774e+06 0 0 0 503.306 3047.89 0 0 65.75
    86 0 7.71477

info.txt
========

Every sync-point, the time-bins, time, timestep and number of active
particles are written into this file, e.g.

.. code-block :: python

    Sync-Point 13327, TimeBin=16, Time: 0.999408, Redshift: 0.000592464,
     Systemstep: 0.000147974, Dloga: 0.000148072, Nsync-grv:      17679,
     Nsync-hyd:       2739

memory.txt 
==========

Arepo uses its own internal memory manager. This means that one large
chunk of memory is reserved initially for Arepo (specified by the
parameter `MaxMemSize`), and the allocation for individual arrays is
then handled internally from this pool.  The reason for introducing
this was to avoid memory fragmentation during runtime on some
machines, but also to have detailed information about how much memory
Arepo actually needs and to terminate if this exceeds a pre-defined
threshold. ``memory.txt`` reports this internal memory usage, and how
much memory is actually needed by the simulation.

.. code-block :: python

    MEMORY:  Largest Allocation = 816.742 Mbyte  |  Largest Allocation W
    ithout Generic = 132.938 Mbyte
    
    -------------------------- Allocated Memory Blocks---- ( Step      
     0 )------------------
    Task    Nr F                  Variable      MBytes   Cumulative  Fun
    ction|File|Linenumber
    --------------------------------------------------------------------
    ----------------------
       0     0 0                               Parameters      0.0040   
        0.0040  read_parameter_file()|src/io/parameters.c|648
       0     1 0                          ParametersValue      0.0157   
        0.0197  read_parameter_file()|src/io/parameters.c|649
       0     2 0                            ParamtersType      0.0001   
        0.0198  read_parameter_file()|src/io/parameters.c|650
       0     3 0                                IO_Fields      0.0052   
        0.0250  init_field()|src/io/io.c|130
       0     4 0                                    RateT      0.1527   
        0.1777  InitCool()|src/cooling/cooling.c|794
       0     5 0                                   PhotoT      0.0058   
        0.1835  ReadIonizeParams()|src/cooling/cooling.c|691
       0     6 0                             slab_to_task      0.0005   
        0.1840  my_slab_based_fft_init()|src/gravity/pm/pm_mpi_fft.c|101
       0     7 0                         slabs_x_per_task      0.0001   
        0.1840  my_slab_based_fft_init()|src/gravity/pm/pm_mpi_fft.c|116
       0     8 0                     first_slab_x_of_task      0.0001   
        0.1841  my_slab_based_fft_init()|src/gravity/pm/pm_mpi_fft.c|119
       0     9 0                         slabs_y_per_task      0.0001   
        0.1841  my_slab_based_fft_init()|src/gravity/pm/pm_mpi_fft.c|122
       0    10 0                     first_slab_y_of_task      0.0001   
        0.1842  my_slab_based_fft_init()|src/gravity/pm/pm_mpi_fft.c|125
       0    11 0                               Exportflag      0.0001   
        0.1843  allocate_memory()|src/utils/allocate.c|46
       0    12 0                              Exportindex      0.0001   
        0.1843  allocate_memory()|src/utils/allocate.c|47
       0    13 0                          Exportnodecount      0.0001   
        0.1844  allocate_memory()|src/utils/allocate.c|48
       0    14 0                                     Send      0.0001   
        0.1844  allocate_memory()|src/utils/allocate.c|50
       0    15 0                                     Recv      0.0001   
        0.1845  allocate_memory()|src/utils/allocate.c|51
       0    16 0                            TasksThatSend      0.0001   
        0.1846  allocate_memory()|src/utils/allocate.c|53
       0    17 0                            TasksThatRecv      0.0001   
        0.1846  allocate_memory()|src/utils/allocate.c|54
       0    18 0                               Send_count      0.0001   
        0.1847  allocate_memory()|src/utils/allocate.c|56
       0    19 0                              Send_offset      0.0001   
        0.1848  allocate_memory()|src/utils/allocate.c|57
       0    20 0                               Recv_count      0.0001   
        0.1848  allocate_memory()|src/utils/allocate.c|58
       0    21 0                              Recv_offset      0.0001   
        0.1849  allocate_memory()|src/utils/allocate.c|59
       0    22 0                         Send_count_nodes      0.0001   
        0.1849  allocate_memory()|src/utils/allocate.c|61
       0    23 0                        Send_offset_nodes      0.0001   
        0.1850  allocate_memory()|src/utils/allocate.c|62
       0    24 0                         Recv_count_nodes      0.0001   
        0.1851  allocate_memory()|src/utils/allocate.c|63
       0    25 0                        Recv_offset_nodes      0.0001   
        0.1851  allocate_memory()|src/utils/allocate.c|64
       0    26 0                          Mesh_Send_count      0.0001   
        0.1852  allocate_memory()|src/utils/allocate.c|66
       0    27 0                         Mesh_Send_offset      0.0001   
        0.1852  allocate_memory()|src/utils/allocate.c|67
       0    28 0                          Mesh_Recv_count      0.0001   
        0.1853  allocate_memory()|src/utils/allocate.c|68
       0    29 0                         Mesh_Recv_offset      0.0001   
        0.1854  allocate_memory()|src/utils/allocate.c|69
       0    30 0                         Force_Send_count      0.0001   
        0.1854  allocate_memory()|src/utils/allocate.c|71
       0    31 0                        Force_Send_offset      0.0001   
        0.1855  allocate_memory()|src/utils/allocate.c|72
       0    32 0                         Force_Recv_count      0.0001   
        0.1855  allocate_memory()|src/utils/allocate.c|73
       0    33 0                        Force_Recv_offset      0.0001   
        0.1856  allocate_memory()|src/utils/allocate.c|74
       0    34 1                                        P     11.2500   
       11.4356  allocate_memory()|src/utils/allocate.c|77
       0    35 1                                     SphP     15.3125   
       26.7481  allocate_memory()|src/utils/allocate.c|80
       0    36 1                  NextActiveParticleHydro      0.1562   
       26.9044  timebins_allocate()|src/time_integration/timestep.c|686
       0    37 1                       NextInTimeBinHydro      0.1562   
       27.0606  timebins_allocate()|src/time_integration/timestep.c|689
       0    38 1                       PrevInTimeBinHydro      0.1562   
       27.2169  timebins_allocate()|src/time_integration/timestep.c|692
       0    39 1                NextActiveParticleGravity      0.3125   
       27.5294  timebins_allocate()|src/time_integration/timestep.c|686
       0    40 1                     NextInTimeBinGravity      0.3125   
       27.8419  timebins_allocate()|src/time_integration/timestep.c|689
       0    41 1                     PrevInTimeBinGravity      0.3125   
       28.1544  timebins_allocate()|src/time_integration/timestep.c|692
       0    42 1                                       DC     14.0267   
       42.1811  voronoi_init_connectivity()|src/mesh/voronoi/voronoi_dyn
       amic_update|709
       0    43 1                          DomainStartList      0.0001   
       42.1812  domain_allocate()|src/domain/domain.c|386
       0    44 1                            DomainEndList      0.0001   
       42.1812  domain_allocate()|src/domain/domain.c|387
       0    45 1                   DomainFirstLocTopleave      0.0001   
       42.1813  domain_allocate()|src/domain/domain.c|388
       0    46 1                     DomainNLocalTopleave      0.0001   
       42.1813  domain_allocate()|src/domain/domain.c|389
       0    47 1                                 TopNodes      0.0020   
       42.1833  domain_allocate()|src/domain/domain.c|390
       0    48 1                               DomainTask      0.0001   
       42.1835  domain_allocate()|src/domain/domain.c|391
       0    49 1               DomainListOfLocalTopleaves      0.0001   
       42.1836  domain_Decomposition()|src/domain/domain.c|144
       0    50 1                      Ngb_DomainNodeIndex      0.0001   
       42.1837  ngb_treeallocate()|src/ngbtree/ngbtree.c|1346
       0    51 1                                Ngb_Nodes      3.0735   
       45.2573  ngb_treeallocate()|src/ngbtree/ngbtree.c|1348
       0    52 1                             ExtNgb_Nodes      0.7684   
       46.0257  ngb_treeallocate()|src/ngbtree/ngbtree.c|1352
       0    53 1                             Ngb_Nextnode      0.1564   
       46.1821  ngb_treeallocate()|src/ngbtree/ngbtree.c|1355
       0    54 1                               Ngb_Father      0.1562   
       46.3383  ngb_treeallocate()|src/ngbtree/ngbtree.c|1356
       0    55 1                               Ngb_Marker      0.2661   
       46.6044  ngb_treeallocate()|src/ngbtree/ngbtree.c|1358
       0    56 1                                       VF     16.0329   
       62.6373  initialize_and_create_first_tetra()|src/mesh/voronoi/vor
    onoi_3d.c|143
       0    57 1                                       DP      4.4672   
       67.1045  initialize_and_create_first_tetra()|src/mesh/voronoi/vor
    onoi_3d.c|145
       0    58 1                                       DT     17.6054   
       84.7099  initialize_and_create_first_tetra()|src/mesh/voronoi/vor
    onoi_3d.c|148
       0    59 1                              ListExports      0.6153   
       85.3252  create_mesh()|src/mesh/voronoi/voronoi.c|169
       0    60 1                              List_InMesh      0.1250   
       85.4502  create_mesh()|src/mesh/voronoi/voronoi.c|172
       0    61 1                                   List_P      0.2500   
       85.7002  create_mesh()|src/mesh/voronoi/voronoi.c|174
       0    62 1                                      DTC     11.7369   
       97.4371  create_mesh()|src/mesh/voronoi/voronoi.c|176
       0    63 1                                 PrimExch      0.0001   
       97.4372  mesh_setup_exchange()|src/mesh/voronoi/voronoi_exchange.
    c|193
       0    64 1                                 GradExch      0.0001   
       97.4373  mesh_setup_exchange()|src/mesh/voronoi/voronoi_exchange.
    c|194
       0    65 0                              Sndpm_count      0.0001   
       97.4373  pmforce_uniform_optimized_prepare_densi()|src/gravity/pm
    /pm_periodic.c|733
       0    66 0                             Sndpm_offset      0.0001   
       97.4374  pmforce_uniform_optimized_prepare_densi()|src/gravity/pm
    /pm_periodic.c|734
       0    67 0                              Rcvpm_count      0.0001   
       97.4374  pmforce_uniform_optimized_prepare_densi()|src/gravity/pm
    /pm_periodic.c|735
       0    68 0                             Rcvpm_offset      0.0001   
       97.4375  pmforce_uniform_optimized_prepare_densi()|src/gravity/pm
    /pm_periodic.c|736
       0    69 0                                   partin      2.0000   
       99.4375  pmforce_uniform_optimized_prepare_densi()|src/gravity/pm
    /pm_periodic.c|880
       0    70 0                                  rhogrid     16.2500   
      115.6875  pmforce_uniform_optimized_prepare_densi()|src/gravity/pm
    /pm_periodic.c|1038
       0    71 0                                forcegrid     16.2500   
      131.9375  pmforce_periodic()|src/gravity/pm/pm_periodic.c|1597
       0    72 0                                  flistin      0.5000   
      132.4375  pmforce_uniform_optimized_readout_force()|src/gravity/pm
    /pm_periodic.c|1328
       0    73 0                                 flistout      0.5000   
      132.9375  pmforce_uniform_optimized_readout_force()|src/gravity/pm
    /pm_periodic.c|1329
       0    74 0                                   scount      0.0001   
      132.9376  myMPI_Alltoallv()|src/mpi_utils/myalltoall.c|67
       0    75 0                                   rcount      0.0001   
      132.9376  myMPI_Alltoallv()|src/mpi_utils/myalltoall.c|68
       0    76 0                                     soff      0.0001   
      132.9377  myMPI_Alltoallv()|src/mpi_utils/myalltoall.c|69
       0    77 0                                     roff      0.0001   
      132.9377  myMPI_Alltoallv()|src/mpi_utils/myalltoall.c|70
    --------------------------------------------------------------------
    ----------------------


sfr.txt 
=======

In case ``USE_SFR`` is active, Arepo will create a ``sfr.txt`` file,
which reports the stars created in every call of the star-formation
routine.

The individual columns are:

 * time (code units or scale factor)
 * total stellar mass to be formed in timestepo prior to stochastic
   sampling (code units),
 * instantaneous star formation rate of all cells (Msun/yr), 
 * instantaneous star formation rate of active cells (Msun/yr), 
 * total mass in stars formed in this timestep (after sampling) (code units), 
 * cumulative stellar mass formed (code units).

Example: 

.. code-block :: python

      4.373019e-01   9.714635e-03   1.100743e+02   1.405136e+02   2.2809
      41e-02   2.752464e+01
      4.373667e-01   4.007455e-04   1.104648e+02   5.795346e+00   0.0000
      00e+00   2.752464e+01
      4.374315e-01   2.009357e-02   1.104276e+02   2.905270e+02   0.0000
      00e+00   2.752464e+01
      4.374962e-01   3.904148e-04   1.103389e+02   5.643836e+00   0.0000
      00e+00   2.752464e+01


timebins.txt 
============

Arepo is optimized for time-integrating both hydrodynamical as well as
gravitational interactions on the largest possible timestep that is
allowed by the timestep criterion and allowed by the binary hierarchy
of time steps.  For each timestep, a linked list of particles on this
particular integration step exists, and their statistics are reported
in `timebins.txt`.  In this file, the number of gas cells and
collisionless particles in each timebin (i.e. integration timestep) is
reported for each sync-point, as well as the cpu time and the fraction
of the total cost contributed by each timebin. A typical block looks
like

.. code-block :: python

    Sync-Point 2658, Time: 0.307419, Redshift: 2.25289, Systemstep: 9.10
    27e-05, Dloga: 0.000296144
    Occupied timebins: gravity      hydro               dt              
    cumul-grav   cumul-sph A D    avg-time  cpu-frac
        bin=19           57159       30466        0.004738310805        
         65419       32200            3.04    100.0%
        bin=17            6854        1222        0.001184577701        
          8260        1734            0.36      1.0%
        bin=16            1221         369        0.000592288851        
          1406         512            0.19      0.0%
     X  bin=15             185         143        0.000296144425        
           185         143 <          0.02      0.3%
                   ------------------------
    Total active:          185         143  


timings.txt
===========

The performance of the gravitational tree algorithm is reported in
`timings.txt` for each sync-point. An example of a single sync-point
looks like the following

.. code-block :: python

    Step(*): 372, t: 0.0455302, dt: 0.000215226, highest active timebin:
     19  (lowest active: 19, highest occupied: 19)
    Nf=    65536  timebin=19  total-Nf=24510464
       work-load balance: 1 (1 0), rel1to2: 1
       number of iterations:  max=1, exported fraction: 0
       part/sec: raw=127058, effective=127058     ia/part: avg=32.2464
       maximum number of nodes: 31091, filled: 0.677246
       avg times: all=0.519064  tree1=0.515797  tree2=0  commwait=1.0013
    6e-05 sec

