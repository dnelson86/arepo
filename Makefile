
# AREPO Makefile
#
# If you add a new system below, also add that systype to Template-Makefile.systype

EXEC   = Arepo
LIBRARY = arepo
CONFIG   = Config.sh
BUILD_DIR = build
SRC_DIR = src

###################
#determine SYSTYPE#
###################
ifdef SYSTYPE
SYSTYPE := "$(SYSTYPE)"
-include Makefile.systype
else
include Makefile.systype
endif

MAKEFILES = Makefile config-makefile
ifeq ($(wildcard Makefile.systype), Makefile.systype)
MAKEFILES += Makefile.systype
endif

$(info Build configuration:)
$(info SYSTYPE: $(SYSTYPE))
$(info CONFIG: $(CONFIG))
$(info EXEC: $(EXEC))
$(info )

PYTHON = python
PERL   = /usr/bin/perl
RESULT     := $(shell CONFIG=$(CONFIG) PERL=$(PERL) BUILD_DIR=$(BUILD_DIR) make -f config-makefile)
CONFIGVARS := $(shell cat $(BUILD_DIR)/arepoconfig.h)
RESULT     := $(shell SRC_DIR=$(SRC_DIR) BUILD_DIR=$(BUILD_DIR) ./git_version.sh)

# Default
MPICH_LIB  = -lmpich
GMP_LIB    = -lgmp
GSL_LIB    = -lgsl -lgslcblas
MATH_LIB   = -lm -lstdc++
HWLOC_LIB = -lhwloc


# e.g. Mac OS using MacPorts modules for openmpi, fftw, gsl, hdf5 and hwloc
ifeq ($(SYSTYPE),"Darwin")
# compiler and its optimization options
CC        =  mpicc   # sets the C-compiler
OPTIMIZE  =  -std=c11 -ggdb -O3 -Wall -Wno-format-security -Wno-unknown-pragmas -Wno-unused-function

# overwrite default:
MPICH_LIB = -lmpi
GSL_INCL  = -I/opt/local/include
GSL_LIB   = -L/opt/local/lib -lgsl -lgslcblas
HWLOC_LIB = -L/opt/local/lib -lhwloc

# libraries that are included on demand, depending on Config.sh options
FFTW_INCL = -I/opt/local/include -I/usr/local/include
FFTW_LIBS = -L/opt/local/lib -I/usr/local/lib
HDF5_INCL = -I/opt/local/include -DH5_USE_16_API 
HDF5_LIB  = -L/opt/local/lib  -lhdf5 -lz
HWLOC_INCL= -I/opt/local/include
endif

# insert the library paths for your system here, similar to SYSTYPE "Darwin" above


ifndef LINKER
LINKER = $(CC)
endif


##########################################
#determine the needed object/header files#
##########################################

SUBDIRS = . \
          debug_md5 \
          domain \
          gitversion \
          gravity \
          gravity/pm \
          hydro \
          init \
          io \
          main \
          mesh \
          mesh/voronoi \
          mpi_utils \
          ngbtree \
          pm \
          star_formation \
          time_integration \
          utils

OBJS =   debug_md5/calc_checksum.o \
         debug_md5/Md5.o \
         domain/domain.o \
         domain/domain_balance.o \
         domain/domain_box.o \
         domain/domain_counttogo.o \
         domain/domain_DC_update.o \
         domain/domain_exchange.o \
         domain/domain_rearrange.o \
         domain/domain_sort_kernels.o \
         domain/domain_toplevel.o \
         domain/domain_vars.o \
         domain/peano.o \
         gravity/accel.o \
         gravity/forcetree.o \
         gravity/forcetree_ewald.o  \
         gravity/forcetree_optimizebalance.o \
         gravity/forcetree_walk.o \
         gravity/grav_external.o \
         gravity/grav_softening.o \
         gravity/gravdirect.o \
         gravity/gravtree.o \
         gravity/gravtree_forcetest.o \
         gravity/longrange.o \
         gravity/pm/pm_periodic2d.o \
         gravity/pm/pm_periodic.o \
         gravity/pm/pm_mpi_fft.o \
         gravity/pm/pm_nonperiodic.o \
         hydro/finite_volume_solver.o \
         hydro/gradients.o \
         hydro/riemann.o \
         hydro/riemann_hllc.o \
         hydro/riemann_hlld.o \
         hydro/scalars.o \
         hydro/update_primitive_variables.o \
         init/begrun.o \
         init/density.o \
         init/init.o \
         io/global.o \
         io/hdf5_util.o \
         io/io.o \
         io/io_fields.o \
         io/logs.o \
         io/parameters.o \
         io/read_ic.o \
         io/restart.o \
         main/allvars.o \
         main/main.o \
         main/run.o \
         mesh/criterion_derefinement.o \
         mesh/criterion_refinement.o \
         mesh/refinement.o \
         mesh/set_vertex_velocities.o \
         mesh/voronoi/voronoi.o \
         mesh/voronoi/voronoi_1d.o \
         mesh/voronoi/voronoi_1d_spherical.o \
         mesh/voronoi/voronoi_3d.o \
         mesh/voronoi/voronoi_check.o \
         mesh/voronoi/voronoi_derefinement.o \
         mesh/voronoi/voronoi_dynamic_update.o \
         mesh/voronoi/voronoi_exchange.o \
         mesh/voronoi/voronoi_ghost_search.o \
         mesh/voronoi/voronoi_gradients_lsf.o \
         mesh/voronoi/voronoi_gradients_onedims.o \
         mesh/voronoi/voronoi_refinement.o \
         mesh/voronoi/voronoi_utils.o \
         mpi_utils/checksummed_sendrecv.o \
         mpi_utils/hypercube_allgatherv.o \
         mpi_utils/mpi_util.o \
         mpi_utils/myalltoall.o \
         mpi_utils/sizelimited_sendrecv.o \
         mpi_utils/pinning.o \
         ngbtree/ngbtree.o \
         ngbtree/ngbtree_search.o \
         ngbtree/ngbtree_walk.o \
         star_formation/sfr_eEOS.o \
         star_formation/starformation.o \
         time_integration/darkenergy.o \
         time_integration/do_gravity_hydro.o \
         time_integration/driftfac.o \
         time_integration/predict.o \
         time_integration/timestep.o \
         time_integration/timestep_treebased.o \
         utils/allocate.o \
         utils/debug.o \
         utils/mpz_extension.o \
         utils/mymalloc.o \
         utils/parallel_sort.o \
         utils/system.o

INCL += debug_md5/Md5.h \
        domain/bsd_tree.h \
        domain/domain.h \
        gitversion/version.h\
        gravity/forcetree.h \
        main/allvars.h \
        main/proto.h \
        mesh/mesh.h \
        mesh/voronoi/voronoi.h \
        time_integration/timestep.h \
        utils/dtypes.h \
        utils/generic_comm_helpers2.h \
        utils/timer.h 

ifeq (TWODIMS,$(findstring TWODIMS,$(CONFIGVARS)))
OBJS    += mesh/voronoi/voronoi_2d.o
endif

ifeq (MYIBARRIER,$(findstring MYIBARRIER,$(CONFIGVARS)))
OBJS    += mpi_utils/myIBarrier.o
INCL    += mpi_utils/myIBarrier.h
endif

ifeq (MHD,$(findstring MHD,$(CONFIGVARS)))
OBJS    += hydro/mhd.o
endif

ifeq (ADDBACKGROUNDGRID,$(findstring ADDBACKGROUNDGRID,$(CONFIGVARS)))
OBJS    += add_backgroundgrid/add_bggrid.o \
           add_backgroundgrid/calc_weights.o \
           add_backgroundgrid/distribute.o
INCL    += add_backgroundgrid/add_bggrid.h
SUBDIRS += add_backgroundgrid
endif

ifeq (COOLING,$(findstring COOLING,$(CONFIGVARS)))
OBJS    += cooling/cooling.o
INCL    += cooling/cooling_vars.h \
           cooling/cooling_proto.h
SUBDIRS += cooling 
endif

ifeq (FOF,$(findstring FOF,$(CONFIGVARS)))
OBJS    += fof/fof.o \
           fof/fof_distribute.o \
           fof/fof_findgroups.o \
           fof/fof_io.o \
           fof/fof_nearest.o \
           fof/fof_sort_kernels.o \
           fof/fof_vars.o
INCL    += fof/fof.h
SUBDIRS += fof
endif

ifeq (SUBFIND,$(findstring SUBFIND,$(CONFIGVARS)))
OBJS	+= subfind/subfind.o \
           subfind/subfind_vars.o \
           subfind/subfind_serial.o \
           subfind/subfind_coll_tree.o \
           subfind/subfind_properties.o \
           subfind/subfind_so.o \
           subfind/subfind_distribute.o \
           subfind/subfind_collective.o \
           subfind/subfind_findlinkngb.o \
           subfind/subfind_nearesttwo.o \
           subfind/subfind_loctree.o \
           subfind/subfind_coll_domain.o \
           subfind/subfind_coll_treewalk.o \
           subfind/subfind_density.o \
           subfind/subfind_io.o \
           subfind/subfind_sort_kernels.o \
           subfind/subfind_reprocess.o \
           subfind/subfind_so_potegy.o
INCL	+= subfind/subfind.h
SUBDIRS += subfind
endif


################################
#determine the needed libraries#
################################

# we only need fftw if PMGRID is turned on, make sure variable is empty otherwise
FFTW_LIB =
ifeq (PMGRID, $(findstring PMGRID, $(CONFIGVARS)))
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
FFTW_LIB = $(FFTW_LIBS) -lfftw3 
else
FFTW_LIB = $(FFTW_LIBS) -lfftw3f 
endif
endif

ifneq (HAVE_HDF5,$(findstring HAVE_HDF5,$(CONFIGVARS)))
HDF5_INCL = 
HDF5_LIB = 
endif

ifneq (IMPOSE_PINNING,$(findstring IMPOSE_PINNING,$(CONFIGVARS)))
HWLOC_INCL = 
HWLOC_LIB = 
endif


##########################
#combine compiler options#
##########################

CFLAGS = $(OPTIMIZE) $(HDF5_INCL) $(GSL_INCL) $(FFTW_INCL) $(HWLOC_INCL) -I$(BUILD_DIR)

LIBS = $(GMP_LIB) $(MATH_LIB) $(MPICH_LIB) $(HDF5_LIB) $(GSL_LIB) $(FFTW_LIB) $(HWLOC_LIB)

FOPTIONS = $(OPTIMIZE)
FFLAGS = $(FOPTIONS)


SUBDIRS := $(addprefix $(BUILD_DIR)/,$(SUBDIRS))
OBJS := $(addprefix $(BUILD_DIR)/,$(OBJS)) $(BUILD_DIR)/compile_time_info.o $(BUILD_DIR)/compile_time_info_hdf5.o $(BUILD_DIR)/version.o
INCL := $(addprefix $(SRC_DIR)/,$(INCL)) $(BUILD_DIR)/arepoconfig.h

TO_CHECK := $(addsuffix .check, $(OBJS) $(patsubst $(SRC_DIR)%, $(BUILD_DIR)%, $(INCL)) )
TO_CHECK +=  $(BUILD_DIR)/Makefile.check
CONFIG_CHECK = $(BUILD_DIR)/$(notdir $(CONFIG)).check

DOCS_CHECK = $(BUILD_DIR)/README.check

################
#create subdirs#
################
RESULT := $(shell mkdir -p $(SUBDIRS)  )


#############
#build rules#
#############

all: check build 

build: $(EXEC)

$(EXEC): $(OBJS)
	$(LINKER) $(OPTIMIZE) $(OBJS) $(LIBS) -o $(EXEC)

lib$(LIBRARY).a: $(filter-out $(BUILD_DIR)/main.o,$(OBJS))
	$(AR) -rcs lib$(LIBRARY).a $(OBJS)

clean:
	@echo Cleaning all build files...
	@rm -f $(OBJS) $(EXEC) lib$(LIBRARY).a
	@rm -f $(BUILD_DIR)/compile_time_info.c $(BUILD_DIR)/compile_time_info_hdf5.c $(BUILD_DIR)/arepoconfig.h
	@rm -f $(BUILD_DIR)/version.c
	@rm -f $(TO_CHECK) $(CONFIG_CHECK)
	@rm -rf $(BUILD_DIR)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c $(INCL) $(MAKEFILES)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/compile_time_info.o: $(BUILD_DIR)/compile_time_info.c $(MAKEFILES)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/compile_time_info_hdf5.o: $(BUILD_DIR)/compile_time_info_hdf5.c $(MAKEFILES)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cu $(INCL) $(MAKEFILES)
	$(NVCC)  -c $< -o $@

# sanity checks:

check: $(CONFIG_CHECK)

check_docs: $(DOCS_CHECK)

$(CONFIG_CHECK): $(TO_CHECK) $(CONFIG) check.py 
	@$(PYTHON) check.py 2 $(CONFIG) $(CONFIG_CHECK) defines_extra $(TO_CHECK)

$(BUILD_DIR)/%.o.check: $(SRC_DIR)/%.c Template-Config.sh defines_extra check.py
	@$(PYTHON) check.py 1 $< $@ Template-Config.sh defines_extra

$(BUILD_DIR)/%.o.check: $(SRC_DIR)/%.F
	touch $@

$(BUILD_DIR)/%.o.check: $(SRC_DIR)/%.f90
	touch $@

$(BUILD_DIR)/%.o.check: $(SRC_DIR)/%.F90
	touch $@

$(BUILD_DIR)/%.o.check: $(SRC_DIR)/%.cc
	touch $@

$(BUILD_DIR)/%.h.check: $(SRC_DIR)/%.h Template-Config.sh defines_extra check.py
	@$(PYTHON) check.py 1 $< $@ Template-Config.sh defines_extra

$(BUILD_DIR)/%.o.check: $(BUILD_DIR)/%.c Template-Config.sh defines_extra check.py
	@$(PYTHON) check.py 1 $< $@ Template-Config.sh defines_extra

$(BUILD_DIR)/%.h.check: $(BUILD_DIR)/%.h Template-Config.sh defines_extra check.py
	@$(PYTHON) check.py 1 $< $@ Template-Config.sh defines_extra

$(BUILD_DIR)/Makefile.check: Makefile Template-Config.sh defines_extra check.py
	@$(PYTHON) check.py 3 $< $@ Template-Config.sh defines_extra

$(BUILD_DIR)/Config.check: Template-Config.sh check.py
	@$(PYTHON) check.py 4 Template-Config.sh $@

