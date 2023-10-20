
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2023  The SHYFEM team
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

DIR     = shyelab

#--------------------------------------------- just for default
FEMDIR = ../..
FEMBIN = $(FEMDIR)/bin
LIBDIR = $(FEMDIR)/lib
SRCDIR = $(FEMDIR)/src
F77    = g77
FFLAGS =  -O -Wall -g
LINKER = $(F77)
LFLAGS:=
INCLUDE = -I$(FEMDIR)/include/
#---------------------------------------------

include $(FEMDIR)/Rules.make

#---------------------------------------------

LIBUTIL         = util
LIBSHYUTIL      = shyutil
LIBSHYMPI       = shympi
LIBSHYUTILMPI   = shyutilmpi

LIBLUT  = -L$(LIBDIR) -l$(LIBUTIL)
LIBFUT  = $(LIBDIR)/lib$(LIBUTIL).a

LIBLSU  = -L$(LIBDIR) -l$(LIBSHYUTIL)
LIBFSU  = $(LIBDIR)/lib$(LIBSHYUTIL).a

LIBLMP  = -L$(LIBDIR) -l$(LIBSHYMPI)
LIBFMP  = $(LIBDIR)/lib$(LIBSHYMPI).a

LIBLUM  = -L$(LIBDIR) -l$(LIBSHYUTILMPI)
LIBFUM  = $(LIBDIR)/lib$(LIBSHYUTILMPI).a

LIBL = $(LIBLUM) $(LIBLMP) $(LIBLSU) $(LIBLUT)
LIBF = $(LIBFUM) $(LIBFMP) $(LIBFSU) $(LIBFUT)

#---------------------------------------------

EXES = shyelab
EXELINKS =      $(FEMBIN) $(EXES)

#----------------------------------------------------------------------------

ifneq ($(NETCDF),true)
  NETCDF_OBJ= \
	$(FEMDIR)/tools/shyelab/elab_nc_dummy.o
else
  LIB_NETCDF = -L$(NETCDFDIR)/lib -lnetcdf -L$(NETCDFFDIR)/lib -lnetcdff
  NETCDF_OBJ=\
	$(FEMDIR)/tools/shyelab/elab_nc.o \
	$(FEMDIR)/tools/2nc/netcdf.o \
	$(FEMDIR)/tools/2nc/netcdf_util.o 
endif

ifneq ($(PARALLEL_MPI),NONE)
  LIB_MPI = -L$(MPI_DIR)/lib -lmpi
  SHYMPI_OBJ=\
	$(FEMDIR)/src/mpi_node.o \
	$(FEMDIR)/src/mpi_node_internal.o \
	$(FEMDIR)/src/mpi_general.o 
else
  SHYMPI_OBJ=\
	$(FEMDIR)/src/mpi_dummy.o \
	$(FEMDIR)/src/mpi_general.o 
endif

#----------------------------------------------------------------------------

LFLAGS:= $(LIB_NETCDF) $(LIB_MPI)

#----------------------------------------------------------------------------

OBJ := \
        $(SHYMPI_OBJ) \
	$(NETCDF_OBJ) \
	../../tools/shyelab/elabtime.o \
	../../tools/shyelab/elab_off.o \
	../../tools/shyelab/elabutil.o \
	../../tools/shyelab/elab_gis.o \
	../../tools/shyelab/elab_expand.o \
	../../tools/shyelab/elab_regular.o \
	../../tools/shyelab/shyelab_output.o \
	../../tools/shyelab/shyelab_util.o \
	../../tools/shyelab/shyelab_area.o \
	../../tools/shyelab/shyelab_average.o \
	../../tools/shyelab/shyelab_nodes.o \
	../../tools/shyelab/shyelab_proj.o \
	../../tools/shyelab/shyelab_extract.o \
	../../tools/shyelab/shyelab1.o \
	../../tools/shyelab/lgrelab1.o \
	../../tools/shyelab/extelab1.o \
	../../tools/shyelab/flxelab1.o \
	../../tools/shyelab/elab_resample.o \
	../../tools/shyelab/elab_check.o \
	../../tools/shyelab/tselab1.o \
	../../tools/shyelab/femelab1.o \
	../../tools/2nc/output_util.o \
	../../src/ioshy.o \
	../../src/mod_offline.o \
	../../src/shy_util.o \
	../../src/zadaptutil.o \
	../../src/depth_util.o \
	../../src/heat_util.o \
	../../src/tsfile.o \
	../../src/ioflux.o \
	../../src/ioext.o \
	../../src/mpi_buffer.o \
	../../src/def_para.o \
	../../src/lagrange_read.o \
	../../src/dates.o \
	../../fem3d/unitutils.o \
	../../fem3d/shyutil.o \
	../../fem3d/subshyutil2.o \

OBSOLETE := \
	../../src/to_linear.o \
	../../src/coordinates.o \
	../../src/projection.o \
	../../src/femfile.o \
	../../src/iso8601.o \
	../../src/defnames.o \
	../../src/geometrical.o \
	../../src/fem_util.o \
	../../src/sort_util.o \
	../../src/sigma.o \
	../../src/util.o \
	../../src/evgeom.o \
	../../src/geom_dynamic.o \
	../../src/basin.o \
	../../src/levels.o \
	../../src/depth.o \
	../../src/clo.o \
	../../src/strings.o \
	../../src/line_util.o \
	../../src/round.o \
	../../src/regular.o \
	../../src/grd.o \
	../../src/chk_NaN.o \
	../../src/vec_util.o \
	../../src/dts.o \
	../../src/utility.o \
	../../src/para.o \
	../../src/nls.o \
	../../src/fil.o \
	../../src/convert.o \
	../../src/version.o \

#----------------------------------------------------------------------------

default: fem
fem: $(EXES) links

links:
	@$(FEMBIN)/symlink $(EXELINKS)

$(EXES): $(OBJ) $(LIBF) $(EXES).o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $^ $(LIBL) -o $@
	cp $(EXES) ../../bin/	

#----------------------------------------------------------------------------

clean: cleanobj cleanlinks cleanexe

cleanobj:
	-rm -f *.o
	-rm -f *.mod
	-rm -f *_genmod.f90

cleanlinks:
	@$(FEMBIN)/symlink -remove $(EXELINKS)

cleanexe:
	-rm -f $(EXES)

cleanall: clean

#----------------------------------------------------------------------------

%.o:%.f
	${F77} $(FFLAGS) $(INCLUDE) -c $< -o $@

%.o:%.f90
	${F77} $(FFLAGS) $(INCLUDE) -c $< -o $@

%.o:%.c
	${CC} $(CFLAGS) $(INCLUDE) -c $< -o $@

#----------------------------------------------------------------------------

