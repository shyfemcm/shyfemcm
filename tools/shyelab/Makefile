
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

LIB_NETCDF =
ifeq ($(NETCDF),true)
  LIB_NETCDF = -L$(NETCDFDIR)/lib -lnetcdf -L$(NETCDFFDIR)/lib -lnetcdff
endif

LIBUTIL         = util
LIBSHYUTIL      = shyutil
LIBSHYMPI       = shympi
LIBSHYUTILMPI   = shyutilmpi
LIBNETCDF       = shynetcdf

LIBLUT  = -L$(LIBDIR) -l$(LIBUTIL)
LIBFUT  = $(LIBDIR)/lib$(LIBUTIL).a

LIBLSU  = -L$(LIBDIR) -l$(LIBSHYUTIL)
LIBFSU  = $(LIBDIR)/lib$(LIBSHYUTIL).a

LIBLMP  = -L$(LIBDIR) -l$(LIBSHYMPI)
LIBFMP  = $(LIBDIR)/lib$(LIBSHYMPI).a

LIBLUM  = -L$(LIBDIR) -l$(LIBSHYUTILMPI)
LIBFUM  = $(LIBDIR)/lib$(LIBSHYUTILMPI).a

LIBLNC  = -L$(LIBDIR) -l$(LIBNETCDF)
LIBFNC  = $(LIBDIR)/lib$(LIBNETCDF).a

LIBL = $(LIBLNC) $(LIBLUM) $(LIBLMP) $(LIBLSU) $(LIBLUT) $(LIB_NETCDF)
LIBF = $(LIBFNC) $(LIBFUM) $(LIBFMP) $(LIBFSU) $(LIBFUT)

#---------------------------------------------

EXES = shyelab
EXELINKS =      $(FEMBIN) $(EXES)

#----------------------------------------------------------------------------

LFLAGS:= $(LIB_NETCDF) $(LIB_MPI)

#----------------------------------------------------------------------------

OBJ := \
	elabutil.o \
	elabtime.o \
	elab_off.o \
	elab_gis.o \
	elab_regular.o \
	elab_resample.o \
	elab_check.o \
	shyelab_output.o \
	shyelab_util.o \
	shyelab_area.o \
	shyelab_average.o \
	shyelab_nodes.o \
	shyelab_proj.o \
	shyelab_extract.o \
	shyelab1.o \
	extelab1.o \
	flxelab1.o \
	tselab1.o \
	femelab1.o \
	lgrelab1.o \
	shy_util2.o \

OBSOLETE := \
	$(NETCDF_OBJ) \
        $(SHYMPI_OBJ) \
	../../tools/2nc/output_util.o \
	../../fem3d/shyutil.o \
	../../fem3d/subshyutil2.o \
	../../fem3d/unitutils.o \
	../../src/def_para.o \
	../../src/dates.o \
	../../src/mpi_buffer.o \
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

