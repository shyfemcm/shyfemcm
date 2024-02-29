
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2023  The SHYFEM team
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

DIR     = shyelab

#--------------------------------------------- just for default
FEMDIR = ../../..
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

LIBL     = -L$(LIBDIR)

LIBL    += -l$(LIBNETCDF)
LIBF    += $(LIBDIR)/lib$(LIBNETCDF).a

LIBL    += -l$(LIBSHYUTILMPI)
LIBF    += $(LIBDIR)/lib$(LIBSHYUTILMPI).a

LIBL    += -l$(LIBSHYMPI)
LIBF    += $(LIBDIR)/lib$(LIBSHYMPI).a

LIBL    += -l$(LIBSHYUTIL)
LIBF    += $(LIBDIR)/lib$(LIBSHYUTIL).a

LIBL    += -l$(LIBUTIL)
LIBF    += $(LIBDIR)/lib$(LIBUTIL).a

LIBL    += $(LIB_NETCDF)

#---------------------------------------------

EXES = shyelab extelab flxelab lgrelab tselab femelab shyproj

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

#----------------------------------------------------------------------------

default: fem
fem: $(EXES) links

links:
	@$(FEMBIN)/symlink $(FEMBIN) $(EXES)

shyelab: $(OBJ) $(LIBF) shyelab.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $(OBJ) $@.o $(LIBL) -o $@

extelab: $(OBJ) $(LIBF) extelab.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $(OBJ) $@.o $(LIBL) -o $@

flxelab: $(OBJ) $(LIBF) flxelab.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $(OBJ) $@.o $(LIBL) -o $@

lgrelab: $(OBJ) $(LIBF) lgrelab.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $(OBJ) $@.o $(LIBL) -o $@

femelab: $(OBJ) $(LIBF) femelab.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $(OBJ) $@.o $(LIBL) -o $@

tselab: $(OBJ) $(LIBF) tselab.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $(OBJ) $@.o $(LIBL) -o $@

shyproj: $(LIBF) shyproj.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $@.o $(LIBL) -o $@

#----------------------------------------------------------------------------

clean: cleanobj cleanlinks cleanexe

cleanobj:
	-rm -f *.o
	-rm -f *.mod
	-rm -f *_genmod.f90

cleanlinks:
	@$(FEMBIN)/symlink -remove $(FEMBIN) $(EXES)

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

