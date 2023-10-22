
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2023  The SHYFEM team
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

DIR     = shyplot

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
LIBPOST         = post

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

LIBLPS  = -L$(LIBDIR) -l$(LIBPOST)
LIBFPS  = $(LIBDIR)/lib$(LIBPOST).a

#LIBL = $(LIBLNC) $(LIBLUM) $(LIBLMP) $(LIBLSU) $(LIBLUT) $(LIB_NETCDF)
LIBL = $(LIBLNC) $(LIBLUM) $(LIBLMP) $(LIBLSU) $(LIBLUT) $(LIBLPS)
LIBF = $(LIBFNC) $(LIBFUM) $(LIBFMP) $(LIBFSU) $(LIBFUT)

#---------------------------------------------

EXES = shyplot
EXELINKS =      $(FEMBIN) $(EXES)

#----------------------------------------------------------------------------

LFLAGS:= $(LIB_NETCDF) $(LIB_MPI)

#----------------------------------------------------------------------------

OBJ := \
	col.o \
	supcust.o \
	mod_hydro_plot.o \
	$(FEMDIR)/tools/shyelab/elabtime.o \
	$(FEMDIR)/tools/shyelab/elabutil.o \
	$(FEMDIR)/src/hydro_admin.o \
	$(FEMDIR)/src/hydro_vel.o \
	$(FEMDIR)/src/hydro_print.o \
	plotutil.o \
	supbas.o \
	supano.o \
	supdep.o \
	supint.o \
	supiso.o \
	suplag.o \
	suplin.o \
	suplin1.o \
	supmpi.o \
	suprgf.o \
	supsim.o \
	suptim.o \
	supver.o \
	supout.o \
	util/suputi.f90 \
	util/supcol.f90 \
	genutil.o \
	def_para_post.o \
	subnsa.o \

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

