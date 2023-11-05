
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2023  The SHYFEM team
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

DIR	= shybas

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

#----------------------------------------------------------------------------

EXES = shybas
EXELINKS =      $(FEMBIN) $(EXES)

#----------------------------------------------------------------------------

OBJ := \
	msh.o \
	basutil.o \
	basutil_smooth.o \
	basutil_resol.o \
	basutil_box.o \
	basutil_custom.o \
	basutil_bathy.o \
	basutil_interp.o \
	basutil_hsigma.o \
	basutil_part.o \

#----------------------------------------------------------------------------

default: fem
fem: $(EXES) links

links:
	@$(FEMBIN)/symlink $(EXELINKS)

$(EXES): $(OBJ) $(LIBF) $(EXES).o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $^ $(LIBL) -o $@

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

