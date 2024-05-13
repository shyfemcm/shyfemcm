
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2023  The SHYFEM team
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

DIR	= shypre

#--------------------------------------------- just for default
FEMDIR = ../../..
FEMBIN = $(FEMDIR)/bin
LIBDIR = $(FEMDIR)/lib
F77    = g77
FFLAGS =  -O -Wall -g
LINKER = $(F77)
LFLAGS:= 
INCLUDE = -I$(FEMDIR)/include/
#---------------------------------------------

include $(FEMDIR)/Rules.make

#---------------------------------------------

LIBUTIL       = util
LIBSHYUTIL    = shyutil
LIBSHYMPI     = shympi
LIBSHYUTILMPI = shyutilmpi

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

EXES = shypre shypre_orig shypre_noopti

OBJ_ORIG := rosen.o cuthill_mckee.o
OBJ_OPTI := \
		$(FEMDIR)/src/shyfem/shypre_subs.o \
		$(FEMDIR)/src/shyfem/shypre_opti.o \

OBJ_OLD := 	shypre_opti_subs.o rosen.o cuthill_mckee.o

OBJ_NOOPTI := \
		$(FEMDIR)/src/shyfem/shypre_subs.o \
		$(FEMDIR)/src/shyfem/shypre_noopti.o

MAKEDEPEND = $(FEMBIN)/mkdp.pl

#----------------------------------------------------------------------------

default: fem
fem: $(EXES) links

links:
	@$(FEMBIN)/symlink $(FEMBIN) $(EXES)

#$(EXES): $(OBJ) $(LIBF) $(EXES).o
#	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $^ $(LIBL) -o $@

shypre_orig: $(OBJ_ORIG) $(LIBF) shypre_orig.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $^ $(LIBL) -o $@

shypre: $(OBJ_OPTI) $(LIBF) shypre.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $^ $(LIBL) -o $@

shypre_noopti: $(OBJ_NOOPTI) $(LIBF) shypre_noopti.o
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $^ $(LIBL) -o $@

#----------------------------------------------------------------------------

clean: cleanvar cleanobj cleanlinks cleanexe

cleanvar:
	-rm -f ggg hhh
	-rm -f *.tmp *.bak

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

depend:
	@-$(MAKEDEPEND) -moddir=$(MODDIR) *.f90

#----------------------------------------------------------------------------

%.o:%.f
	${F77} $(FFLAGS) $(INCLUDE) -c $< -o $@

%.o:%.f90
	${F77} $(FFLAGS) $(INCLUDE) -c $< -o $@

%.o:%.c
	${CC} $(CFLAGS) $(INCLUDE) -c $< -o $@

#----------------------------------------------------------------------------

# DO NOT DELETE THIS LINE -- make depend depends on it.

../../../lib/mod/mod_shypre.mod: shypre.o
../../../lib/mod/mod_shypre.mod: shypre_noopti.o
../../../lib/mod/mod_shypre.mod: shypre_orig.o
shypre.o: ../../../lib/mod/clo.mod
shypre_noopti.o: ../../../lib/mod/clo.mod
shypre_opti.o: ../../../lib/mod/basin.mod ../../../lib/mod/grd.mod
shypre_opti_subs.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/grd.mod 
shypre_orig.o: ../../../lib/mod/basin.mod ../../../lib/mod/clo.mod \
		../../../lib/mod/evgeom.mod \
		../../../lib/mod/grd.mod \
		../../../lib/mod/mod_sort.mod 

