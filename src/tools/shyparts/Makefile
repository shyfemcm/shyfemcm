
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2023  The SHYFEM team
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

#--------------------------------------------- 

DIR	= shyparts

#--------------------------------------------- just for default
FEMDIR = ../../..
FEMBIN = $(FEMDIR)/bin
LIBDIR = $(FEMDIR)/lib
F77    = g77
FFLAGS =  -O -Wall -g
LINKER = $(F77)
#---------------------------------------------

include $(FEMDIR)/Rules.make

#--------------------------------------------- libraries

ifeq ($(PARTS),$(filter $(PARTS),METIS PARMETIS))
  LIBG_METIS = -L$(METISDIR)/lib -lmetis
endif

LIBUTIL     = util
LIBSHYUTIL  = shyutil
LIBSHYMPI   = shympi

LIBL	 = -L$(LIBDIR)

LIBL	+= -l$(LIBSHYMPI)
LIBF	+= $(LIBDIR)/lib$(LIBSHYMPI).a

LIBL	+= -l$(LIBSHYUTIL)
LIBF	+= $(LIBDIR)/lib$(LIBSHYUTIL).a

LIBL	+= -l$(LIBUTIL)
LIBF	+= $(LIBDIR)/lib$(LIBUTIL).a

LIBL	+= $(LIBG_METIS)

#--------------------------------------------- macros

EXES	=	shyparts

XXOBJ	=	mod_adj_grade.o \
		adj4el.o adj5el.o adj7el.o adj575.o \
		adjplo.o adjgrd.o adjneu.o adjvar.o

ALL	=	$(OBJ)

SPECIAL =       Makefile TODO README
SAVEFILES =     *.f90 $(SPECIAL)

MAKEDEPEND = $(FEMBIN)/mkdp.pl

#--------------------------------------------- EXES

default: fem
fem:	shyparts links
all:	$(EXES)

links:
	@$(FEMBIN)/symlink $(FEMBIN) $(EXES)

shyparts: $(ALL) $(LIBF) shyparts.o
	$(LINKER) $(LFLAGS) $(ALL) $@.o $(LIBL) -o $@

#--------------------------------------------- clean

cleanall: clean cleanzip cleangrd

clean: cleanobj cleanexe cleanvar cleantex cleanlinks

cleanobj:
	-rm -f *.o
	-rm -f *.stb
	-rm -f *.mod
	-rm -f *_genmod.f90

cleanlinks:
	@$(FEMBIN)/symlink -remove $(FEMBIN) $(EXES)

cleanexe:
	-rm -f $(EXES)
	-rm -f *.exe

cleantex:
	-rm -f *.dvi *.aux

cleanvar:
	-rm -f a.out core
	-rm -f *.bak *.tmp
	-rm -f ggg hhh
	-rm -f mon.out nohup.out
	-rm -f fort.*
	-rm -f errout.dat
	-rm -f gmon.out
	-rm -f out.ps plot.ps

cleanzip:
	-rm -f $(DIR).zip
	-rm -f save.zip

cleangrd:
	-rm -f *.grd

#--------------------------------------------- various

depend:
	@$(MAKEDEPEND) -moddir=$(MODDIR) *.f90

list:
	pwd; ls

save:   zip
zip:	cleanall
	-zip $(DIR) $(SAVEFILES)

strip:
	-strip $(EXES)

#--------------------------------------------- rules

.SUFFIXES: .F .f90

.c.o:
	$(CC) -c $(CFLAGS) $*.c

.f.o:
	$(F77) -c $(FFLAGS) $<

.f90.o:
	$(F77) -c $(FFLAGS) $<

.F.o:
	$(F77) -c -cpp $(FFLAGS) $<

.for.obj:
	$(FLC) $< , $(<,B).OBJ $(FLFLAGS)

#--------------------------------------------- makedepend

# DO NOT DELETE THIS LINE -- make depend depends on it.

../../../lib/mod/mod_shyparts.mod: shyparts.o
shyparts.o: ../../../lib/mod/basin.mod ../../../lib/mod/clo.mod \
		../../../lib/mod/evgeom.mod \
		../../../lib/mod/grd.mod \
		../../../lib/mod/mod_geom.mod \
		../../../lib/mod/shympi.mod 

