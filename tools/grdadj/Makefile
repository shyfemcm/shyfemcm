
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

#--------------------------------------------- 

DIR	= gridadj

#--------------------------------------------- just for default
FEMDIR = ../..
FEMBIN = $(FEMDIR)/bin
DIRLIB = $(FEMDIR)/lib
F77    = g77
FFLAGS =  -O -Wall -g
LINKER = $(F77)
#---------------------------------------------

include $(FEMDIR)/Rules.make

#--------------------------------------------- libraries

LIBFEM  = fem
LIBPOST = post

LIBFS   = -L$(DIRLIB) -l$(LIBFEM)
LIBIFS  = $(DIRLIB)/lib$(LIBFEM).a
LIBFS   = 
LIBIFS  = 

LIBPS   = -L$(DIRLIB) -l$(LIBPOST)
LIBIPS  = $(DIRLIB)/lib$(LIBPOST).a

#--------------------------------------------- macros

OBJS =		basin.o depth.o clo.o convert.o vec_util.o \
		geometrical.o fil.o geom_dynamic.o \
		grd.o grdutil.o version.o \
		sort_util.o

ADJOBJS	=	mod_adj_grade.o \
		adj4el.o adj5el.o adj7el.o adj575.o \
		adjplo.o adjgrd.o adjneu.o adjvar.o

ALLOBJS =	$(OBJS) $(ADJOBJS) shyadj.o

EXES	=	shyadj
EXELINKS = 	$(FEMBIN) $(EXES)

SPECIAL =       Makefile TODO COMMIT VERSION
SAVEFILES =     *.[fh] $(SPECIAL)

MAKEDEPEND = $(FEMBIN)/mkdp.pl

#--------------------------------------------- EXES

default: fem
fem:	shyadj links
all:	$(EXES)

links:
	@$(FEMBIN)/symlink $(EXELINKS)

shyadj: $(ALLOBJS) $(LIBIFS) $(LIBIPS)
	$(LINKER) $(LFLAGS) $(ALLOBJS) $(LIBFS) $(LIBPS) -o $@

#--------------------------------------------- clean

cleanall: clean cleanzip cleangrd cleandiff

clean: cleanobj cleanexe cleanvar cleantex cleanlinks

cleanobj:
	-rm -f *.o
	-rm -f *.stb
	-rm -f *.mod
	-rm -f *_genmod.f90

cleanlinks:
	@$(FEMBIN)/symlink -remove $(EXELINKS)

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
	-rm -f .memory
	-rm -fr cvsdiff_tmp

cleandiff:
	-rm -f DIFFS

cleanzip:
	-rm -f $(DIR).zip
	-rm -f save.zip

cleangrd:
	-rm -f *.grd

#--------------------------------------------- various

param:

diff:
	@$(FEMBIN)/cvsdiff -norcs -ifdiff * > DIFFS
	@echo "differences written to file DIFFS"

depend:
	@$(MAKEDEPEND) -moddir=$(MODDIR) *.f *.f90

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

