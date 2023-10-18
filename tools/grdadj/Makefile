
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2023  The SHYFEM team
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

LIBPOST     = post
LIBUTIL     = util
LIBSHYUTIL  = shyutil

LIBLPS  = -L$(DIRLIB) -l$(LIBPOST)
LIBFPS  = $(DIRLIB)/lib$(LIBPOST).a

LIBLUT  = -L$(DIRLIB) -l$(LIBUTIL)
LIBFUT  = $(DIRLIB)/lib$(LIBUTIL).a

LIBLSU  = -L$(DIRLIB) -l$(LIBSHYUTIL)
LIBFSU  = $(DIRLIB)/lib$(LIBSHYUTIL).a

LIBL = $(LIBLPS) $(LIBLUT) $(LIBLSU)
LIBF = $(LIBFPS) $(LIBFUT) $(LIBFSU)

#--------------------------------------------- macros

OBJS =

ADJOBJS	=	mod_adj_grade.o \
		adj4el.o adj5el.o adj7el.o adj575.o \
		adjplo.o adjgrd.o adjneu.o adjvar.o

ALLOBJS =	$(OBJS) $(ADJOBJS) shyadj.o

EXES	=	shyadj
EXELINKS = 	$(FEMBIN) $(EXES)

SPECIAL =       Makefile TODO README
SAVEFILES =     *.f90 $(SPECIAL)

MAKEDEPEND = $(FEMBIN)/mkdp.pl

#--------------------------------------------- EXES

default: fem
fem:	shyadj links
all:	$(EXES)

links:
	@$(FEMBIN)/symlink $(EXELINKS)

shyadj: $(ALLOBJS) $(LIBF)
	$(LINKER) $(LFLAGS) $(ALLOBJS) $(LIBL) -o $@

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

