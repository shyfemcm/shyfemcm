
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2023  The SHYFEM team
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

#--------------------------------------------- 

DIR	= bfm

#--------------------------------------------- just for default
FEMDIR = ../../../..
FEMBIN = $(FEMDIR)/bin
LIBDIR = $(FEMDIR)/lib
F77    = g77
FFLAGS =  -O -Wall -g
LINKER = $(F77)
#---------------------------------------------

include $(FEMDIR)/Rules.make

#--------------------------------------------- libraries

LIBPOST     = post
LIBUTIL     = util
LIBSHYUTIL  = shyutil

LIBL	 = -L$(LIBDIR)

LIBL	+= -l$(LIBPOST)
LIBF	+= $(LIBDIR)/lib$(LIBPOST).a

LIBL	+= -l$(LIBUTIL)
LIBF	+= $(LIBDIR)/lib$(LIBUTIL).a

LIBL	+= -l$(LIBSHYUTIL)
LIBF	+= $(LIBDIR)/lib$(LIBSHYUTIL).a

#--------------------------------------------- macros

EXES	=	

OBJ	=	bfm_statevars.o bfm.o bfm_internal.o 
MOD_BFM =  $(FEMDIR)/lib/mod/mod_bfm.mod $(FEMDIR)/lib/mod/mod_bfm_internal.mod
LIBG_BFM     = -L$(BFMDIR)/lib/ -l$(LIBBFM)
LIBF_BFM     = $(BFMDIR)/lib/lib$(LIBBFM).a
INCL_BFM     = -I$(BFMDIR)/include
FFLAGS          += -fpp $(INCL_BFM) $(LIBG_BFM) $(LIBF_BFM)

ALL	=	$(OBJ)

SPECIAL =       Makefile TODO README
SAVEFILES =     *.f90 $(SPECIAL)

MAKEDEPEND = $(FEMBIN)/mkdp.pl

#--------------------------------------------- EXES

default: fem
fem:	$(OBJ) 
all:	$(EXES)

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
	@-$(MAKEDEPEND) -moddir=$(MODDIR) *.f90

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

../../../lib/mod/mod_adj_grade.mod: mod_adj_grade.o
../../../lib/mod/mod_adj_static.mod: mod_adj_grade.o
adj4el.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/mod_adj_grade.mod 
adj575.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/mod_adj_grade.mod 
adj5el.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/mod_adj_grade.mod 
adj7el.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/mod_adj_grade.mod 
adjgrd.o: ../../../lib/mod/mod_adj_static.mod
adjneu.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/mod_adj_grade.mod \
		../../../lib/mod/mod_depth.mod 
adjplo.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/mod_adj_grade.mod 
adjvar.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/mod_adj_grade.mod \
		../../../lib/mod/mod_adj_static.mod 
shyadj.o: ../../../lib/mod/basin.mod ../../../lib/mod/clo.mod \
		../../../lib/mod/mod_adj_grade.mod \
		../../../lib/mod/mod_depth.mod 

