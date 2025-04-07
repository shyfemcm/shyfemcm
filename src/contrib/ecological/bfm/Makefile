
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2023  The SHYFEM team
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

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

OBJ = bfm_dummy.o
ifeq ($(ECOLOGICAL),BFM)
  OBJ	     =	bfm_statevars.o bfm.o bfm_internal.o 
  LIBG_BFM     = -L$(BFMDIR)/lib/ -l$(LIBBFM)
  INCL_BFM     = -I$(BFMDIR)/include
  FFLAGS      += $(INCL_BFM)
  LFLAGS      += $(LIBG_BFM)
endif

BFMOBJ := $(foreach item,$(OBJ),$(PWD)/$(item))

#$(info ------------------------------ )
#$(info LIBG_BFM = $(LIBG_BFM) )
#$(info INCL_BFM = $(INCL_BFM) )
#$(info ------------------------------ )

ALL	=	$(OBJ)

SPECIAL =       Makefile TODO README
SAVEFILES =     *.f90 $(SPECIAL)

MAKEDEPEND = $(FEMBIN)/mkdp.pl

#--------------------------------------------- EXES

default: fem
fem:	$(BFMPRIVATE) $(OBJ) nothing

nothing:
	@true

objs:
	@echo $(BFMOBJ)

info:
	@echo $(PWD)
	@echo $(BFM_CONTRIB_DIR)
	@echo $(OBJ)
	@echo $(BFMOBJ)

#--------------------------------------------- private stuff

BFMPRIVATE = bfm_statevars.f90 bfm.f90 bfm_internal.f90

$(BFMPRIVATE):
	cp -f arc/$@ .

bfmcopy:
	cp -f arc/*.f90 .

bfmsave:
	cp -f $(BFMPRIVATE) arc

bfmdiff:
	for file in $(BFMPRIVATE) ; do \
	  echo $$file ; \
	  diff arc $$file ; \
	done

bfmzip:
	zip bfm.zip $(BFMPRIVATE)

bfmclean: 
	-rm -f $(BFMPRIVATE)
	-rm -f bfm.zip

#--------------------------------------------- clean

cleanall: clean cleanzip cleangrd

clean: cleanobj cleanexe cleanvar cleanlinks cleanbfm

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

cleanbfm: bfmdelete

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

../../../../lib/mod/mod_bfm.mod: bfm_dummy.o
bfm_dummy.o: ../../../../lib/mod/levels.mod \
		../../../../lib/mod/mod_diff_visc_fric.mod \
		../../../../lib/mod/mod_layer_thickness.mod \
		../../../../lib/mod/mod_meteo.mod \
		../../../../lib/mod/mod_ts.mod 

