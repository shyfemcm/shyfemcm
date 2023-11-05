
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

#------------------------------------------------------------

DIR	= femplot

#------------------------------------------------------------ just for default
FEMDIR = ../../../..
FEMBIN = $(FEMDIR)/bin
DIRLIB = $(FEMDIR)/lib
F77    = g77
FFLAGS =  -O
LINKER = $(F77)
#------------------------------------------------------------

include $(FEMDIR)/Rules.make

#------------------------------------------------------------

OBJ := \
	../col.o \
	../color_table.o \
	../plotgeneric.o \

LIBUTIL         = util
LIBSHYUTIL      = shyutil
LIBSHYMPI       = shympi
LIBSHYUTILMPI   = shyutilmpi
LIBPOST         = post

LIBLUT  = -L$(LIBDIR) -l$(LIBUTIL)
LIBFUT  = $(LIBDIR)/lib$(LIBUTIL).a

LIBLSU  = -L$(LIBDIR) -l$(LIBSHYUTIL)
LIBFSU  = $(LIBDIR)/lib$(LIBSHYUTIL).a

LIBLMP  = -L$(LIBDIR) -l$(LIBSHYMPI)
LIBFMP  = $(LIBDIR)/lib$(LIBSHYMPI).a

LIBLUM  = -L$(LIBDIR) -l$(LIBSHYUTILMPI)
LIBFUM  = $(LIBDIR)/lib$(LIBSHYUTILMPI).a

LIBLPS  = -L$(LIBDIR) -l$(LIBPOST)
LIBFPS  = $(LIBDIR)/lib$(LIBPOST).a

LIBL = $(LIBLNC) $(LIBLUM) $(LIBLMP) $(LIBLSU) $(LIBLUT) $(LIBLPS)
LIBF = $(LIBFNC) $(LIBFUM) $(LIBFMP) $(LIBFSU) $(LIBFUT)

#------------------------------------------------------------

SPECIAL = Makefile VERSION COMMIT

EXES    = 
EXESAUX = plotcolb windrose testdash plotbug lagxitest
EXESALL = $(EXES) $(EXESAUX)
EXESLINKS = $(EXES)

MAKEDEPEND = $(FEMBIN)/mkdp.pl

#------------------------------------------------------------

default: all

fem: all

all: $(EXESALL) links

help:
	@echo "help         this screen"
	@echo "all          compile utility programs"

plotcolb: $(OBJ) plotcolb.o $(LIBF)
	$(LINKER) $(LFLAGS) $@.o $(OBJ) $(LIBL) -o $@

windrose: windrose.o $(LIBF)
	$(LINKER) $(LFLAGS) $@.o $(LIBL) -o $@

testdash: testdash.o  $(LIBF)
	$(LINKER) $(LFLAGS) $@.o $(LIBL) -o $@

plotbug: plotbug.o $(LIBF)
	$(LINKER) $(LFLAGS) $@.o $(LIBL) -o $@

lagxitest: lagxitest.o $(LIBF)
	$(LINKER) $(LFLAGS) $@.o $(LIBL) -o $@

links:

#------------------------------------------------------------

cleanall: clean cleanzip cleandiff

clean: cleanobj cleanexe cleanvar

cleanobj:
	-rm -f *.o
	-rm -f *.stb
	-rm -f *.mod
	-rm -f *_genmod.f90

cleanexe:
	-rm -f $(EXESALL)
	-rm -f *.exp
	-rm -f *.exe

cleanvar:
	-rm -f a.out core
	-rm -f *.bak *.tmp
	-rm -f ggg hhh
	-rm -f mon.out nohup.out gmon.out
	-rm -f fort.*
	-rm -f *.0
	-rm -f *.ps
	-rm -f .memory
	-rm -f *~
	-rm -rf cvsdiff_tmp

cleandiff:
	-rm -f DIFFS

cleanzip:
	-rm -f $(DIR).zip
	-rm -f save.zip

cleanlinks:
	@$(FEMBIN)/symlink -remove $(FEMBIN) $(EXESALL)

#------------------------------------------------------------

param:
	
diff:
	@$(FEMBIN)/cvsdiff -norcs -ifdiff * > DIFFS
	@echo "differences written to file DIFFS"

depend:
	@$(MAKEDEPEND) -moddir=$(MODDIR) *.f
	
list:
	pwd; ls

rcs:
	ci -l *.f *.h *.c $(SPECIAL)

save:   cleanall
	@echo $(DIR)
	-zip $(DIR) *.[fhc]
	-zip $(DIR) $(SPECIAL)

zip:	save

strip:	$(EXES)
	-strip $(EXES)

#------------------------------------------------------------
# Rules -----------------------------------------------------
#------------------------------------------------------------

.SUFFIXES: .f90

.c.o:
	$(CC)  -c $(CFLAGS) $<

.f.o:
	$(F77) -c $(FFLAGS) $<

.f90.o:
	$(F77) -c $(FFLAGS) $<

#------------------------------------------------------------

# DO NOT DELETE THIS LINE -- make depend depends on it.

../../femlib/mod/mod_testdash.mod: testdash.o

