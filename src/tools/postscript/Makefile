
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

#-----------------------------------------------------------

DIR    = postscript

#----------------------------------------------------------- just for default
FEMDIR = ../../..
FEMBIN = $(FEMDIR)/bin
LIBDIR = $(FEMDIR)/lib
CC     = gcc
CFLAGS =  -O -Wall -pedantic
F77    = g77
FFLAGS =  -O -Wall
#-----------------------------------------------------------

include $(FEMDIR)/Rules.make

#-----------------------------------------------------------

LIBPOST     = post
LIBCALP     = calp

LIBLPS  = -L$(LIBDIR) -l$(LIBPOST)
LIBFPS  = $(LIBDIR)/lib$(LIBPOST).a

LIBLCP  = -L$(LIBDIR) -l$(LIBCALP)
LIBFCP  = $(LIBDIR)/lib$(LIBCALP).a

LIBL = $(LIBLPS) $(LIBLCP)
LIBF = $(LIBFPS) $(LIBFCP)

#--------------------------------

EXES = democalp demopost psdemo colorbar pf pc layout

SPECIAL = Makefile CALCOMP pspost.tex ps/cnr.ps

MAKEDEPEND = $(FEMBIN)/mkdp.pl

#---------------------------

default: fem

fem: libs demopost nothing

all: libs $(EXES) nothing

nothing:
	@true

democalp: democalp.o pgraph.o general.o pcalp.o
	$(F77) $@.o pgraph.o general.o pcalp.o -o $@

demopost: demopost.o pgraph.o general.o
	$(CC) $@.o pgraph.o general.o -o $@

psdemo: psdemo.o $(LIBF)
	$(F77) $@.o $(LIBL) -o $@

colorbar: colorbar.o $(LIBF)
	$(F77) $@.o $(LIBL) -o $@

pf: pf.o $(LIBF)
	$(F77) $@.o $(LIBL) -o $@

pc: pc.o $(LIBF)
	$(CC) $@.o $(LIBL) -o $@

layout: layout.o
	$(CC) $@.o -o $@

#---------------------------

libs: $(LIBFPS) $(LIBFCP)

$(LIBFPS): general.o gustd.o psgraph.o psgraphf.o colorutil.o
	ar rvs $(LIBFPS) *.o

$(LIBFCP): general.o pgraph.o pcalp.o pcalpf.o
	ar rvs $(LIBFCP) *.o

#----------------------------- documentation

docs: tex
latex: tex
tex:
	latex pspost
	latex pspost
	latex pspost
	pdflatex pspost

ps:
	dvips pspost

#----------------------------- dependencies

diff:
	@$(FEMBIN)/cvsdiff -norcs -ifdiff -quiet * > DIFFS
	@echo "differences written to file DIFFS"

depend:
	@$(MAKEDEPEND) *.[cf]

#---------------------------

cleanall : clean cleanlib cleanzip cleandiff

clean: cleanobj cleanexe cleanvar cleantex

cleanlib:
	-rm -f $(LIBFPS)
	-rm -f $(LIBFCP)

cleanexe:
	-rm -f $(EXES)
	-rm -f *.exe

cleanobj:
	-rm -f *.o
	-rm -f *.stb
	-rm -f *.mod
	-rm -f *_genmod.f90

cleanvar:
	-rm -f a.out core
	-rm -f *.bak
	-rm -f ggg hhh
	-rm -f mon.out nohup.out
	-rm -f fort.* *.tmp
	-rm -f plot.ps plot.*.ps demo.ps
	-rm -f plot.pdf plot.*.pdf
	-rm -f out.ps out.*.ps
	-rm -f *.stackdump
	-rm -fr cvsdiff_tmp

cleantex:
	-rm -f pspost.aux pspost.log pspost.dvi pspost.ps pspost.toc pspost.pdf

cleandiff:
	-rm -f DIFFS

cleanzip:
	-rm -f *.zip

#----------------------------

cvschange:
	cvschange -d

list:
	pwd; ls

zip: save
save:   cleanall
	@echo $(DIR)
	-zip $(DIR) *.[fch] $(SPECIAL)

#----------------------------

.c.o:
	$(CC)  -c -g $(CFLAGS) $<
		 
.f.o:
	$(F77) -c -g $(FFLAGS) $<

#----------------------------

# DO NOT DELETE THIS LINE -- make depend depends on it.

colorutil.o: general.h
demopost.o: general.h pgraph.h
general.o: general.h
gustd.o: general.h gustd.h
pc.o: pspost.h
pcalp.o: general.h pgraph.h
pgraph.o: general.h pgraph.h
psgraph.o: colorutil.h general.h gustd.h psgraph.h
psgraphf.o: general.h psgraph.h

