
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

#--------------------------------------------- 

DIR	= shyutil

#--------------------------------------------- just for default
FEMDIR = ../../..
FEMBIN = $(FEMDIR)/bin
F77    = g77
FFLAGS =  -O -Wall -g
LINKER = $(F77)
#--------------------------------------------- include Rules.make

include $(FEMDIR)/Rules.make

#--------------------------------------------- libraries

LIBNAME = $(DIR)
LIBDIR  = $(FEMDIR)/lib

INCLUDE = -I$(FEMDIR)/include
FFLAGS += $(INCLUDE)

#--------------------------------------------- macros

MOD  = \

OBJ  = \
	version.o \
	basin.o \
	depth.o \
	grdutil.o \
	strings.o \
	fem_util.o \
	femfile_util.o \
	defnames.o \
	geom_dynamic.o \
	regular.o \
	expand.o \
	geom.o \
	levels.o \
	links.o \
	lnku.o \
	bnd.o \
	bnd_aux.o \
	bnd_dynamic.o \
	bnd_geom.o \
	sigma.o \
	line_util.o \
	simul.o \
	io_shy.o \
	io_ext.o \
	io_flx.o \
	io_lgr.o \
	io_off.o \
	zadaptutil.o \
	submask.o \
	topological.o \
	pkonst.o \
	mkonst.o \
	heat_const.o \

ALL = $(MOD) $(OBJ)

SPECIAL =       Makefile

MAKEDEPEND = $(FEMBIN)/mkdp.pl

#--------------------------------------------- EXES

LIB = $(LIBDIR)/lib$(LIBNAME).a

default: lib
all:	lib
fem: lib
lib:	$(LIB) nothing

nothing:
	@true

$(LIB): $(ALL)
	ar rvs $(LIB) $(ALL)

LIBUTIL     = util

LIBLUT  = -L$(LIBDIR) -l$(LIBUTIL)
LIBFUT  = $(LIBDIR)/lib$(LIBUTIL).a

LIBL = $(LIBLUT)
LIBF = $(LIBFUT)

#--------------------------------------------- test compile

test: lib main_dummy

main_dummy:
	@echo " write(6,*) 'all dependencies ok'\n end" > main_dummy.f90
	$(F77) $(FFLAGS) -c $@.f90 
	$(F77) $(FFLAGS) $(INCLUDE) $(LFLAGS) $(ALL) $(LIBL) $@.o -o $@
	@echo "...running main_dummy"
	@./main_dummy
	@-rm -f main_dummy main_dummy.f90 main_dummy.o

#--------------------------------------------- clean

cleanall: clean cleanlib

clean: cleanobj cleanvar

cleanobj:
	-rm -f *.o
	-rm -f *.mod
	-rm -f *_genmod.f90

cleanvar:
	-rm -f main_dummy main_dummy.f90
	-rm -f a.out core
	-rm -f *.bak *.tmp
	-rm -f ggg hhh
	-rm -f mon.out nohup.out
	-rm -f fort.*
	-rm -f gmon.out
	-rm -f out.ps plot.ps

cleanlib:
	-rm -f $(LIB)

#--------------------------------------------- various

depend:
	@-$(MAKEDEPEND) -moddir=$(MODDIR) *.f90

check:
	~/shyfem/fem3d/bin/compile.sh -single *.f90
	~/shyfem/fem3d/bin/usemod.pl *.f90

list:
	pwd; ls

save:   zip
zip:	cleanall
	-zip $(DIR) $(SAVEFILES)

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

../../../lib/mod/basin.mod: basin.o
../../../lib/mod/default_names.mod: defnames.o
../../../lib/mod/extfile.mod: io_ext.o
../../../lib/mod/fem_util.mod: femfile_util.o
../../../lib/mod/flxfile.mod: io_flx.o
../../../lib/mod/heat_const.mod: heat_const.o
../../../lib/mod/levels.mod: levels.o
../../../lib/mod/link_structure.mod: links.o
../../../lib/mod/mkonst.mod: mkonst.o
../../../lib/mod/mod_bnd.mod: bnd.o
../../../lib/mod/mod_bnd_aux.mod: bnd_aux.o
../../../lib/mod/mod_bound_dynamic.mod: bnd_dynamic.o
../../../lib/mod/mod_bound_geom.mod: bnd_geom.o
../../../lib/mod/mod_depth.mod: depth.o
../../../lib/mod/mod_geom.mod: geom.o
../../../lib/mod/mod_geom_dynamic.mod: geom_dynamic.o
../../../lib/mod/mod_offline.mod: io_off.o
../../../lib/mod/pkonst.mod: pkonst.o
../../../lib/mod/regular.mod: regular.o
../../../lib/mod/shyfem_strings.mod: strings.o
../../../lib/mod/shyfem_version.mod: version.o
../../../lib/mod/shyfile.mod: io_shy.o
../../../lib/mod/sigma.mod: sigma.o
../../../lib/mod/simul.mod: simul.o
../../../lib/mod/zadapt.mod: zadaptutil.o
fem_util.o: ../../../lib/mod/basin.mod
grdutil.o: ../../../lib/mod/basin.mod ../../../lib/mod/grd.mod
line_util.o: ../../../lib/mod/basin.mod
links.o: ../../../lib/mod/basin.mod ../../../lib/mod/mod_geom.mod
lnku.o: ../../../lib/mod/basin.mod ../../../lib/mod/mod_geom.mod
regular.o: ../../../lib/mod/basin.mod
submask.o: ../../../lib/mod/basin.mod
topological.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/mod_geom.mod 
zadaptutil.o: ../../../lib/mod/basin.mod \
		../../../lib/mod/levels.mod 

