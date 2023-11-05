
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2023  The SHYFEM team
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

#--------------------------------------------- 

DIR	= util

#--------------------------------------------- just for default
FEMDIR = ../..
F77    = g77
FFLAGS =  -O -Wall -g
LINKER = $(F77)
#--------------------------------------------- include Rules.make

include $(FEMDIR)/Rules.make

#--------------------------------------------- libraries

LIBNAME = $(DIR)
LIBDIR  = $(FEMDIR)/lib
INCDIR = $(FEMDIR)/include
INCLUDE = -I$(INCDIR)

#--------------------------------------------- macros

MOD  = \
	para.o \

OBJ  = \
	chk_NaN.o \
	clo.o \
	convert.o \
	dts.o \
	file.o \
	geometrical.o \
	intp.o \
	iso8601.o \
	mat_util.o \
	nls.o \
	projection.o \
	random.o \
	random_gen.o \
	round.o \
	grd.o \
	strings.o \
	subsort.o \
	sort_util.o \
	util.o \
	utility.o \
	vec_util.o \
	femfile.o \
	color.o \
	to_linear.o \
	dates.o \
	unitutils.o \
	tsfile.o \
	heat_util.o \
	expand.o \
	connections.o \
	time_util.o \
	vert_intp.o \
	int_to_cart.o \
	bit10.o \
	nodelist.o \

ALL = $(MOD) $(OBJ)

SPECIAL =       Makefile

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

LIBL =
LIBF =

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

check:
	~/shyfem/fem3d/bin/compile.sh -single *.f90
	~/shyfem/fem3d/bin/compile.sh *.f90
	~/shyfem/fem3d/bin/usemod.pl *.f90

check_utils:
	@~/shyfem/fem3d/bin/check_utils.sh $(ALL)

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
	$(F77) -c $(INCLUDE) $(FFLAGS) $<

.F.o:
	$(F77) -c -cpp $(FFLAGS) $<

.for.obj:
	$(FLC) $< , $(<,B).OBJ $(FLFLAGS)

#--------------------------------------------- makedepend

# DO NOT DELETE THIS LINE -- make depend depends on it.

