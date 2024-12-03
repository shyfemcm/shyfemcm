
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2023  The SHYFEM team
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

DIR	= aquabc

#--------------------------------------------- just for default
FEMDIR = ../../../..
FEMBIN = $(FEMDIR)/bin
LIBDIR = $(FEMDIR)/lib
F77    = g77
FFLAGS =  -O -Wall -g
LINKER = $(F77)
#---------------------------------------------

include $(FEMDIR)/Rules.make

#--------------------------------------------- macros

EXES	=	

ifeq ($(ECOLOGICAL),AQUABC)
  OBJ	     =	\
		aquabc_II_svindex.o \
		aquabc_II_apara.o \
		aquabc_II_ini.o \
		aquabc_II_aout.o \
		aquabc_II.o \
		aquabc_II_co2sys.o \
		aquabc_II_fem_interface.o \
		aquabc_II_pelagic_lib.o \
		aquabc_II_pelagic_model.o \
		aquabc_II_sediment_lib.o \
		aquabc_II_sediment_model_1.o \
		aquabc_II_util.o \

endif

AQUABC_OBJ := $(foreach item,$(OBJ),$(PWD)/$(item))

ALL	=	$(OBJ)

SPECIAL =       Makefile TODO README
SAVEFILES =     *.f90 $(SPECIAL)

MAKEDEPEND = $(FEMBIN)/mkdp.pl

#--------------------------------------------- EXES

default: fem
fem:	$(OBJ) nothing

nothing:
	@true

objs:
	@echo $(AQUABC_OBJ)

#--------------------------------------------- clean

cleanall: clean cleanmod cleanzip cleangrd

clean: cleanobj cleanexe cleanvar cleanlinks

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

cleanvar:
	-rm -f a.out core
	-rm -f *.bak *.tmp
	-rm -f ggg hhh
	-rm -f mon.out nohup.out
	-rm -f fort.*
	-rm -f errout.dat
	-rm -f gmon.out
	-rm -f out.ps plot.ps

cleanmod:

	-rm -f $(LIBDIR)/mod/aquabc_ii_*.mod

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

../../../../lib/mod/aquabc_ii_aout.mod: aquabc_II_aout.o
../../../../lib/mod/aquabc_ii_bsed_model_constants.mod: \
		aquabc_II_sediment_model_1.o 
../../../../lib/mod/aquabc_ii_co2sys_cdiac.mod: aquabc_II_co2sys.o
../../../../lib/mod/aquabc_ii_global.mod: aquabc_II.o
../../../../lib/mod/aquabc_ii_layers_data.mod: aquabc_II_ini.o
../../../../lib/mod/aquabc_ii_para_aqua.mod: aquabc_II_apara.o
../../../../lib/mod/aquabc_ii_pel_state_var_indexes.mod: \
		aquabc_II_svindex.o 
../../../../lib/mod/aquabc_ii_pelagic_model_constants.mod: \
		aquabc_II_pelagic_model.o 
../../../../lib/mod/aquabc_ii_sed_ini.mod: aquabc_II_ini.o
../../../../lib/mod/aquabc_ii_vars.mod: aquabc_II_ini.o
../../../../lib/mod/aquabc_ii_vector_matrix_utils.mod: \
		aquabc_II_co2sys.o 
../../../../lib/mod/aquabc_ii_wc_ini.mod: \
		aquabc_II_pelagic_model.o 
aquabc_II.o: ../../../../lib/mod/aquabc_ii_layers_data.mod \
		../../../../lib/mod/aquabc_ii_para_aqua.mod \
		../../../../lib/mod/aquabc_ii_pel_state_var_indexes.mod \
		../../../../lib/mod/aquabc_ii_sed_ini.mod \
		../../../../lib/mod/basin.mod \
		../../../../lib/mod/levels.mod \
		../../../../lib/mod/mod_debug.mod 
aquabc_II_aout.o: ../../../../lib/mod/aquabc_ii_vars.mod
aquabc_II_co2sys.o: ../../../../lib/mod/aquabc_ii_global.mod
aquabc_II_fem_interface.o: \
		../../../../lib/mod/aquabc_ii_para_aqua.mod \
		../../../../lib/mod/aquabc_ii_sed_ini.mod \
		../../../../lib/mod/aquabc_ii_vars.mod \
		../../../../lib/mod/basin.mod \
		../../../../lib/mod/evgeom.mod \
		../../../../lib/mod/femtime.mod \
		../../../../lib/mod/levels.mod \
		../../../../lib/mod/mod_debug.mod \
		../../../../lib/mod/mod_diff_visc_fric.mod \
		aquabc_II_aout.h 
aquabc_II_ini.o: ../../../../lib/mod/aquabc_ii_para_aqua.mod \
		../../../../lib/mod/aquabc_ii_pel_state_var_indexes.mod \
		../../../../lib/mod/basin.mod \
		../../../../lib/mod/femtime.mod \
		../../../../lib/mod/levels.mod \
		../../../../lib/mod/mod_sediment.mod 
aquabc_II_pelagic_lib.o: ../../../../lib/mod/aquabc_ii_global.mod \
		../../../../lib/mod/mod_debug.mod 
aquabc_II_pelagic_model.o: \
		../../../../lib/mod/aquabc_ii_co2sys_cdiac.mod \
		../../../../lib/mod/aquabc_ii_global.mod \
		../../../../lib/mod/aquabc_ii_para_aqua.mod \
		../../../../lib/mod/aquabc_ii_pel_state_var_indexes.mod \
		../../../../lib/mod/mod_debug.mod 
aquabc_II_sediment_lib.o: ../../../../lib/mod/aquabc_ii_global.mod \
		../../../../lib/mod/mod_debug.mod 
aquabc_II_sediment_model_1.o: \
		../../../../lib/mod/aquabc_ii_co2sys_cdiac.mod \
		../../../../lib/mod/aquabc_ii_global.mod \
		../../../../lib/mod/aquabc_ii_para_aqua.mod \
		../../../../lib/mod/mod_debug.mod 
aquabc_II_util.o: ../../../../lib/mod/basin.mod \
		../../../../lib/mod/femtime.mod \
		../../../../lib/mod/flux.mod \
		../../../../lib/mod/levels.mod \
		../../../../lib/mod/simul.mod 

