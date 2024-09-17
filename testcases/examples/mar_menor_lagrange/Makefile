
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------


DIR = mar_menor_lagrange
BASE_DIR = mar_menor

GRID  = MM_coarse
GRID1  = MM_coarse_smooth

#-----------------------------------------------------------

default:
	@echo "possible targets:"
	@echo "  init      initializes the bas file"
	@echo "  run_all   this runs all the simulations"
	@echo "  clean     cleans directory"
	@echo "  cleanall  cleans everyting in directory"

fem: nothing

nothing:
	@true

list:
	pwd; ls

depend: nothing

#-----------------------------------------------------------

init: links
	cd GRID; shypre $(GRID).grd
	cd GRID; shypre $(GRID1).grd

links:
	ln -fs ../$(BASE_DIR)/GRID
	ln -fs ../$(BASE_DIR)/INPUT
	ln -fs ../$(BASE_DIR)/util.sh
	ln -fs ../$(BASE_DIR)/apn_vel.str
	ln -fs ../$(BASE_DIR)/apnbath.str

run_all: init
	./make_all_sims.sh

#-----------------------------------------------------------

clean:
	-rm -f rkdiff.nos dist.nos bndo_info.dat
	-rm -f *.tmp *.bak
	-rm -f *.[0123]d.[0-9]*
	-rm -f *.[23]d.txt
	-rm -f out*.ps plot*.ps
	-rm -f zcomp.2d.*
	-rm -f out.fem out.ext out.flx
	-rm -f fort.*
	-rm -f errout.dat gmon.out
	-rm -f .memory

cleanall: clean cleansims cleanzip cleanfigs cleanbas cleanlinks

cleanbas:
	-rm -f GRID/*.bas
	-rm -f GRID/errout.dat

cleanfigs:
	-rm -f *.ps
	-rm -f *.eps
	-rm -f *.jpg

cleansims:
	-rm -f mm_lgr*.shy
	-rm -f mm_lgr*.inf
	-rm -f mm_lgr*.ext
	-rm -f mm_lgr*.flx
	-rm -f out.shy

cleanzip:
	-rm -f $(DIR).zip

cleanlinks:
	-rm -f GRID INPUT util.sh
	-rm -f apn_vel.str apnbath.str

#-----------------------------------------------------------

zip:
	#-zip -r $(DIR).zip *
	-zip  $(DIR).zip Makefile *.sh *.str GRID/* INPUT/*

#-----------------------------------------------------------

