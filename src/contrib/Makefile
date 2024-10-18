
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

default: fem

FEMDIR = ../../
include ../../Rules.make

ifeq ($(ECOLOGICAL),BFM)
	compile=compile_bfm
endif

shyfem: nothing

fem: $(compile) nothing

compile_bfm:
	cd ecological/bfm ; make

regress: nothing

status: nothing

list:
	pwd; ls

depend: nothing

clean: nothing

cleanall: nothing

save: nothing

zip: nothing

nothing:
	@true

