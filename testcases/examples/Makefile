
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------


DIR = examples

FEMDIR = ..
FEMBIN = ../bin

SUBDIRS = mar_menor

#-----------------------------------------------------------

default: nothing

fem: nothing

nothing:
	@true

list:
	@pwd; ls

depend: nothing

clean:
	$(FEMBIN)/recursivemake $@ $(SUBDIRS)

cleanall:
	$(FEMBIN)/recursivemake $@ $(SUBDIRS)

#-----------------------------------------------------------

