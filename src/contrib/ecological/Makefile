
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

FEMDIR = ../../..

BINDIR = $(FEMDIR)/bin
SUBDIRS = bfm mercury aquabc weutro

#------------------------------------------------------------------------

default: fem

shyfem: fem

fem: nothing
	@$(BINDIR)/recursivemake $@ $(SUBDIRS)

test: nothing
	@$(BINDIR)/recursivemake $@ $(SUBDIRS)

regress: nothing

status: nothing

list:
	pwd; ls

depend: nothing
	@$(BINDIR)/recursivemake $@ $(SUBDIRS)

clean: nothing
	@$(BINDIR)/recursivemake $@ $(SUBDIRS)

cleanall: nothing
	@$(BINDIR)/recursivemake $@ $(SUBDIRS)

save: nothing

zip: nothing

nothing:
	@true

#------------------------------------------------------------------------

