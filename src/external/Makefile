
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2023  The SHYFEM team
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

FEMDIR = ../..
BINDIR = $(FEMDIR)/bin

SUBDIRS = gotm sparskit

#------------------------------------------------------------------------

default: shyfem

fem: shyfem
shyfem: nothing
	@$(BINDIR)/recursivemake $@ $(SUBDIRS)

test: nothing
	@$(BINDIR)/recursivemake $@ shynetcdf

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

