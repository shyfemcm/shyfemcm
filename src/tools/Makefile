
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2023  The SHYFEM team
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

FEMDIR = ../..
BINDIR = $(FEMDIR)/bin

SUBDIRS = \
	grid mesh postscript hcbs \
	shyadj shybas \
	shynetcdf shyelab \
	shyparts shyplot shypre \

#------------------------------------------------------------------------

default: fem

fem: nothing
	@$(BINDIR)/recursivemake $@ $(SUBDIRS)

test: nothing
	@$(BINDIR)/recursivemake $@ shynetcdf

regress: nothing

status: nothing

list:
	pwd; ls

depend: nothing

clean: nothing
	@$(BINDIR)/recursivemake $@ $(SUBDIRS)

cleanall: nothing
	@$(BINDIR)/recursivemake $@ $(SUBDIRS)

save: nothing

zip: nothing

nothing:
	@true

