
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

FEMDIR = ..
BINDIR = $(FEMDIR)/bin
SUBDIRS = \
		grid mesh postscript hcbs grdadj \
		shynetcdf shypre shybas shyelab \
		shyplot

#------------------------------------------------------------------------

default: fem

fem: nothing
	$(BINDIR)/recursivemake $@ $(SUBDIRS)

test: nothing
	$(BINDIR)/recursivemake $@ shynetcdf

regress: nothing

status: nothing

list:
	pwd; ls

depend: nothing

clean: nothing
	$(BINDIR)/recursivemake $@ $(SUBDIRS)

cleanall: nothing
	$(BINDIR)/recursivemake $@ $(SUBDIRS)

save: nothing

zip: nothing

nothing:
	@true

