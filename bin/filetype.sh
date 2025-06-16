#!/bin/sh

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

FEMDIR=${SHYFEMDIR:=$HOME/shyfem}

elabdir=$FEMDIR/src/tools/shyelab

$elabdir/filetype $*


