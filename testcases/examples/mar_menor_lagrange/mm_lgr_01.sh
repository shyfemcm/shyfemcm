#!/bin/sh
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
#--------------------------------------------------
. ./util.sh
#--------------------------------------------------

sim=mm_lgr_01

CleanFiles $sim.hydro.shy 

Run $sim
freq=8

CheckFiles $sim.hydro.shy
PlotMapVel apn_vel
PlotMapVel apnbath

#--------------------------------------------------

