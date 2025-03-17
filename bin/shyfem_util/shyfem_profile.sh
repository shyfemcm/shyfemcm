#!/bin/bash
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# shyfem_profile.sh
#
# sets path and aliases
#
# is called from within .bashrc, .bash_profile, .profile
#
#------------------------------------------------------

dist=shyfemcm

FEMDIR=${SHYFEMDIR:=$HOME/$dist}
FEMBIN=$FEMDIR/bin

FEMDIR_INSTALL=${SHYFEM_INSTALL:=$HOME/$dist}
FEMBIN_install=$FEMDIR_INSTALL/bin

# set PATH ----------------------------------------

binutil=$FEMBIN_install
[ -d $FEMBIN_install/shyfem_util ] && binutil=$FEMBIN_install/shyfem_util

path=$( $binutil/shyfem_path.pl $PATH )
export PATH=$path:$FEMBIN

# set aliases ----------------------------------------

alias shyfemdir=". $binutil/shyfem_dir.sh"
alias shyfeminstall=". $binutil/shyfem_dir.sh $FEMDIR_INSTALL"
alias shyfemcd="cd $FEMDIR"
alias shypath=". $binutil/shypath.sh"

# end of routine ----------------------------------------

