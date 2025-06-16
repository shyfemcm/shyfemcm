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
# shyfem_dir.sh
#
# resets path and and environment variables
#
# called without argument returns info on actual SHYFEM distribution
# called with argument sets the new shyfem directory
#
# in order to work this script must be called through shyfemdir
# or must be sourced: ". shyfem_dir.sh" or "source shyfem_dir.sh"
#
# Usage: 
#    . shyfem_dir.sh [femdir]
#    source shyfem_dir.sh [femdir]
#    shyfemdir [femdir]
#
#------------------------------------------------------

dist=shyfemcm

FEMDIR=${SHYFEMDIR:=$HOME/$dist}
FEMBIN=$FEMDIR/bin

#FEMDIR_INSTALL=${SHYFEM_INSTALL:=$HOME/$dist}
#FEMDIR_INSTALL=${SHYFEM_INSTALL:=$FEMDIR}
FEMDIR_INSTALL=${FEMDIR}
FEMBIN_install=$FEMDIR_INSTALL/bin

GetVersion()
{
  local dir=$1

  if [ ! -x $shyutil/shyfem_version.sh ]; then
    echo "*** cannot find shyfem_version.sh ...aborting" >&2
    version=""
  else
    version=$( $shyutil/shyfem_version.sh $dir )
  fi

  echo $version
}

shyutil=$FEMBIN_install
if [ -d $FEMBIN_install/shyfem_util ]; then
  shyutil=$FEMBIN_install/shyfem_util
fi

# command line options ----------------------------

write="normal"
if [ "$1" = "-quiet" ]; then
  write="quiet"
  shift
elif [ "$1" = "-debug" ]; then
  write="debug"
  shift
fi

# change path and environment variables ---------------

[ $write = "debug" ] && echo "debug message: FEMDIR = $FEMDIR"

if [ -n "$1" ]; then

  dir=`readlink -f $1`	# get full path name

  [ $write = "debug" ] && echo "debug message: using dir as $dir"

  version=$( GetVersion $dir )
  if [ -z "$version" -o "$version" = "unknown" ]; then
    echo "*** cannot get version for $dir ... aborting" 1>&2
  else
    export SHYFEMDIR=$dir
    FEMDIR=${SHYFEMDIR:=$HOME/$dist}
    FEMBIN=$FEMDIR/bin

    path=$( $shyutil/shyfem_path.pl $PATH )
    export PATH=$FEMBIN:$path
  fi

fi

# show path and environment variables ---------------

version=$( GetVersion )

if [ $write != "quiet" ]; then

  echo "SHYFEM version: $version"
  echo "SHYFEM directory: $FEMDIR"
  echo "SHYFEM install directory: $FEMDIR_INSTALL"

  if [ $write = "debug" ]; then
    echo "SHYFEM path: $PATH"
  fi
fi

# end of routine ----------------------------------------

