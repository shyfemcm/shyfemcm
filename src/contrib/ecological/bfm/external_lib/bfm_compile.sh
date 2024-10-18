#!/bin/bash

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2018  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

#---------------------------------------------------------------
#
# this script sets up the BFM model for SHYFEM
#
# this script has to be executed in the main directory of SHYFEM
#
#---------------------------------------------------------------

CleanDir()
{
  local dir=$1

  [ $# -eq 0 ] && return

  if [ -d "$dir" ]; then
    echo "cleaning dir $dir"
    rm -rf $dir/*
  fi
}

RemoveDir()
{
  local dir=$1

  [ $# -eq 0 ] && return

  if [ -d "$dir" ]; then
    [ $verbose = "YES" ] && echo "removing dir $dir"
    rm -rf $dir
  fi
}

RemoveFile()
{
  local file=$1

  [ $# -eq 0 ] && return

  if [ -f "$file" ]; then
    [ $verbose = "YES" ] && echo "removing file $file"
    rm -f $file
  fi
}

RestoreSave()
{
  local file=$1

  [ $# -eq 0 ] && return

  if [ -f $file.save ]; then
    [ $verbose = "YES" ] && echo "restoring $file"
    mv -f $file.save $file
  fi
  rm -f $file.*
}

RestoreSavesDir()
{
  local dir=$1

  [ $# -eq 0 ] && return

  local actdir=$( pwd )
  cd $dir

  for file in *.save
  do
    name=$( basename $file .save )
    RestoreSave $name
  done

  cd $actdir
}

SaveFile()
{
  local file=$1

  [ $# -eq 0 ] && return

  if [ ! -f $file.save ]; then
    [ $verbose = "YES" ] && echo "saving $file"
    cp $file $file.save
  fi
}

SetNetcdf()
{
  # this is not yet used, but maybe needed with new BFM version

  netcdfdir=$shydir/var/netcdf
  netcdfscript=$netcdfdir/set_netcdf.sh

  export NETCDFF_LIB=$( $netcdfscript libnetcdff )
  export NETCDF_LIB=$( $netcdfscript libnetcdf )
  export NETCDFF_INC=$( $netcdfscript netcdf.mod )

  echo "NETCDFF_LIB $NETCDFF_LIB"
  echo "NETCDF_LIB  $NETCDF_LIB"
  echo "NETCDFF_INC $NETCDFF_INC"
}

#---------------------------------------------------------------

verbose="NO"

if [ ! -f VERSION ]; then
  echo "*** you must be in the main directory of SHYFEM to run this script"
  exit 1
fi

clean="NO"
if [ "$1" = "-clean" ]; then	# this cleans BFMDIR
  clean="YES"
  shift
fi

BFMDIR=$1
if [ -z "$BFMDIR" ]; then
  echo "*** BFMDIR is not given - please specify in Rules.make"
  echo "    Please see Rules.make for more details"
  exit 3
fi

if [ ! -d $BFMDIR ]; then
  echo "*** BFMDIR is not exisiting - please correct in Rules.make"
  exit 5
fi

if [ $clean = "YES" ]; then
  verbose="YES"
  echo "cleaning BFM directory: $BFMDIR"
  cd $BFMDIR
  CleanDir build/tmp
  CleanDir build/Logs
  CleanDir bin
  CleanDir run
  RemoveDir lib
  RemoveDir build/configurations/BFMLIB
  RemoveDir src/shyfem
  RestoreSave build/bfm_configure.sh
  RestoreSavesDir compilers
  RemoveFile compilers/shyfem.extra
  RemoveFile src/share/init_var_bfm.F90
  RemoveFile src/BFM/include/INCLUDE.h
  RemoveFile src/BFM/General/AllocateMem.F90
  RemoveFile src/BFM/General/ModuleMem.F90
  RemoveFile src/BFM/General/set_var_info_bfm.F90
  exit 0
fi

shydir=$( pwd )
echo "using BFM directory: $BFMDIR"
echo "this creates the BFM library to be linked with SHYFEM"
echo "this has to be done only once"
echo "you also have to recompile if you change compiler"

#---------------------------------------------------------------
# compile
#---------------------------------------------------------------
PL_ARCH=x86_64
PL_OS=LINUX
PL_COMPILER=intel
PL_DEBUG=       # this is the choice for production flags 
#PL_DEBUG=.dbg   # this is the one for debug flags
INC_FILE=${PL_ARCH}.${PL_OS}.${PL_COMPILER}${PL_DEBUG}.inc
#export MODULEFILE=/g100_work/OGS23_PRACE_IT/SHYFEM_BFM/shyfem-bfm/v.00.veg/shyfem-bfm-sgr.v00.config

#source $MODULEFILE

cd $BFMDIR

sed -i "s/.*ARCH.*/        ARCH    = '$INC_FILE'  /"  build/configurations/OGS_PELAGIC/configuration

cd $BFMDIR/build

mkdir -p $BFMDIR/lib

./bfm_configure.sh -gcv -o ../lib/libbfm.a -p OGS_PELAGIC

if [ $? -ne 0 ] ; then  echo  ERROR in code generation; exit 1 ; fi

echo "...the library is ready in $BFMDIR/lib"

#---------------------------------------------------------------
# end of routine
#---------------------------------------------------------------

