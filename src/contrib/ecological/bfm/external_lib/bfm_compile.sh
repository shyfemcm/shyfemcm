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

RemoveExtension()
{
  local dir=$1
  local ext=$2
  local actdir=$( pwd )

  [ $# -lt 2 ] && return

  if [ -d $dir ]; then
    [ $verbose = "YES" ] && echo "removing files with extension $ext in $dir"
    cd $dir
    #ls  *.$ext
    rm -f *.$ext
    cd $actdir
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

DeletFromLib()
{
  local object=$1
  local lib=$2

  [ $# -lt 2 ] && return

  if [ -f $lib ]; then
    [ $verbose = "YES" ] && echo "removing $object from $lib"
    ar -dv $lib $object
    echo ""
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

FORTRAN_COMPILER=$2
if [ -z "$FORTRAN_COMPILER" -a $clean = "NO" ]; then
  echo "*** FORTRAN_COMPILER is not given - please specify in Rules.make"
  echo "    Please see Rules.make for more details"
  exit 7
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
  RemoveExtension include mod
  RemoveFile include/BFM_module_list.h
  RemoveFile include/BFM_var_list.h
  exit 0
fi

shydir=$( pwd )
echo "using BFM directory: $BFMDIR"
echo "this creates the BFM library to be linked with SHYFEM"
echo "this has to be done only once"
echo "you also have to recompile if you change compiler"

#---------------------------------------------------------------
# prepare compilation parameters
#---------------------------------------------------------------

PL_DEBUG="YES"
PL_DEBUG="NO"		# this is for production runs
PL_DEBUG_EXT=""

if [ $FORTRAN_COMPILER = "GNU_GFORTRAN" ]; then
  PL_COMPILER=gfortran
  [ $PL_DEBUG = "YES" ] && PL_DEBUG_EXT="_debug"
elif [ $FORTRAN_COMPILER = "INTEL" ]; then
  PL_COMPILER=intel
  [ $PL_DEBUG = "YES" ] && PL_DEBUG_EXT=".dbg"
else
  echo "fortran compiler not recognized: $FORTRAN_COMPILER"
  exit 9
fi

PL_ARCH=x86_64
PL_OS=LINUX

if [ $PL_COMPILER = "intel" ]; then
  INC_FILE=${PL_ARCH}.${PL_OS}.${PL_COMPILER}${PL_DEBUG_EXT}.inc
elif [ $PL_COMPILER = "gfortran" ]; then
  INC_FILE=${PL_COMPILER}${PL_DEBUG_EXT}.inc
fi

echo "using .inc file: $INC_FILE"

#export MODULEFILE=/g100_work/OGS23_PRACE_IT/SHYFEM_BFM/shyfem-bfm/v.00.veg/shyfem-bfm-sgr.v00.config

#source $MODULEFILE

#---------------------------------------------------------------
# compile
#---------------------------------------------------------------

cd $BFMDIR

sed -i "s/.*ARCH.*/        ARCH    = '$INC_FILE'  /"  \
		build/configurations/OGS_PELAGIC/configuration

mkdir -p $BFMDIR/lib
mkdir -p $BFMDIR/build/tmp/OGS_PELAGIC
echo "copying BFM_module_list.proto.h to $BFMDIR/build/tmp/OGS_PELAGIC"
cp include/BFM_module_list.proto.h $BFMDIR/build/tmp/OGS_PELAGIC

cd $BFMDIR/build

./bfm_configure.sh -gcv -o ../lib/libbfm.a -p OGS_PELAGIC

if [ $? -ne 0 ] ; then  echo  ERROR in code generation; exit 1 ; fi

echo "deleting standalone_main.o from library"
DeletFromLib standalone_main.o $BFMDIR/lib/libbfm.a

echo "...the library is ready in $BFMDIR/lib"

#---------------------------------------------------------------
# end of routine
#---------------------------------------------------------------

