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
# compiles with different available compilers
#
# need the follwing routines:
#	check_server.sh
#	subst_make.pl
#
#--------------------------------------------------------

compilers="GNU_GFORTRAN INTEL PGI"
#compilers="GNU_GFORTRAN"
#compilers="INTEL"
#compilers="PGI"

femdir=$( pwd )

rules_arc_dir=$femdir/arc/rules
rules_dist_dir=$femdir/var/rules

rules_save=$rules_arc_dir/Rules.save
rules_dist=$rules_dist_dir/Rules.dist

subst_make=$femdir/bin/subst_make.pl
check_server=$femdir/var/servers/check_server.sh

export FEMDIR=$femdir

debug="YES"
debug="NO"
catch_warnings="NO"
catch_warnings="YES"

#--------------------------------------------------------

#trap Clean_up SIGHUP SIGINT SIGTERM
trap Clean_up 1 2 15

#--------------------------------------------------------

Error()
{
  echo "*** $*"
  exit 1
}

ExitOnError()
{
  [ $? -ne 0 ] && Error $*
}

Clean_up()
{
  echo "Cleaning up before exiting..."
  Clean_after
  exit
}

Clean_before()
{
  rm -f *.out *.tmp
  mkdir -p $rules_arc_dir
  mv --backup=numbered ./Rules.make $rules_save
  cp $rules_dist ./Rules.make
  make cleanall > /dev/null 2>&1
  [ -f allstdout.txt ] && rm allstdout.txt
  [ -f allstderr.txt ] && rm allstderr.txt
}

Clean_after()
{
  rm -f *.tmp
  rm -f stdout.out stderr.out
  #cp ./Rules.make ./Rules.last		#save last Rules.make for inspection
  mv -f $rules_save ./Rules.make
  make cleanall > /dev/null 2>&1
  [ -f allstdout.txt ] && mv allstdout.txt allstdout.tmp
  [ -f allstderr.txt ] && mv allstderr.txt allstderr.tmp
}

SetUp()
{
  . bin/colors.sh

  mkdir -p $rules_arc_dir $rules_dist_dir
  if [ $? -ne 0 ]; then
    echo "Cannot create directory arc... aborting"
    exit 1
  fi
  [ -f $rules_dist ] || cp ./Rules.make $rules_dist
}

WrapUp()
{
  lines=0
  if [ -f allstderr.tmp ]; then
    if [ $catch_warnings = "YES" ]; then
      lines=`cat allstderr.tmp | grep -v 'setting macros' | wc -l`
    else
      lines=`cat allstderr.tmp | grep -i Error | wc -l`
    fi
  fi
  echo "================================="
  echo "Final message on all compilations: "
  if [ $lines -ne 0 ]; then
    echo "${red}  *** some errors occured in compilation...${normal}"
    echo "  (see allstderr.tmp for more details)"
  else
    echo "${green}  no compilation errors${normal}"
  fi
  echo "================================="
}

#--------------------------------------------------------

CompTest()
{
  echo "running CompTest"
  $check_server -show
  Comp "NETCDF=false PARALLEL_OMP=true"
}

CompAll()
{
  host=$( hostname )
  echo "hostname = $host"
  if [ "$host" = "tide" ]; then
    METISDIR=$HOME/georg/lib/metis
  else
    METISDIR=$HOME/lib/metis
  fi

  INTEL_VERSION="IFORT"
  #INTEL_VERSION="IFX"

  Comp "INTEL_VERSION=$INTEL_VERSION \
	ECOLOGICAL=NONE GOTM=true NETCDF=false \
	SOLVER=SPARSKIT \
	PARALLEL_OMP=false PARALLEL_MPI=NONE"
  #Comp "ECOLOGICAL=EUTRO GOTM=false SOLVER=PARDISO"
  Comp "ECOLOGICAL=EUTRO"
  Comp "ECOLOGICAL=NONE MERCURY=true"
  #Comp "ECOLOGICAL=ERSEM GOTM=true NETCDF=true SOLVER=GAUSS"
  Comp "ECOLOGICAL=BFM NETCDF=true SOLVER=SPARSKIT"
  Comp "ECOLOGICAL=AQUABC NETCDF=false PARALLEL_OMP=true"
  Comp "ECOLOGICAL=NONE GOTM=true NETCDF=false SOLVER=SPARSKIT \
	PARALLEL_OMP=false PARALLEL_MPI=NODE \
	PARTS=METIS METISDIR=$METISDIR"
  Comp "SOLVER=PETSC"

  [ "$regress" = "NO" ] && return

  Rules "ECOLOGICAL=NONE GOTM=true NETCDF=false SOLVER=SPARSKIT \
		PARALLEL_OMP=false PARALLEL_MPI=NONE"

  Comp "COMPILER_PROFILE=SPEED PARALLEL_OMP=true"
  Comp "COMPILER_PROFILE=CHECK PARALLEL_OMP=false"
}

Comp()
{
  echo "-----------------"
  Rules "$1"

  echo "start compiling in" `pwd`
  rm -f stdout.out stderr.out
  touch stdout.out stderr.out

  if [ "$dry_run" = "NO" ]; then
    make cleanall > tmp.tmp 2> tmp.tmp
    make fem > stdout.out 2> stderr.tmp
  fi

  [ -f stderr.tmp ] && cat stderr.tmp | grep -v "ar: creating" > stderr.out

  lines=0
  if [ $catch_warnings = "YES" ]; then
    lines=`cat stderr.out | wc -l`	# also catches warnings
  else
    lines=`cat stderr.out | grep -i Error | wc -l`
  fi

  if [ $lines -ne 0 ]; then
    echo "${red}*** compilation errors...${normal}"
    cat stderr.out
  else
    echo "${green}successfull compilation${normal}"
  fi
  echo "-----------------"

  cat stdout.out >> allstdout.txt
  cat stderr.out >> allstderr.txt

  Regress
}

RulesReset()
{
  # resets variables in Rules.make

  echo "resetting Rules.make file..."
  cp $rules_dist ./Rules.make
}

Rules()
{
  # sets variables in Rules.make

  $subst_make -quiet -first "$1" Rules.make > tmp.tmp
  ExitOnError "error executing subst_make.pl ...aborting"
  mv tmp.tmp Rules.make

  echo "setting macros: $1"
  echo "setting macros: $1" >> allstdout.txt
  echo "setting macros: $1" >> allstderr.txt
}

PrintAll()
{
  echo "setting macros: $1"
  echo "setting macros: $1" >> allstdout.txt
  echo "setting macros: $1" >> allstderr.txt
}

Regress()
{
  [ "$regress" = "NO" ] && return

  echo "running regression test..."
  make regress >> allstdout.txt
  cd femregress
  make status
  cd ..
  echo "finished running regression test..."
}

SetCompiler()
{
  local comp=$1

  local compiler
  local basedir=$femdir/var/servers
  local script=$check_server
  [ ! -x $script ] && Error "script not existing: $script"
  local server=$( $script -server )

  [ "$debug" = "YES" ] && echo "SetCompiler: $comp $hostname"

  if [ "$comp" = "GNU_GFORTRAN" ]; then
    compiler=gfortran
  elif [ "$comp" = "INTEL" ]; then
    compiler=intel
  elif [ "$comp" = "PGI" ]; then
    compiler=pgi
  else
    echo "*** (SetCompiler) no such compiler: $comp"
    return 1
  fi

  [ "$debug" = "YES" ] && echo "SetCompiler: $compiler"

  if [ -n "$server" ]; then
    echo "executing script $script to load settings"
    source $basedir/$script -load $compiler
  fi
}

#--------------------------------------------------------------------
# start routine
#--------------------------------------------------------------------

if [ ! -f VERSION ]; then
  echo "must be in base directory to run this script... aborting"
  exit 1
fi

regress="NO"
[ "$1" = "-regress" ] && regress="YES"
dry_run="NO"
[ "$1" = "-dry_run" ] && dry_run="YES"

SetUp
Clean_before

for comp in $compilers
do

  PrintAll "================================="
  PrintAll "compiling with $comp"
  PrintAll "================================="

  SetCompiler $comp
  [ $? -ne 0 ] && continue

  RulesReset
  Rules "FORTRAN_COMPILER=$comp"

  make compiler_version > /dev/null 2>&1
  [ $? -ne 0 ] && echo "*** compiler $comp is not available..." && continue

  #CompTest
  #continue

  CompAll
done

Clean_after
WrapUp

#--------------------------------------------------------------------
# end of routine
#--------------------------------------------------------------------

