#!/bin/sh

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

log=CHECKLOG.tmp
rm -f $log

missing=""

#---------------------------------------------------

CheckCommand()
{
  name=$1
  command=$2
  error=$3
  options=$4

  if [ -z "$error" ]; then
    error=0
  fi

  ($command) >> $log 2>&1

  status=$?

  if [ $status -eq $error ]; then
    echo "... $name is ${green}installed${normal}"
  else
    if [ "$options" != "quiet" ]; then
      echo "*** $name is ${flag}not installed${normal}"
    fi
    missing="$missing $name"
  fi
}

CreateInputFiles()
{

cat > test_nc.f <<EOI
	include 'netcdf.inc'
        integer ncid,retval
        retval = nf_create('out.nc', nf_clobber, ncid)
        write(6,*) 'Error code: ', nf_strerror(retval)
        end
EOI

cat > test_x11.c <<EOI
#include <stdio.h>
#include <X11/Xlib.h>

int main( void )
{
  printf("Hello world.\n");
  return 0;
}
EOI

echo "quit" > quit.tmp
}

CleanUp()
{
	[ -f test_nc.f ] && rm -f test_nc.f
	[ -f test_x11.c ] && rm -f test_x11.c
	[ -f quit.tmp ] && rm -f quit.tmp
	[ -f a.out ] && rm -f a.out
}

Exists() 
{
  command -v "$1" >/dev/null 2>&1
}

#---------------------------------------------------

CheckMpiCompiler()
{
  local missing_save=$missing

  CheckCommand mpi_compiler "mpif90 -v" "" "quiet"
  [ $status -eq 0 ] && mpi_available="$mpi_available mpif90"

  if [ -n "$mpi_available" ]; then
    echo "... the following Mpi compilers are available:"
    echo "          ${green}$mpi_available${normal}"
  else
    echo "*** ${red}No Mpi compiler found${normal}"
    echo "    ... please install a Mpi compiler if you want to run in MPI"
    echo "    (on debian the packages may be: openmpi-bin libopenmpi-dev)"
    missing_save="$missing_save mpif90"
  fi
}

CheckFortranCompiler()
{
  #local fortran_available=""
  local missing_save=$missing

  #CheckCommand g77 "g77 -v" "" "quiet"
  #[ $status -eq 0 ] && fortran_available="$fortran_available g77"
  CheckCommand gfortran "gfortran -v" "" "quiet"
  [ $status -eq 0 ] && fortran_available="$fortran_available gfortran"
  CheckCommand ifort "ifort -v" "" "quiet"
  [ $status -eq 0 ] && fortran_available="$fortran_available ifort"

  #fortran_available=""		# fake error

  if [ -n "$fortran_available" ]; then
    echo "... the following Fortran compilers are available:"
    echo "          ${green}$fortran_available${normal}"
    echo "... please set COMPILER in Rules.make to one of the compilers above"
  else
    echo "*** ${red}No Fortran compiler found${normal}"
    echo "    ... please install a fortran compiler"
    echo "    (on Linux you can easily install gfortran)"
    missing_save="$missing_save gfortran ifort"
  fi

  fortran=`echo $fortran_available | sed -e 's/ .*//'`	#first fortran found
  missing=$missing_save
}

CheckX11()
{
  CheckCommand "X11-develop" \
		"gcc -L/usr/X11/lib -L/usr/X11R6/lib -lXt -lX11 test_x11.c"

  if [ $status -ne 0 ]; then
    RecommendPackageFull "X11 development package" \
	libx11 libx11-common libx11-dev libxt-dev x11proto-core-dev
    Comment "Please note that you can still compile the model without graphics"
    Comment "by running \"make nograph\" instead of running \"make fem\""
  fi
}

CheckNetcdf()
{
  netcdf=`GetMacro NETCDF`
  netcdfdir=`GetMacro NETCDFDIR`

  local debug=1
  local nfconfig=0
  local ncdump=0

  #echo "netcdf: $netcdf $netcdfdir    $fortran"

  [ -z "$fortran" ] && return		# no fortran compiler
  #[ "$netcdf" != "true" ] && return	# no netcdf requested
  Exists nf-config && nfconfig=1
  Exists ncdump && ncdump=1

  if [ $debug -eq 1 ]; then
    echo "checking netcdf: $nfconfig $ncdump"
    ncdump=0
  fi

  if [ -f $netcdfdir/lib/libnetcdff.a ]; then
    netcdflib=$netcdfdir/lib
  elif [ -f $netcdfdir/lib/x86_64-linux-gnu/libnetcdff.a ]; then
    netcdflib=$netcdfdir/lib/x86_64-linux-gnu/
  else
    netcdflib=$netcdfdir
  fi

  if [ $nfconfig -eq 1 ]; then
    CheckCommand netcdf \
	"$fortran test_nc.f  $(nf-config --flibs) $(nf-config --fflags)"
  elif [ $ncdump -eq 1 ]; then
    CheckCommand netcdf \
	"$fortran -L$netcdflib -I$netcdfdir/include -lnetcdff test_nc.f"
  else
    CheckCommand netcdf false
  fi

  if [ $status -ne 0 ]; then
    RecommendPackageFull "netcdf package" \
	libnetcdf-dev libnetcdff-dev netcdf-bin
    echo "      If you have installed netcdf maybe the libraries cannot be found."
    echo "      Please also try one of the following commands:"
    echo "        ldconfig -p | grep libnetcdff"
    echo "        whereis libnetcdff"
    echo "        locate libnetcdff.a"

  fi
}

Comment()
{
    echo "      $1"
}

RecommendPackageFull()
{
  what=$1
  shift

  echo "*** $what seems not to be installed"
  echo "      Names of the software packages are not always standard"
  echo "      please try with the following: "
  echo "        $*"
  echo "      try to install one package at a time and then check"
  echo "      the status of the installation again: make check_software"
}

GetMacro()	# gets macro from Rules.make file
{
  what=$1
  rules="../Rules.make"
  [ $# -gt 1 ] && rules="$2"

  macro=$( cat $rules | sed -e 's/ *//g' | grep "^$what=" \
	| tail -1 | sed -e 's/.*=//')

  echo "$macro"
}

SetColors()
{
  ncolors=$(tput colors)
  if test -n "$ncolors" && test $ncolors -ge 8; then
        bold="$(tput bold)"
        underline="$(tput smul)"
        standout="$(tput smso)"
        normal="$(tput sgr0)"
        black="$(tput setaf 0)"
        red="$(tput setaf 1)"
        green="$(tput setaf 2)"
        yellow="$(tput setaf 3)"
        blue="$(tput setaf 4)"
        magenta="$(tput setaf 5)"
        cyan="$(tput setaf 6)"
        white="$(tput setaf 7)"
  fi
}

#---------------------------------------------------

SetColors
CreateInputFiles
flag=$red

echo
echo "... ${bold}checking base applications... (needed)${normal}"

CheckCommand make "make -v"
CheckCommand bash "bash --version"
CheckCommand perl "perl -v"

echo
echo "... ${bold}checking Fortran compilers (needed)${normal}"

CheckFortranCompiler
CheckMpiCompiler

echo
echo "... ${bold}checking c compiler and X11 (needed)${normal}"

CheckCommand gcc "gcc -v"
CheckCommand g++ "g++ -v"
CheckX11

flag=$magenta

echo
echo "... ${bold}checking commodity programs (recommended)${normal}"
CheckCommand dialog  "dialog --version"

echo
echo "... ${bold}checking graphical routines (recommended)${normal}"

CheckCommand "ghostview (gv)" "gv --version"
CheckCommand "ghostscript (gs)" "gs -v"
CheckCommand gnuplot "gnuplot quit.tmp"
CheckCommand "ImageMagick (imagemagick)" "mogrify -version"
CheckCommand gifsicle "gifsicle --version"
CheckNetcdf

echo
echo "... ${bold}checking additional routines (not urgently needed)${normal}"

CheckCommand "latex (texlive)" "latex -v"
CheckCommand dvips "dvips -v"
CheckCommand python "python -V"

CheckCommand ssh "ssh -V"
CheckCommand at "at -l"

#CheckCommand dnotify "dnotify --version"
#CheckCommand "Acrobat Reader" "acroread -help"

if [ -n "$missing" ]; then
  echo ""
  echo "${bold}The following programs seem not be installed: ${normal}"
  echo ""
  echo    "${red}$missing${normal}"
  echo ""
  echo "Please see the messages above to find out if these programs"
  echo "are indispensable to run the model."
  echo "It might be a good idea to install these programs anyway."
  echo "How exactly this is done depends a lot on the type of"
  echo "distribution you are using. You could use a software management"
  echo "tool that comes with your distribution."
  echo "Search for the keyword of the missing file and then install"
  echo "the files on your harddisk."
  echo ""
else
  echo ""
  echo    "${green}All programs installed${normal}"
  echo ""
fi

#CleanUp

# xanim

