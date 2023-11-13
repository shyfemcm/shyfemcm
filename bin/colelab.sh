#!/bin/sh
#
# elab columns
#
#---------------------------------------------------------------

Usage()
{
  echo "Usage: colelab.sh [-h|-help] [-options] file(s)"
}

FullUsage()
{
  echo ""
  Usage

  echo ""
  echo "Available options:"
  echo "  -h|-help         this help"
  echo "  -verbose         be verbose"
  echo "  -quiet           be quiet"
  echo "  -split           splits file into single columns"
  echo "  -join            joins single files into one file"
  echo "  -extract #       extracts column #"
  echo "  -scale f,c       scales column(s) with factor f and constant c"
  echo "  -col #           operates on column #, otherwise all columns"
  echo "  -name name       use name for created files with -split"
  echo ""
}

ErrorOption()
{
  echo "No such option : $1"
}

#---------------------------------------------------------------

bindir=$( dirname $0 )
#echo "script run from directory $bindir"

what="none"
newname=""
name="col"
col=""
quiet="NO"
verbose="NO"
options=""

while [ -n "$1" ]
do
   case $1 in
        -quiet)         quiet="YES"; options="$options -quiet";;
        -verbose)       verbose="YES"; options="$options -verbose";;
        -split)         what="split";;
        -join)          what="join";;
        -extract)       what="extract"; col=$2; shift;;
        -scale)         what="scale"; scale=$2; shift;;
        -col)           col=$2; shift;;
        -name)          name=$2; shift;;
        -h|-help)       FullUsage; exit 0;;
        -*)             ErrorOption $1; exit 1;;
        *)              break;;
   esac
   shift
done

#---------------------------------------------------------------

[ -n "$name" ] && options="$options -newfile=$name"
[ -n "$col" ] && options="$options -col=$col"
[ -n "$col" ] && options="$options -col=$col"

if [ $# -eq 0 -o $what = "none" ]; then
  Usage
  exit 0
elif [ $what = "split" ]; then
  $bindir/colelab.pl $options -split $1
  [ $quiet = "NO" ] && echo "file split into files $name.n.txt"
elif [ $what = "extract" ]; then
  $bindir/colelab.pl $options -extract=$col $1
  [ $quiet = "NO" ] && echo "extracted column is in extract.$col.txt"
elif [ $what = "scale" ]; then
  $bindir/colelab.pl $options -scale=$scale $1 > scale.txt
  [ $quiet = "NO" ] && echo "scaled file is in scale.txt"
elif [ $what = "join" ]; then
  $bindir/joincol.pl $options $* > join.txt
  [ $quiet = "NO" ] && echo "joined columns are in join.txt"
else
  Usage
  exit 0
fi

#---------------------------------------------------------------

