#!/bin/sh
#
# finds file or pattern in source directories
#
#-----------------------------------------------

base="$HOME/shyfemcm/shyfemcm"

src="src src/util src/shyutil src/shmpi src/shyutilmpi"
tools="shybas shyelab shyplot shypre"

alldir="$src $tools"

actdir=$( pwd )

#-----------------------------------------------

do_file=YES
do_pattern=NO
if [ $# -gt 0 -a $1 = "-file" ] && do_file=YES && shift
if [ $# -gt 0 -a $1 = "-pattern" ] && do_pattern=YES && shift

if [ $do_file = "YES" ]; then
  for dir in alldir
  do
    cd $base/$dir
    ls $*
    cd $actdir
  done
elif [ $do_pattern = "YES" ]; then
  :
else
  echo "do not know what to do..."
fi

#-----------------------------------------------

