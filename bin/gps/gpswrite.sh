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
# writes extra comments into ps/eps files
#
#--------------------------------------------------------------

script=$(realpath $0)
FEMBIN=$(dirname $script)

gpswrite=$FEMBIN/gpswrite.pl

#--------------------------------------------------------------

NoFile()
{
  echo "No such file: $1"
  Usage
}

Usage()
{
  echo "Usage: gpswrite \"(text) x y points\" file(s)"
  exit 1
}

#--------------------------------------------------------------

string=$1
shift

#--------------------------------------------------------------

if [ $# -eq 0 ]; then
  Usage
fi

echo "String used: $string"

for file
do
  [ -f $file ] || NoFile $file
  echo $file
  $gpswrite "$string AW" $file > tmp.tmp
  mv $file $file.bak
  mv tmp.tmp $file
done

#--------------------------------------------------------------

