#!/bin/sh

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

# make scripts executable

dirs="femanim femdoc grid femcheck femplot fem3d fem3d/bin \
	femersem femspline"

CHX()
{
  dir=$1
  echo "making scripts executable in directory $dir"
  shift
  chmod +x $* 2> /dev/null
  #chmod +x $*
}

for dir in $dirs
do
  CHX $dir $dir/*.sh $dir/*.pl
done

dir=lib/perl
CHX $dir $dir/*.sh $dir/*.pl $dir/*.pm

dir=lib/python
CHX $dir $dir/*.sh $dir/*.py

exit 0

