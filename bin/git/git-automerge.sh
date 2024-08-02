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
# tries to automerge some files
#
#--------------------------------------------------------------------

srcdir=fem3d
srcdir=src/utils/shyutil
version=version.f90

git-automerge.pl COMMIT > COMMIT.tmp
git-automerge.pl VERSION > VERSION.tmp
git-automerge.pl $srcdir/$version > $srcdir/$version.tmp

echo "the following files have been merged:"

echo "    COMMIT         -> COMMIT.tmp"
echo "    VERSION        -> VERSION.tmp"
echo "    $srcdir/$version -> $srcdir/$version.tmp"

if [ "$1" = "-write" ]; then
  mv -f COMMIT.tmp COMMIT
  mv -f VERSION.tmp VERSION
  mv -f $srcdir/$version.tmp $srcdir/$version
  echo "files have been written"
elif [ "$1" = "-diff" ]; then
  echo "------------- diffing COMMIT -----------------"
  diff COMMIT.tmp COMMIT
  echo "------------- diffing VERSION ----------------"
  diff VERSION.tmp VERSION
  echo "---------- diffing $srcdir/$version-----------"
  diff $srcdir/$version.tmp $srcdir/$version
  echo "----------------------------------------------"
elif [ "$1" = "-tkdiff" ]; then
  tkdiff COMMIT.tmp COMMIT
  tkdiff VERSION.tmp VERSION
  tkdiff $srcdir/$version.tmp $srcdir/$version
else
  echo "in order to write changes use option -write"
fi

#--------------------------------------------------------------------

