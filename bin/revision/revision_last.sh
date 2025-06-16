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
# shows revisions after last VERSION change
#
#------------------------------------------------------------

version_file="VERSION"
compare_file=$version_file

bindir=$SHYFEMDIR/bin/revision

tmpfile=tmp0.tmp
tmpfile1=tmp1.tmp
tmpfile2=tmp2.tmp

log="NO"
date=""

#------------------------------------------------------------
# definition of functions
#------------------------------------------------------------

Usage()
{
  echo "Usage: revision_last.sh [ -h | -help ] [ options ]"
}

FullUsage()
{
  echo ""
  Usage

  echo ""
  echo "Available options:"
  echo "  -h|-help         this help"
  echo "  -log             insert output into LOG file"
  echo "  -file file       use file for compare date"
  echo "  -date date       use date for compare date (YYYY-MM-DD)"
  echo ""
}

ErrorOption()
{
  echo "No such option : $1"
}

ErrorExtra()
{
  echo "Extra information on command line : $1"
}

#------------------------------------------------------------

date2num()
{
  echo $1 | sed 's/\(.*\)-\(.*\)-\(.*\)/\3\2\1/'
}

GetDate()
{
  if [ -n "$date" ]; then
    comparedate=`echo $date | sed -e 's/-//g'`
  elif [ $compare_file  = "VERSION" ]; then
    firstline=`head -1 $version_file | sed 's/ \{1,\}/ /g'`
    actdate=`echo $firstline | cut -d" " -f3`
    comparedate=`date2num $actdate`
  else
    actdate=`ls -l --time-style=full-iso $compare_file | cut -f 6 -d " "`
    comparedate=`echo $actdate | sed -e 's/-//g'`
  fi
}

#------------------------------------------------------------
# read in command line options
#------------------------------------------------------------

while [ -n "$1" ]
do  
   case $1 in
        -log)           log="YES";;
        -file)          compare_file=$2; shift;;
        -date)          date=$2; shift;;
        -h|-help)       FullUsage; exit 0;;
        -*)             ErrorOption $1; exit 1;;
        *)              ErrorExtra $1; exit 1;;
   esac
   shift
done

if [ ! -f $version_file ]; then
  echo "No file $version_file"
  echo "You must run this command in the main shyfem directory."
  exit 1
fi

#------------------------------------------------------------
# get compare date
#------------------------------------------------------------

GetDate
echo "using compare date $comparedate"

#------------------------------------------------------------
# find files
#------------------------------------------------------------

files=`find . -name "*.f90"`
#files=`find . -name "*.[cfFh]"`
#echo "looking for files: $files"

$bindir/revisionlog.sh -after $comparedate $files > $tmpfile
$bindir/revisionlog_adjust.pl $tmpfile $version_file > $tmpfile2

cat $tmpfile2

#------------------------------------------------------------
# insert into LOG file
#------------------------------------------------------------

if [ "$log" = "YES" ]; then	#add output to LOG
  echo "adding revision log to file LOG"
  date=`date +"%d.%m.%Y"`
  cp LOG LOG.tmp
  echo "" 		 > tmp.tmp
  echo "$date"		>> tmp.tmp
  echo "" 		>> tmp.tmp
  cat $tmpfile2		>> tmp.tmp
  echo "" 		>> tmp.tmp
  cat LOG		>> tmp.tmp
  mv tmp.tmp LOG
fi

#------------------------------------------------------------
# clean up
#------------------------------------------------------------

#rm -f $tmpfile $tmpfile1 $tmpfile2

#------------------------------------------------------------
# end of routine into LOG file
#------------------------------------------------------------

