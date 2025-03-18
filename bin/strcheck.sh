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
# checks STR file and its files
#
#-----------------------------------------------------------------

FEMDIR=${SHYFEMDIR:=$HOME/shyfem}
FEMBIN=$FEMDIR/bin
elabdir=$FEMDIR/src/tools/shyelab/

if [ $# -ne 1 ]; then
  echo "Usage: strcheck.sh str-file"
  exit 0
fi

str=$1

#-----------------------------------------------------------------

ElabFile()
{
  first=$( grep 'first time record' tmp2.tmp | sed -e 's/.*: //' )
  last=$( grep 'last time record' tmp2.tmp | sed -e 's/.*:  //' )
  dt=$( grep 'regular time step' tmp2.tmp | sed -e 's/.*:  * //' )

  [ -z "$dt" ] && dt="irreg"

  echo "$first - $last - $dt - $file"
}

ElabSTR()
{
  itanf=$( echo $1 | sed -E "s/'//g" )
  itend=$( echo $2 | sed -E "s/'//g" )
  idt=$( echo $3 | sed -E "s/'//g" )
  
  echo "$itanf - $itend - $idt - STR"
}

#-----------------------------------------------------------------

itanf=`$FEMBIN/strparse.pl -quiet -value=itanf $str`
itend=`$FEMBIN/strparse.pl -quiet -value=itend $str`
date=`$FEMBIN/strparse.pl -quiet -value=date $str`
idt=`$FEMBIN/strparse.pl -quiet -value=idt $str`
[ "$date" = "" ] && date=0

ElabSTR $itanf $itend $idt

$FEMBIN/strparse.pl -quiet -files $str > tmp.tmp

while read line
do
  where=$( echo $line | sed -e 's/ :.*//' )
  what=$( echo $line | sed -e 's/^.* : *//' | sed -e 's/ = .*//' )
  file=$( echo $line | sed -e 's/.*= //' | sed -E "s/'//g" )

  #echo "$where - $what - $file"
  [ $what = "grid" ] && continue
  [ $what = "gotmpa" ] && continue

  $elabdir/shyelab -quiet $file > tmp2.tmp

  ElabFile

  #echo "checking file $file"
  #$FEMBIN/strcheck.pl $itanf $itend $date tmp2.tmp

done < tmp.tmp

rm -f tmp*.tmp

