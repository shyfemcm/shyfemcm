#!/bin/bash

bindir=~/shyfemcm/bin
compdir=~/shyfem/fem3d
dplus="$bindir/d+.pl"

files=$( ls mpi_*.f90 )

for file in $files
do
  name=$( basename $file .f90 )
  skelname=$( echo $name | sed -e 's/mpi_//' )
  #echo "$file $name $skelname"

  if [ -f $compdir/submpi_$skelname.f ]; then
    compfile=$compdir/submpi_$skelname.f
  elif [ -f $compdir/shympi_$skelname.f ]; then
    compfile=$compdir/shympi_$skelname.f
  elif [ -f $compdir/mod_shympi_$skelname.f ]; then
    compfile=$compdir/mod_shympi_$skelname.f
  else
    echo "no such file for $skelname"
    exit 1
  fi

  writefile=$( basename $compfile )
  #echo "$file with $compfile $writefile"

  $bindir/fem3d/f902f.pl -old $file > $writefile
  
  difference=$( diff -wn $compfile $writefile | $dplus )
  echo $writefile: $difference
done

