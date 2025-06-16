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
# handle postscript files
#
# needed programs:
#
# parr (with option -A)
# gps2eps
# psresetbb.sh
# gifcrop (needs netpbm routines)
#
# epstopdf	(/usr/bin/epstopdf)
#
# convert	(Imagemagick)
# gifsicle	(gifsicle)
#
# 02.05.2013    ggu     allow for single pages to be extracted
# 11.01.2025    ggu     some minor changes for usability
#
######################################################### copyright

copy2="gps - shell for generation of plot files"
copy3="Copyright (c) 1998-2025 Georg Umgiesser - ISMAR-CNR"
copy4="e-mail: georg.umgiesser@cnr.it"

######################################################### help - usage

script=$(realpath $0)
scriptdir=$(dirname $script)

gps2eps=$scriptdir/gps2eps.pl
psresetbb=$scriptdir/psresetbb.sh
parr=$scriptdir/parr.pl
ps2pdf=/usr/bin/ps2pdf

#--------------------------------------------------------

CheckExec()
{
  for r in $*
  do
    if [ ! -x $r ]; then
      echo "No such executable: $r"
      exit 1
    fi
  done
}

#--------------------------------------------------------

Copyright()
{
  echo "$copy2"
  echo "$copy3"
}

Usage()
{
  echo "Usage: gps.sh [ -h | -help ] [ -options ] file[s]"
}

FullUsage()
{
  echo ""

  Copyright
  echo ""
  Usage

  echo ""
  echo "Available options:"
  echo "  -h|-help         this help"
  echo "  -quiet           be quiet"
  echo "  -verbose         be verbose"
  echo "  -debug           write some debug messages"
  echo "  -split           split ps file into single pages"
  echo "  -pages range     extract single pages (range: 1,3,4-8)"
  echo "  -eps             convert ps files into eps files"
  echo "  -pdf             convert (e)ps files into pdf files"
  echo "  -gif             convert files into gif files"
  echo "  -jpg             convert files into jpg files"
  echo "  -png             convert files into png files"
  echo "  -bmp             convert files into bmp files"
  echo "  -dpi #           for raster files use this resolution in dpi"
  echo ""
}

Debug()
{
  echo "SHYFEMDIR: $SHYFEMDIR"
  echo "scriptname: $script"
  echo "scriptdir: $scriptdir"
}

ErrorOption()
{
  echo "No such option : $1"
}

Extension()
{
  file=$1
  ext=`echo $file | sed -e 's/^.*\././'`
  echo $ext
}

Basename()
{
  file=$1
  ext=`Extension $file`
  base=`basename $file $ext`
  echo $base
}

CleanExtra()
{
  # cleans up files like "plot.png.1" etc...

  name=$1
  if [ -f $name.0 ]; then
    mv $name.0 $name
    rm -f $name.[1-9]*
  fi
}

MustRotate()
{
  name=$1
  ext=$2

  if [ $ext = ".eps" -o $ext = ".ps" ]; then
    grep "^%%Orientation: Landscape" $name > /dev/null
    status=$?
    #echo "eps or ps: $name   -  $status"
    if [ $status -eq 0 ]; then
      echo "YES"
    else
      echo "NO"
    fi
  else
    echo "NO"
  fi
}

RotateEps()
{
  name=$1
  ext=$2

  if [ $ext = ".eps" -o $ext = ".ps" ]; then
    grep "^%%Orientation: Landscape" $name > /dev/null
    status=$?
    #echo "eps or ps: $name   -  $status"
    if [ $status -eq 0 ]; then
      echo "-rotate 90"
    else
      echo ""
    fi
  else
    echo ""
  fi
}

ReduceGif()
{
  # this has to be done since convert creates too big GIF files
  # this also eliminates the second (useless) frame from a gif image

  newext=".gif"

  for file in $*
  do
    base=`Basename $file`
    target=$base$newext

    [ $quiet = "NO" ] && echo "reducing  $target"

    #gifcrop -o $target > tmp.gif
    #[ $? -eq 0 ] && mv tmp.gif $target
    #gifsicle $target '#0' > tmp.gif
    #[ $? -eq 0 ] && mv tmp.gif $target
  done
}

ConvertJpg()
{
  [ -n "$dpi" ] && dpi_opt="-r$dpi"

  for file in $*
  do
    base=`Basename $file`
    ext=`Extension $file`
    source=$base$ext
    target=$base.jpg

    [ $quiet = "NO" ] && echo "$source  ->  $target"
    #echo "options: $dpi $rotate $format"

#   rotate="-c '<</Orientation 3>> setpagedevice'"	#these are not working
#   rotate='-c \"<</Orientation 3>> setpagedevice\"'
#   rotate='-c "<</Orientation 3>> setpagedevice"'

    rotate=`RotateEps $source $ext`	# this rotates Landscape ps/eps files

    gs -sDEVICE=jpeg -dEPSCrop -dNOPAUSE -dBATCH -dSAFER  $dpi_opt \
	-sOutputFile=$target \
	-f $source > /dev/null

    if [ -n "$rotate" ]; then
      mogrify $rotate -trim -border 20 -bordercolor white $target
    fi
  done
}

Convert()
{
  # in $format is the file type if needed (only for jpeg...)

  newext=$1
  shift

  [ -n "$dpi" ] && dpi_opt="-units PixelsPerInch -density $dpi"

  for file in $*
  do
    base=`Basename $file`
    ext=`Extension $file`
    source=$base$ext
    target=$base$newext

    [ $quiet = "NO" ] && echo "$source  ->  $target"
    #echo "options: $dpi $rotate $format"

    rotate=`RotateEps $source $ext`	# this rotates Landscape ps/eps files

    CleanExtra $target
    convert -trim $dpi_opt $rotate $source $format$target
    CleanExtra $target
  done
}

ConvertPDF()
{
  newext=$1
  shift

  for file in $*
  do
    base=`Basename $file`
    ext=`Extension $file`
    source=$base$ext
    target=$base$newext

    [ $quiet = "NO" ] && echo "$source  ->  $target"

    $ps2pdf $source
  done
}

######################################################### end

#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################

######################################################### defaults

quiet="NO"
verbose="NO"
debug="NO"
split="NO"
pages=""
eps="NO"
pdf="NO"
gif="NO"
jpg="NO"
png="NO"
bmp="NO"
dpi=""

######################################################### read options

while [ -n "$1" ]
do
   case $1 in
	-quiet)		quiet="YES";;
	-verbose)	verbose="YES";;
	-debug)		debug="YES";;
	-split)		split="YES";;
	-pages)		pages=$2; shift;;
	-eps)		eps="YES";;
	-pdf)		pdf="YES";;
	-gif)		gif="YES";;
	-jpg)		jpg="YES";;
	-png)		png="YES";;
	-bmp)		bmp="YES";;
	-dpi)		dpi=$2; shift;;
	-h|-help)	FullUsage; exit 0;;
	-*)		ErrorOption $1; exit 1;;
	*)		break;;
   esac
   shift
done

######################################################### no file -> write help

CheckExec $gps2eps $psresetbb $parr $ps2pdf

[ $debug = "YES" ] && Debug

if [ $# -le 0 ]; then
  Usage
  exit 1;
fi

######################################################### split files

newfiles=""

if [ $split = "YES" -o -n "$pages" ]; then
  for file in $*
  do
    filename=`basename $file .ps`
    rm -f $filename.[0-9]*.ps
    if [ $split = "YES" ]; then
      $parr -A $file
    else
      $parr -A -o $pages $file
    fi
    newfile=`ls -m $filename.[0-9]*.ps`
    [ $quiet = "NO" ] && echo "$file -> $newfile"
    newfiles="$newfiles $newfile"
  done
else
  newfiles=`ls -m $*`
fi

newfiles=`echo $newfiles | sed -e 's/,/ /g' | sed -e 's/  */ /g'`

[ $verbose = "YES" ] && echo "$newfiles"

######################################################### look for 0 size files

auxfiles=""
for file in $newfiles
do
  if [ -s $file ]; then
    auxfiles="$auxfiles $file"
  else
    echo "file $file has size 0 ... ignoring"
  fi
done
newfiles=$auxfiles

######################################################### ps -> eps

if [ $eps = "YES" ]; then
  auxfiles=""
  for file in $newfiles
  do
    basefile=`Basename $file`
    ext=`Extension $file`
    psfile=$basefile$ext
    epsfile=$basefile.eps

    [ $quiet = "NO" ] && echo "$psfile  ->  $epsfile"

    #ps2eps $psfile -o$epsfile  2> /dev/null
    $gps2eps $psfile > $epsfile
    $psresetbb $epsfile
    auxfiles="$auxfiles $epsfile"
  done
  newfiles=$auxfiles
fi

# $newfiles 	->	files to convert

######################################################### 

if [ $gif = "YES" ]; then
  Convert ".gif" $newfiles
  ReduceGif $newfiles
fi

if [ $jpg = "YES" ]; then
  format="jpeg:"
  #Convert ".jpg" $newfiles
  format=""
  ConvertJpg $newfiles
fi

if [ $png = "YES" ]; then
  Convert ".png" $newfiles
fi

if [ $bmp = "YES" ]; then
  Convert ".bmp" $newfiles
fi

if [ $pdf = "YES" ]; then
  ConvertPDF ".pdf" $newfiles
fi

######################################################### end

