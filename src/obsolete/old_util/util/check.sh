#!/bin/sh
#
# checks if all files are available
#
#-----------------------------------------

copy=NO
[ "$1" = "-copy" ] && copy=YES

dir=..
new=utils/generic

files=$( ls *.f90 )

for file in $files
do
  if [ -f $dir/$file ]; then
    cmp $file $dir/$file > /dev/null
    status=$?
    if [ $status -eq 0 ]; then
      :
      #echo "file exists and is identical: $file"
    else
      echo "*** file exists but is different: $file"
    fi
  else
    echo "*** file not existing: $file"
  fi
done

if [ $copy = NO ]; then
  echo "not copying to $dir/$new ... use -copy to copy"
  exit 1
else
  echo "copying to $new"
fi

cd $dir
echo "new: $dir/$new"
ls $files
#cp $files $new
git mv $files $new

#-----------------------------------------



