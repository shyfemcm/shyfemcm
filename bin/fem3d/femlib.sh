#!/bin/sh
#
#---------------------------------------------

change="NO"
[ "$1" = "-change" ] && change="YES"

#files=$( ls *.pl *.pm )
files=$( findf '*.pl' )

for file in $files
do
  grep femlib $file > /dev/null
  if [ $? -eq 0 ]; then
    if [ $change = "YES" ]; then
      echo "changing  $file"
      ./femlib.pl $file
    else
      echo "found  $file"
    fi
  fi
done

if [ $change = "NO" ]; then
  echo "  nothing changed"
  echo "  use -change to actually change files in place"
fi

