#!/bin/sh
#
#------------------------------------------------------------

if [ $# -eq 0 ]; then
  echo "Usage: wait_to_finish.sh program [secs]"
  exit 1
fi

defsecs=10

program=$1
secs=$2
[ -z "$secs" ] && secs=$defsecs

#------------------------------------------------------------

while [ : ]
do
  nprogs=$( ps aux | grep $program | grep -v grep | grep -v wait_to_finish | wc -l )
  [ $nprogs -eq 0 ] && break
  echo "still $nprogs processes running... sleeping for $secs seconds"
  sleep $secs
done

echo "process finished"

#------------------------------------------------------------

