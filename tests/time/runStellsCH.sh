#!/bin/bash
STELLS=/beegfs/home/truszkowsk/STELLSv2.1.0/stells
#STELLS= /home/jakub/stells/STELLS2-master/stells-v2-1-0-linux64
gtname=$1
stname=$2
outpname=$1.outch
timename=$1.timech

echo $outpname
echo $timename
/usr/bin/time -o $timename $STELLS -g $gtname -s $stname -B -f > $outpname
