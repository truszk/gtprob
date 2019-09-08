#!/bin/bash
STELLS=/beegfs/home/truszkowsk/STELLSv2.1.0/stells
gtname=$1
stname=$2
outpname=$1.out
timename=$1.time

echo $outpname
echo $timename
/usr/bin/time -o $timename $STELLS -g $gtname -s $stname -B -O > $outpname
