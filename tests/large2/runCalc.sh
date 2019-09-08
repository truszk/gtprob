#!/bin/bash
source /beegfs/home/pulicani/anaconda2/bin/activate benj
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
CALC=/beegfs/home/truszkowsk/calcProbConcordant.py
#CALC=~/speciestrees/calcProbConcordant.py
gtname=$1 #here without nonums
stname=$2
outpname=$1.out
timename=$1.time

echo $outpname
echo $timename
/usr/bin/time -o $timename python $CALC $gtname $stname > $outpname
