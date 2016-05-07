#!/bin/bash
# CWD=`pwd`

for i in `seq 1 6`
do
  for j in `seq 1 5`
  do
    L=$((i * 5))
    bwsubmit exec launch-ebf.sh $L $j
    bwsubmit exec launch-gsk.sh $L $j
  done
done
exit 0
