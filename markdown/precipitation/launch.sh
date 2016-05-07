#!/bin/bash
# CWD=`pwd`

for i in `seq 1 2`
do
  for j in `seq 1 2`
  do
    L=$((i * 5))
    bwsubmit exec launch-ebf.sh $L $j
  done
done
exit 0
