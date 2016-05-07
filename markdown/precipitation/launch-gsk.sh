#!/bin/bash
# CWD=`pwd`
RBIN=$HOME/packages/R/lib64/R/bin/R
PREC=$HOME/repos-git/extreme-decomp/markdown/precipitation
L=$1
CV=$2

$RBIN CMD BATCH --vanilla --no-save $PREC/fit-gsk-$L-$CV.R $PREC/fit-gsk-$L-$CV.out 2>&1

exit 0
