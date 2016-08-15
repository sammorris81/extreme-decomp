#!/bin/bash
# CWD=`pwd`
RBIN=$HOME/packages/R/lib64/R/bin/R
PREC=$HOME/repos-git/extreme-decomp/markdown/precipitation
L=$1
CV=$2

echo "Knots: $L, CV: $CV"
$RBIN CMD BATCH --vanilla --no-save $PREC/fit-ebf-$L-$CV.R $PREC/fit-ebf-$L-$CV.out 2>&1

exit 0
