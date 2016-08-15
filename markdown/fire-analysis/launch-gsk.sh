#!/bin/bash
# CWD=`pwd`
RBIN=$HOME/packages/R/lib64/R/bin/R
FIRE=$HOME/repos-git/extreme-decomp/markdown/fire-analysis
L=$1
CV=$2

echo "Knots: $L, CV: $CV"
$RBIN CMD BATCH --vanilla --no-save $FIRE/fit-gsk-$L-$CV.R $FIRE/fit-gsk-$L-$CV.out 2>&1

exit 0
