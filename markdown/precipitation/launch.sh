#!/bin/bash
# CWD=`pwd`
RBIN=$HOME/packages/R/lib64/R/bin/R
PREC=$HOME/repos-git/extreme-decomp/markdown/precipitation

bwsubmit exec $RBIN CMD BATCH --vanilla --no-save fit-ebf-5-1.R fit-ebf-5-1.out 2>&1

exit 0
