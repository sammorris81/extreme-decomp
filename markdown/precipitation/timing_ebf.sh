#!/bin/bash
# CWD=`pwd`
RBIN=$HOME/packages/R/lib64/R/bin/R
PREC=$HOME/repos-git/extreme-decomp/markdown/precipitation

$RBIN CMD BATCH --vanilla --no-save $PREC/get_timing_ebf.R $PREC/get_timing_ebf.out 2>&1

exit 0
