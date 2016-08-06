#!/bin/bash
# CWD=`pwd`
RBIN=$HOME/packages/R/lib64/R/bin/R
FIRE=$HOME/repos-git/extreme-decomp/markdown/fire-analysis

$RBIN CMD BATCH --vanilla --no-save $FIRE/get_timing_gsk.R $FIRE/get_timing_gsk.out 2>&1

exit 0
