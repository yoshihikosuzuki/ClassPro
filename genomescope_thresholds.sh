#!/bin/bash

OUT_PREFIX=$1

# 1. Run genomescope.R with the `--fitted_hist` option.
#genomescope.R -i ${OUT_HIST} -o ${OUT_PREFIX} -p ${PLOIDY} -k ${K} --fitted_hist

# 2. Find thresholds from the lookup table
awk -F',' 'prev != $1 {print NR-1 "\t" $0} {prev = $1}' ${OUT_PREFIX}/lookup_table.txt
