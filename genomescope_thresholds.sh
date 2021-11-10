#!/bin/bash
# Find the count thresholds between Haplo/Diplo/Repeat based on the GenomeScope's fitting to the count histogram.
# Usage: ./genomescope_thresholds.sh <genomescope-out-dir>
# [NOTE] GenomeScope must be run with the `--fitted_hist` option; e.g. `$ genomescope.R -i ${OUT_HIST} -o ${OUT_PREFIX} -p ${PLOIDY} -k ${K} --fitted_hist`

OUT_PREFIX=$1

awk -F',' 'prev != $1 {print NR-1 "\t" $0} {prev = $1}' ${OUT_PREFIX}/lookup_table.txt
