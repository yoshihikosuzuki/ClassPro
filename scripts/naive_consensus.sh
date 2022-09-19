#!/bin/bash
# Compute the consistency among the classifications for the same k-mer.
# Usage: ./naive_consensus.sh <classification>.class <fastk-root> <num-threads>

IN_CLASS=$1
IN_FASTK_PREFIX=$2
N_THREADS=$3

OUT_PREFIX=${IN_CLASS%.*}
OUT_KMER=${OUT_PREFIX}.kmers
OUT_SORT_KMER=${OUT_PREFIX}.sorted.kmers
OUT_KMER_CNT=${OUT_SORT_KMER}.cnt

~/develop/ClassPro/class2cns ${IN_CLASS} ${IN_FASTK_PREFIX} >${OUT_KMER}
sort --parallel ${N_THREADS} ${OUT_KMER} >${OUT_SORT_KMER}
uniq -c ${OUT_SORT_KMER} >${OUT_KMER_CNT}
python3 ~/develop/ClassPro/agg2cons.py ${OUT_KMER_CNT}
