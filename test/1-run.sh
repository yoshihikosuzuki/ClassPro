#!/bin/bash
set -e

IN_READS=mhc_reads.fasta.gz
IN_GENOME=mhc_genome.fasta.gz

READS_PREF=${IN_READS%.fasta.gz}
GENOME_PREF=${IN_GENOME%.fasta.gz}

# 1. Run ClassPro
FastK -v -k40 -t1 -p ${IN_READS}
ClassPro -v ${IN_READS}

# 2-1. Compute ground-truth classification
FastK -v -k40 -t1 -p ${IN_GENOME}
FastK -v -k40 -T1 -p:${GENOME_PREF} -N${READS_PREF}.truth ${IN_READS}
prof2class ${READS_PREF}.truth.prof ${IN_READS}

# 2-2. ClassPro's accuracy
class2acc ${READS_PREF}.class ${READS_PREF}.truth.class

# 3. GenomeScope classification
kmc -k40 -m30 -ci1 -cs1000 -t1 -fm ${IN_READS} ${READS_PREF}.gs .
kmc_tools transform ${READS_PREF}.gs -cx1000 histogram ${READS_PREF}.gs.hist
genomescope.R -i ${READS_PREF}.gs.hist -o ${READS_PREF}.gs -p 2 -k 40 --fitted_hist -l 20
genomescope_thresholds.sh ${READS_PREF}.gs >${READS_PREF}.thres
EH_THRES=$(cut -d' ' -f1 ${READS_PREF}.thres)
let EH_THRES++
HD_THRES=$(cut -d' ' -f2 ${READS_PREF}.thres)
let HD_THRES++
DR_THRES=$(cut -d' ' -f3 ${READS_PREF}.thres)
let DR_THRES++
gunzip -c ${IN_READS} >${READS_PREF}.fasta
fasta2DAM ${READS_PREF} ${READS_PREF}.fasta
DBsplit -s500 ${READS_PREF}
ClassGS ${READS_PREF} ${EH_THRES} ${HD_THRES} ${DR_THRES}

# 4.GenomeScope's accuracy
class2acc ${READS_PREF}.GS.class ${READS_PREF}.truth.class
