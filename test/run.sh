#!/bin/bash



FastK -v -t1 -p -T4 hifi.fastq
ClassPro -v -T4 hifi.fastq

FastK -v -t1 -p -T4 ref.fasta
prof2class
class2acc
