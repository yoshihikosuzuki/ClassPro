import argparse
from collections import Counter
from scipy.stats import hmean
import kmer_profiler as kp


def calc_mcf(counter):
    return counter.most_common()[0][1] / sum(counter.values())


def calc_consistency(class_table):
    mcfs = []
    counter = Counter()
    prev_kmer = ""
    i = 0
    with open(class_table, 'r') as f:
        for line in f:
            cnt, kmer, c = line.strip().split()
            cnt = int(cnt)
            if prev_kmer != "" and prev_kmer != kmer:
                mcf = calc_mcf(counter)
                mcfs.append(mcf)
                counter = Counter()
            counter[c] = cnt
            prev_kmer = kmer
            i += 1
            if i % 100000 == 0:
                print(f"{i} (Consistency so far = {hmean(mcfs)})")
        mcf = calc_mcf(counter)
        mcfs.append(mcf)
    print(f"Overall consistency = {hmean(mcfs)}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute consistency from classification tables per k-mer")
    parser.add_argument(
        "kmer_cnts",
        type=str,
        help="File name of k-mer classification table.")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    calc_consistency(args.kmer_cnts)
