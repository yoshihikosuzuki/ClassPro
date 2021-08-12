# ClassPro: A K-mer classifier for HiFi reads

A k-mer classification program using the [FASTK](https://github.com/thegenemyers/FASTK) software suite developed by Gene Myers.

```text
ClassPro [-vs] [-c<int>] [-r<int(20000)>] [-N<fastk_root>]
         [-P<tmp_dir(/tmp)>] [-T<int(4)>]
         <source>[.db|.f[ast][aq][.gz]] ...
```

Mandatory arguments:

- `source`: Sequence file(s) with which FASTK was executed. When multiple files are give, the order of those files must be same as that for FASTK (This is checked during the classification) and all files must have the same extension. Only a single file is accepted for `.db`.

Options:

- `-N`: Base name of the FASTK's output files given the input sequence files. If not specified, it is automatically inferred by the first sequence file name. Both `<fastk_root>.hist` and `<fastk_root>.prof` (and the hidden files associated with them) must exist.
- `-P`: Path to a temporary directory where intermediate files are generated.
- `-v`: Run in verbose mode.
- `-T`: Number of threads used.
- `-s`: If specified, ClassPro also finds seeds for alignment based on the classification result.
- `-c`: Average sequencing coverage (per diploid; = total bases in reads / haploid genome size). If not specified, it is automatically estimated from the k-mer count histogram, which works well in most cases.
- `-r`: Average read length in bases. A rough estimate is fine.

Output:

- A fastq-like file `<fastk_root>.class` where
  - Sequences of {E,H,D,R} instead of quality score values indicates classification results; and
  - Nucleotides in capital case indicate alignment seeds (if `-s` is specified).
- If the input sequence file is `.db`, then the followings are additionally generated:
  - A DAZZ_DB track `.<source>.class.[anno|data]` for classification results (0: Error, 1: Repeat, 2: Haplo, 3: Diplo); and
  - A DAZZ_DB track `.<source>.seed.[anno|data]` for alignment seeds (if `-s` is specified).
