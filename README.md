# ClassPro

A k-mer classification program using the [FASTK](https://github.com/thegenemyers/FASTK) library by Gene Myers.

```bash
ClassPro [-vs] [-T<int(4)>] [-c<int>] [-r<int(20000)>] <seq_file> <fastk_root>
```

Mandatory arguments:

- `seq_file`: A `.fasta`, `.fastq`, or `.db` file (TODO: fofn?)
- `fastk_root`: Base name of the output files of FASTK for `seq_file`. Both `<fastk_root>.hist` and `<fastk_root>.prof` must exist.

Options:

- `-v`: Run in verbose mode
- `-s`: Find seeds for alignment
- `-T`: Number of threads
- `-c`: Average sequencing coverage (per diploid; = total bases in reads / haploid genome size)
  - (if available; in most cases it is automatically estimated from the k-mer count histogram)
- `-r`: Average read length (in bases)
  - (rough guess is ok)

Output:

- If `seq_file` is `.fasta` or `.fastq`:
  - A fastq-like file `<source_root>.class` where
    - Sequences of {E,H,D,R} instead of quality score values indicates classification results; and
    - Nucleotides in capital case indicate alignment seeds
- Otherwise (i.e. `seq_file` is `.db`):
  - A DAZZ_DB track `.<seq_file>.class.[anno|data]` for classification results (0: Error, 1: Repeat, 2: Haplo, 3: Diplo); and
  - A DAZZ_DB track `.<seq_file>.seed.[anno|data]` for alignment seeds
