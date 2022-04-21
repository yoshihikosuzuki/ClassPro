# ClassPro: A K-mer classifier for HiFi reads


A k-mer classification program based on read profiles (a.k.a. k-mer count profiles) computed by the [FASTK](https://github.com/thegenemyers/FASTK) k-mer counter developed by Gene Myers.

The input data should be HiFi reads. The sequencing coverage is roughly assumed to be not too small (e.g. <10x) nor not too high (e.g. >100x). We storngly assume that the ploidy of the underlying genome is diploid.

We also provide a tool for interactive visualization of read profiles (plus k-mer count histogram and k-mer classification result) here: https://github.com/yoshihikosuzuki/kmer-profile.

## How to install 

``` bash
$ git clone https://github.com/yoshihikosuzuki/ClassPro
$ cd ClassPro; make
```

## How to run

```text
ClassPro [-vs] [-c<int>] [-r<int(20000)>] [-N<fastk_root>]
         [-P<tmp_dir(/tmp)>] [-T<int(4)>]
         <source>[.db|.f[ast][aq][.gz]] ...
```

**Mandatory arguments**:

- `<source>`: A single sequence file of HiFi reads with which FASTK was executed. Currently the input file must be exactly a single file (NOTE: multiple input files could be accepted in the future, although the priority is currently low). BAM/CRAM files are currently not supported. FASTK must be executed with `-t1` (all k-mer counts) `-p` (count profile) options. 

**Options**:

- `-N`: Base name of the output files of FASTK executed with `<source>`. If not specified, it is automatically inferred by the file name of `<source>`. There must exist both `<fastk_root>.hist` and `<fastk_root>.prof` (and the hidden files associated with them).
- `-P`: Path to the temporary directory where intermediate files are generated.
- `-v`: Run in verbose mode.
- `-T`: Number of threads used. Too large number (e.g. 100) causes IO bottleneck (with the current implementation, which might be improved in the future).
- `-s`: If specified, ClassPro also finds seeds for alignment based on the classification result. [**Under development**]
- `-c`: Average sequencing coverage (per diploid; = total bases in reads / haploid genome size). If not specified, it is automatically estimated from the k-mer count histogram, which typically works well.
- `-r`: Average read length in bases. A rough estimate is fine.

**Output**:

- A fastq-like file `<fastk_root>.class` where the classification results (which is a sequence of {`E`,`H`,`D`,`R`} preceded by the first (K-1) `N`) are written instead of the quality score values.

---

Other than `ClassPro`, the following programs are compiled for evaluation etc.:

- `ClassGS`: Classification based on the global k-mer count thresholds computed by GenomeScope.
  - `$ genomescope_thresholds.sh <path-to-genomescope-out-dir>` will extract the values of the thresholds (**NOTE**: GenomeScope must be run with the `--fitted_hist` option. An example command is written in the shell script).
- `prof2class`: Given a read dataset and a relative profile of the reads against the underlying genome, generate the classification file based on the relative k-mer counts. If the relative profile is generated with the complete genome from which the reads are sampled (e.g. simulation), then the classification can be used as the ground truth of the k-mer classification.
- `class2acc`: Given a classification file (by ClassPro, ClassGS, etc.) and the ground-truth classification file, calculate the accuracy of the classification.
- `class2cns`: Given a classification file, generate the consensus of classifications for ecah distinct k-mer. [**Under development**]