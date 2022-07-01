# ClassPro: A K-mer classifier for HiFi reads

A k-mer classification program based on read profiles (a.k.a. k-mer count profiles) computed by the [FASTK](https://github.com/thegenemyers/FASTK) k-mer counter developed by Gene Myers.

- We storngly assume that the ploidy of the underlying genome is **diploid** (Extension to polyploid genomes is a future work).
- The input data should be HiFi reads (for now).
- The global (diploid) sequencing coverage is roughly assumed to be not too small (e.g. <10x) nor not too high (e.g. >100x).

We also provide a tool for interactive visualization of read profiles (plus k-mer count histogram and k-mer classification result) [here](https://github.com/yoshihikosuzuki/kmer-profile), which will be integrated to this repository in the near future.

## How to install

### From a release tarball

```bash
wget XXX
cd ClassPro-XXX
make
```

### From the latest source

``` bash
git clone https://github.com/yoshihikosuzuki/ClassPro
cd ClassPro
make
```

## Workflow example

In addition to ClassPro, you need to install FASTK (and also DAZZ_DB if you want to use `.db`/`.dam` files as input).

We provide `mhc_reads.fasta` a simulation 40x HiFi reads of the human MHC region generated with [HIsim]() using a currently available complete diploid assembly (Chin et al, Nat. Comm., 2020).

```bash
FastK -v -k40 -t1 -p mhc_reads.fasta
ClassPro -v mhc_reads
```

This is a test dataset with the ground-truth classification, and we can calculate the accuracy:

```bash
FastK -v -k40 -t1 -p mhc_haps.fasta
FastK -v -k40 -p:mhc_haps -Nmhc_reads.truth mhc_reads.fasta
prof2class mhc_reads.truth XXX
class2acc mhc_reads.class mhc_reads.truth.class
```

GenomeScope-based classification:

```bash
GenomeScope
calc_thres
ClassGS X Y Z mhc_reads.fasta
class2acc mhc_reads.GS.class mhc_reads.truth.class
```

## How to run

```text
ClassPro [-vs] [-c<int>] [-r<int(20000)>] [-N<fastk_root>]
         [-P<tmp_dir(/tmp)>] [-T<int(4)>]
         <source>[.db|.f[ast][aq][.gz]]
```

**Mandatory arguments**:

- `<source>`: A single sequence file of HiFi reads with which FASTK was executed. Currently the input file must be exactly a single file (NOTE: multiple input files could be accepted in the future, although the priority is currently low). BAM/CRAM files are currently not supported. FASTK must be executed with `-t1` (all k-mer counts) `-p` (count profile) options.

**Options**:

- `-v`: Run in verbose mode.
- `-T`: Number of threads used. Typically there will be almost no gain with a number larger than 16 given a single SSD because of IO bottleneck (with the current implementation, which might be improved in the future).
- `-P`: Path to the temporary directory where intermediate files are generated.
- `-N`: Base name of the output files of FASTK executed with `<source>`. If not specified, it is automatically inferred by the file name of `<source>`. There must exist both `<fastk_root>.hist` and `<fastk_root>.prof` (and the hidden files associated with them).
- `-c`: Average sequencing coverage (per diploid; = total bases in reads / haploid genome size). If not specified, it is automatically estimated from the k-mer count histogram, which typically works well.
- `-r`: Average read length in bases. A rough estimate is fine.
- `-s`: If specified, ClassPro also finds seeds for alignment based on the classification result. [Preliminary feature]

**Output**:

- A fastq-like file `<fastk_root>.class` where the classification results (which is a sequence of {`E`,`H`,`D`,`R`} preceded by the first (K-1) `N`) are written instead of the quality score values.

---

Other than `ClassPro`, the following programs are compiled for evaluation etc.:

- `ClassGS`: Classification based on the global k-mer count thresholds computed by GenomeScope.
  - `$ genomescope_thresholds.sh <path-to-genomescope-out-dir>` will extract the values of the thresholds (**NOTE**: GenomeScope must be run with the `--fitted_hist` option. An example command is written in the shell script).
- `prof2class`: Given a read dataset and a relative profile of the reads against the underlying genome, generate the classification file based on the relative k-mer counts. If the relative profile is generated with the complete genome from which the reads are sampled (e.g. simulation), then the classification can be used as the ground truth of the k-mer classification.
- `class2acc`: Given a classification file (by ClassPro, ClassGS, etc.) and the ground-truth classification file, calculate the accuracy of the classification.
- `class2cns`: Given a classification file, generate the consensus of classifications for ecah distinct k-mer. [**Under development; Currently very inefficient implementation**]

## Limitations

- ClassPro cannot be accurate for noisy reads (at least in the current implementation) because it leverages the coherence principle of read profiles for k-mer classificaion. High noise violates this principle.
- Low-copy (>2) repeat-mers are one of the most dificult types for ClassPro. For example, if all the k-mers except error-mers in a single read profile are three-copy repeat-mers with an average count of 60 given the average diploid sequencing coverage is 40, then it is very possible that those k-mers are diplo-mers rather than repeat-mers, and vice versa. Thus, the classification accuracy tends to be low around such regions. This is not so much harmful for read alignment and genome assembly, but is for direct identification of haplotype-specific k-mers from the k-mer classification.

## TODOs

- For an initial separation between diplo-mers and repeat-mers, ClassPro uses a loose threshold $\theta_{\text R}$ between diplo-mers and repeat-mers, where $\theta_{\text R} := c_{\text D}$ + `N_SIGMA_RCOV` * $\sqrt{c_{\text D}}$ and `N_SIGMA_RCOV` is a hyperparameter defined in `const.c`, and ClassPro does not find a wall inside repeat regions. When k-mer counts are gradually increasing from $<\theta_{\text R}$ to $>\theta_{\text R}$, the current implementation does not recognize intervals around the boundary well.
- Integrate the profile visualizer into this repository.
