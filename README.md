# ClassPro: A K-mer classifier (currently for HiFi reads)

ClassPro is a k-mer classification program based on read profiles (a.k.a. k-mer count profiles) computed by the [FASTK](https://github.com/thegenemyers/FASTK) k-mer counter developed by Gene Myers.
ClassPro classifies every k-mer in each read into one of the four types: **Erroneous, Haploid, Diploid, and Repeat**.

- We storngly assume that the ploidy of the underlying genome is **diploid**.
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

After generating binary executables, in some way you need to make the executables and shell scripts seen via `$PATH`. One choice is to use `$ make install` after edting `INSTALL_DIR` in `Makefile`.

In addition to ClassPro, you need to install [FASTK](https://github.com/thegenemyers/FASTK) to generate input files for ClassPro (i.e. count histogram and count profiles).
You also need to install [DAZZ_DB](https://github.com/thegenemyers/DAZZ_DB) if you want to use `.db`/`.dam` files as input (instead of `.fastx[.gz]` files).

## Example work flow

In the `test/` directory, we put an example bash script to run ClassPro (and other related commands) for a simulation dataset.

First move to `test/` and download the data files by:

```bash
. 0-download.sh
```

and then `mhc_genome.fasta.gz` (a currently available complete diploid assembly of the human MHC region by [Chin et al](https://www.nature.com/articles/s41467-020-18564-9)) and `mhc_reads.fasta` (a simulation 40x HiFi reads generated from the assembly) are downloaded.

The other script, `1-run.sh`, contains commands to 1) run ClassPro, 2) evaluate ClassPro's classification using the ground-truth haplotype sequences in the assembly, 3) run GenomeScope, and 4) evaluate GenomeScope's classification.

For example, the commands to run ClassPro are as follows:

```bash
FastK -v -k40 -t1 -p mhc_reads.fasta.gz
ClassPro -v mhc_reads
```

which generate the output file named `mhc_reads.class`. This is a fastq-like file where each character in the quality field represents the classification of the k-mer ending at the position and is one of `E`(rror), `H`(aploid), `D`(diploid), `R`(epeat), or `N` (for the first k-1 bases where k-mers are undefined):

```
@Sim 1 1 + 42 18909
GGACAACAGGCTTGCGCCATCACTGGAGCTGTTCTTAAATTTTTTTAGAGTT...
+
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNHHHHHHHHHHHHH...
...
```

Note that you need to install [GenomeScope2.0](https://github.com/tbenavi1/genomescope2.0) and [DAZZ_DB](https://github.com/thegenemyers/DAZZ_DB) to run GenomeScope part of the script.

## `ClassPro`: The main command

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

## Other programs and scripts

Other than `ClassPro`, the following programs are compiled for evaluation etc.:

- `ClassGS`: Classification based on the global k-mer count thresholds computed by GenomeScope.
  - `$ genomescope_thresholds.sh <path-to-genomescope-out-dir>` will extract the values of the thresholds (**NOTE**: GenomeScope must be run with the `--fitted_hist` option. An example command is written in the shell script).
- `prof2class`: Given a read dataset and a relative profile of the reads against the underlying genome, generate the classification file based on the relative k-mer counts. If the relative profile is generated with the complete genome from which the reads are sampled (e.g. simulation), then the classification can be used as the ground truth of the k-mer classification.
- `class2acc`: Given a classification file (by ClassPro, ClassGS, etc.) and the ground-truth classification file, calculate the accuracy of the classification.
- `class2cns`: Given a classification file, generate the consensus of classifications for ecah distinct k-mer. [**Under development; Currently very inefficient implementation**]

## Limitations

- In the current implementation ClassPro assumes only HiFi reads as input and the same performance is not guaranteed for noisy reads such as ONT reads. We are currently working on this.
- We do not accept polyploid genomes. Generalization to polyploid genomes is one future direction, although it will take some time.

## TODOs

- For an initial separation between diplo-mers and repeat-mers, ClassPro uses a loose threshold $\theta_{\text R}$ between diplo-mers and repeat-mers, where $\theta_{\text R} := c_{\text D}$ + `N_SIGMA_RCOV` * $\sqrt{c_{\text D}}$ and `N_SIGMA_RCOV` is a hyperparameter defined in `const.c`, and ClassPro does not find a wall inside repeat regions. When k-mer counts are gradually increasing from $<\theta_{\text R}$ to $>\theta_{\text R}$, the current implementation does not recognize intervals around the boundary well. Fix this.
- Integrate the profile visualizer into this repository.
- Accept fastx input in ClassGS
- Efficient generation of consensus k-mer classifications.
