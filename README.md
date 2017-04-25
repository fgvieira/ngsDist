# ngsDist

`ngsDist` is a program to estimate pairwise genetic distances directly, taking the uncertainty of genotype's assignation into account. It does so by avoiding genotype calling and using genotype likelihoods or posterior probabilities.

### Citation

`ngsDist` was published in 2015 at [Biological Journal of the Linnean Society](http://onlinelibrary.wiley.com/doi/10.1111/bij.12511/abstract), so please cite it if you use it in your work:

    Vieira FG, Lassalle F, Korneliussen TS, Fumagalli M
    Improving the estimation of genetic distances from Next-Generation Sequencing data
    Biological Journal of the Linnean Society (2015) doi: 10.1111/bij.12511

### Installation

`ngsDist` can be easily installed but has some external dependencies:

* `gcc`: >= 4.9.25 tested on Debian 7.8 (wheezy)
* `zlib`: v1.2.7 tested on Debian 7.8 (wheezy)
* `gsl` : v1.15 tested on Debian 7.8 (wheezy)
* `md5sum`: only needed for `make test`

To install the entire package just download the source code:

    % git clone https://github.com/fgvieira/ngsDist.git

To install these tools just run:

    % cd ngsDist
    % make
    % make test

Executables are built into the main directory. If you wish to clean all binaries and intermediate files:

    % make clean

### Usage

    % ./ngsDist [options] --geno /path/to/input/file --n_ind INT --n_sites INT --out /path/to/output/file

#### Parameters
* `--geno FILE`: input file with genotypes, genotype likelihoods or genotype posterior probabilities.
* `--n_ind INT`: sample size (number of individuals).
* `--n_sites INT`: number of sites in input file.
* `--tot_sites INT`: total number of sites in dataset.
* `--labels FILE`: labels, one per line, of the input sequences.
* `--probs`: is the input genotype probabilities (likelihoods or posteriors)?
* `--log_scale`: Ii the input in log-scale?.
* `--call_geno`: call genotypes before running analyses.
* `--N_thresh DOUBLE`: minimum threshold to consider site; missing data if otherwise (assumes -call_geno) 
* `--call_thresh DOUBLE`: minimum threshold to call genotype; left as is if otherwise (assumes -call_geno)
* `--pairwise_del`: pairwise deletion of missing data.
* `--avg_nuc_dist`: use average number of nucleotide differences as distance (by default, `ngsDist` uses genotype distances based on allele frequency differences). Only pairs of heterozygous positions are actually affected when using this option, with their distance being 0.5 (instead of 0 by default).
* `--indep_geno`: assume independence between genotypes?
* `--n_boot_rep INT`: number of bootstrap replicates [0].
* `--boot_block_size INT`: block size for bootstrapping [1].
* `--out FILE`: output file name.
* `--n_threads INT`: number of threads to use. [1]
* `--version`: prints program version and exits.
* `--verbose INT`: selects verbosity level. [1]
* `--seed INT`: random number generator seed (only for the bootstrap analysis).

### Input data
As input, `ngsDist` accepts both genotypes, genotype likelihoods (GP) or genotype posterior probabilities (GP). Genotypes must be input as gziped TSV with sites on rows, individuals on columns (__n_sites\*n_ind__) and genotypes coded as [-1, 0, 1, 2]. The file can have a header and an arbitrary number of columns preceeding the actual data (that will all be ignored), much like the Beagle file format ([link](http://faculty.washington.edu/browning/beagle/beagle.html)).
As for GL and GP, `ngsDist` accepts both gzipd TSV and binary formats, but with 3 columns per individual (__3\*n_sites\*n_ind__) and, in the case of the binary, the GL/GP coded as doubles

### Bootstrap Trees
If you want branch support values on your tree, you can use `ngsDist` with the option `--n_boot_rep` and `--boot_block_size` to bootstrap the input data. `ngsDist` will output one distance matrix (the first) for the input dataset, plus `--n_boot_rep` matrices for each of the replicates. After, you can use any program (e.g. [FastME](http://atgc.lirmm.fr/fastme/)) to infer a tree for each of the matrices:

    fastme -T 20 -i testA_8B.dist -s -D 6 -o testA_8B.nwk

split the input dataset tree from the bootstraped ones:

    head -n 1 testA_8B.nwk > testA_8B.main.nwk
    tail -n +2 testA_8B.nwk > testA_8B.boot.nwk

and use [RAxML](https://github.com/stamatak/standard-RAxML) to place supports on the main tree:

    raxmlHPC -f i -t testA_8B.main.nwk -z testA_8B.boot.nwk -m GTRCAT -n testA_8B

### Thread pool
The thread pool	implementation was adapted from Mathias Brossard's and is freely available from:
https://github.com/mbrossard/threadpool
