# ngsDist

`ngsDist` is a program to estimate pairwise genetic distances directly, taking the uncertainty of genotype's assignation into account. It does so by avoiding genotype calling and using genotype likelihoods or posterior probabilities.

### Citation

`ngsDist` was published in 2015 at [Biological Journal of the Linnean Society](http://onlinelibrary.wiley.com/doi/10.1111/bij.12511/abstract), so please cite it if you use it in your work:

    Vieira FG, Lassalle F, Korneliussen TS, Fumagalli M
    Improving the estimation of genetic distances from Next-Generation Sequencing data
    Biological Journal of the Linnean Society (2015) doi: 10.1111/bij.12511

### Installation

`ngsDist` can be easily installed but has some external dependencies:

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

    % ./ngsDist [options] --geno glf/in/file --n_ind INT --n_sites INT --out_prefix output/file

#### Parameters
* `--geno FILE`: Input GL file.
* `--n_ind INT`: Sample size (number of individuals).
* `--n_sites INT`: Total number of sites.
* `--labels FILE`: Labels, one per line, of the input sequences.
* `--probs`: is the input genotype probabilities (likelihoods or posteriors)?
* `--log_scale`: Is the input in log-scale?.
* `--call_geno`: Call genotypes before running analyses.
* `--N_thresh DOUBLE`: minimum threshold to consider site; missing data if otherwise (assumes -call_geno) 
* `--call_thresh DOUBLE`: minimum threshold to call genotype; left as is if otherwise (assumes -call_geno)
* `--het_dist`: Use alternative heterozygote distance [0].
* `--indep_geno`: Assume independence between genotypes?
* `--n_boot_rep INT`: Number of bootstrap replicates [0].
* `--boot_block_size INT`: Block size for bootstrapping [1].
* `--out_prefix FILE`: Output file name.
* `--n_threads INT`: Number of threads to use. [1]
* `--version`: Prints program version and exits.
* `--verbose INT`: Selects verbosity level. [1]
* `--seed INT`: Random number generator seed.

### Input data
As input `ngsDist` needs a Genotype Likelihood (GL) file, formatted as __3\*n_ind\*n_sites__, either as gziped TSV or binary doubles.

### Thread pool
The thread pool	implementation was adapted from Mathias Brossard's and is freely available from:
https://github.com/mbrossard/threadpool
