# ngsDist

`ngsDist` is a program to estimate pairwise genetic distances directly from low genotype likelihoods, taking the uncertainty of genotype's assignation into account.


### Installation

To install the entire package just download the source code:

    % git clone https://github.com/fgvieira/ngsDist.git

To install these tools just run:

    % cd ngsDist
    % make
    % make test

Executables are built into the main directory. If you wish to clean all binaries and intermediate files:

    % make clean

### Usage

    % ./ngsDist [options] -n_ind INT -s INT -g glf/in/file -o output/file

#### Parameters
* `-geno FILE`: Input GL file.
* `-n_ind INT`: Sample size (number of individuals).
* `-n_sites INT`: Total number of sites.
* `-labels FILE`: Labels, one per line, of the input sequences.
* `-probs`: is the input genotype probabilities (likelihoods or posteriors)?
* `-log_scale`: Is the input in log-scale?.
* `-call_geno`: Call genotypes before running analyses.
* `-N_thresh DOUBLE`: minimum threshold to consider site; missing data if otherwise (assumes -call_geno) 
* `-call_thresh DOUBLE`: minimum threshold to call genotype; left as is if otherwise (assumes -call_geno)
* `-het_dist`: Use alternative heterozygote distance [0].
* `-n_boot_rep INT`: Number of bootstrap replicates [0].
* `-boot_block_size INT`: Block size for bootstrapping [1].
* `-out_prefix FILE`: Output file name.
* `-n_threads INT`: Number of threads to use. [1]
* `-version`: Prints program version and exits.
* `-verbose INT`: Selects verbosity level. [1]
* `-seed INT`: Random number generator seed.

### Input data
As input `ngsDist` needs a Genotype Likelihood (GL) file, formatted as 3*n_ind*n_sites, either in text (gziped) of binary doubles.

### Thread pool
The thread pool	implementation was adapted from Mathias Brossard's and is freely available from:
https://github.com/mbrossard/threadpool
