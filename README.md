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
* `-log_scale`: Is the input in log-scale?.
* `-labels FILE`: Labels, one per line, of the input sequences.
* `-n_ind INT`: Sample size (number of individuals).
* `-n_sites INT`: Total number of sites.
* `-n_boot_rep INT`: Number of bootstrap replicates [0].
* `-call_geno`: Call genotypes before running analyses.
* `-out_prefix FILE`: Output file name.
* `-n_threads INT`: Number of threads to use. [1]
* `-version`: Prints program version and exits.
* `-verbose INT`: Selects verbosity level. [1]
* `-seed INT`: Random number generator seed.

### Input data
As input `ngsDist` needs a Genotype Likelihood (GL) file, formatted as 3*n_ind*n_sites, either in text (gziped) of binary doubles.
