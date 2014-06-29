#pragma once

#include "shared.hpp"

using namespace std;

extern bool SIG_COND;

const uint64_t N_GENO = 3;
const uint64_t BUFF_LEN = 100000;

// Struct to store all input arguments //GZIP
typedef struct {
  char *in_geno;
  bool in_bin;
  double ***in_geno_lkl;
  double ***geno_lkl;
  uint64_t n_ind;
  uint64_t n_sites;
  char *in_labels;
  char **ind_labels;
  bool in_probs;
  bool in_logscale;
  bool call_geno;
  double score[N_GENO][N_GENO];
  uint64_t n_boot_rep;
  uint64_t boot_block_size;
  char* out_prefix;
  uint n_threads;
  bool version;
  uint verbose;
  uint seed;
  gsl_rng *rnd_gen;
} params;


// Pthread structure
typedef struct {
  params *pars;
  double **dist_matrix;
  uint64_t i1;
  uint64_t i2;
} pth_struct;


double gen_dist(params*, uint64_t, uint64_t);
void gen_dist_slave(void*);

// parse_args.cpp
void init_pars(params* );
int parse_cmd_args(int, char**, params*);
