#pragma once

#include "read_data.hpp"
#include "threadpool.h"

extern bool SIG_COND;
extern char const* version;


// Struct to store all input arguments
typedef struct {
  char *in_geno;
  bool in_bin;
  bool in_probs;
  bool in_logscale;
  uint64_t n_ind;
  uint64_t n_sites;
  uint64_t tot_sites;
  char *in_labels;
  bool in_labels_header;
  char *in_pos;
  bool in_pos_header;
  bool call_geno;
  double N_thresh;
  double call_thresh;
  bool pairwise_del;
  double score[N_GENO][N_GENO];
  uint64_t evol_model;
  bool indep_geno;
  uint64_t n_boot_rep;
  uint64_t boot_block_size;
  char* out;
  uint n_threads;
  uint verbose;
  uint seed;
  gsl_rng *rnd_gen;

  char **ind_labels;     // n_ind * BUFF_LEN
  double ***in_geno_lkl; // n_ind * n_sites+1 * N_GENO
  double ***geno_lkl;    // n_ind * n_sites+1 * N_GENO
  char ***pos;           // n_sites+1 * N_TSV_FIELDS * BUFF_LEN

  threadpool_t *thread_pool;
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
void parse_cmd_args(params*, int, char**);
