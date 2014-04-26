#pragma once

#include <pthread.h>
#include <semaphore.h>

#include "shared.hpp"

using namespace std;

extern bool SIG_COND;

const uint64_t N_GENO = 3;
const uint64_t BUFF_LEN = 100000;

// Struct to store all input arguments //GZIP
typedef struct {
  char* in_geno;
  bool in_bin;
  bool in_lkl;
  bool in_log;
  double*** geno_lkl;
  char* in_labels;
  char** ind_labels;
  uint64_t n_ind;
  uint64_t n_sites;
  bool call_geno;
  double score[N_GENO][N_GENO];
  char* out_prefix;
  uint n_threads;
  bool version;
  uint verbose;
  uint seed;
  sem_t pth_sem;
} params;


// Pthread structure
typedef struct {
  pthread_t id;
  pthread_attr_t attr;
  params* pars;
  double** dist_matrix;
  uint64_t i1;
  uint64_t i2;
} pth_struct;


double gen_dist(params*, uint64_t, uint64_t);
void* gen_dist_slave(void*);

// parse_args.cpp
void init_pars(params* );
int parse_cmd_args(int, char**, params*);

// read_data.cpp
int read_geno(params*);
int read_labels(params*);
