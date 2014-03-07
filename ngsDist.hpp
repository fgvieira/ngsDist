#pragma once

#include <algorithm>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <zlib.h>
#include <sys/mman.h>

#include "shared.hpp"

using namespace std;

extern bool SIG_COND;

const uint64_t N_GENO = 3;
const uint64_t BUFF_LEN = 100000;


// Struct to store all input arguments //GZIP
typedef struct {
  char* in_geno;
  bool in_bin;
  bool in_log;
  double*** post_prob;
  uint64_t n_ind;
  uint64_t n_sites;
  bool call_geno;
  double score[N_GENO][N_GENO];
  char* out_prefix;
  uint n_threads;
  uint n_chunks;
  uint max_chunk_size;
  bool version;
  uint verbose;
  uint seed;
} params;

double gen_dist(params*, uint64_t, uint64_t);

// parse_args.cpp
void init_pars(params* );
int parse_cmd_args(int, char**, params*);

// read_data.cpp
int read_geno(params*);
uint64_t read_chunk(double**, params*, uint64_t);
