#include <getopt.h>
#include "ngsDist.hpp"


void init_pars(params *pars) {
  pars->in_pp = NULL;
  pars->in_bin = false;
  pars->in_log = false;
  pars->post_prob = NULL;
  pars->in_labels = NULL;
  pars->ind_labels = NULL;
  pars->n_ind = 0;
  pars->n_sites = 0;
  pars->call_geno = false;
  // Distance score matrix based on Eq 8.1 from Del Vecchyo et al. 2013
  pars->score[0][0] = pars->score[1][1] = pars->score[2][2] = 0; 
  pars->score[0][1] = pars->score[1][0] = pars->score[1][2] = pars->score[2][1] = 0.5; 
  pars->score[0][2] = pars->score[2][0] = 1; 
  pars->out_prefix = NULL;
  pars->n_threads = 1;
  pars->n_chunks = 0;
  pars->max_chunk_size = 10000;
  pars->version = false;
  pars->verbose = 1;
  pars->seed = time(NULL);
}


// Parses command line args and stores them into struct params
int parse_cmd_args(int argc, char** argv, params* pars) {
  
  static struct option long_options[] =
    {
      {"geno", required_argument, NULL, 'g'},
      {"log", no_argument, NULL, 'L'},
      {"labels", required_argument, NULL, 'l'},
      {"n_ind", required_argument, NULL, 'n'},
      {"n_sites", required_argument, NULL, 's'},
      {"call_geno", no_argument, NULL, 'G'},
      {"nuc_diff", no_argument, NULL, 'd'},
      {"out", required_argument, NULL, 'o'},
      {"n_threads", required_argument, NULL, 'x'},
      {"chunk_size", required_argument, NULL, 'c'},
      {"version", no_argument, NULL, 'v'},
      {"verbose", required_argument, NULL, 'V'},
      {"seed", required_argument, NULL, 'S'},
      {0, 0, 0, 0}
    };
  
  int c = 0;
  while ( (c = getopt_long_only(argc, argv, "g:Ll:n:s:Gdo:x:c:vV:S:", long_options, NULL)) != -1 )
    switch (c) {
    case 'g':
      pars->in_pp = optarg;
      break;
    case 'L':
      pars->in_log = true;
      break;
    case 'l':
      pars->in_labels = optarg;
      break;
    case 'n':
      pars->n_ind = atoi(optarg);
      break;
    case 's':
      pars->n_sites = atoi(optarg);
      break;
    case 'G':
      pars->call_geno = true;
      break;
    case 'd':
      pars->score[1][1] = 0.5;
      break;
    case 'o':
      pars->out_prefix = optarg;
      break;
    case 'x':
      pars->n_threads = atoi(optarg);
      break;
    case 'c':
      pars->max_chunk_size = atoi(optarg);
      break;
    case 'v':
      pars->version = true;
      break;
    case 'V':
      pars->verbose = atoi(optarg);
      break;
    case 'S':
      pars->seed = atoi(optarg);
      break;
    default:
      exit(-1);
    }

  return 0;
}
