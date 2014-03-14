#include <getopt.h>
#include "ngsDist.hpp"


void init_pars(params *pars) {
  pars->in_geno = NULL;
  pars->in_bin = false;
  pars->in_lkl = false;
  pars->in_log = false;
  pars->geno_lkl = NULL;
  pars->in_labels = NULL;
  pars->ind_labels = NULL;
  pars->n_ind = 0;
  pars->n_sites = 0;
  pars->call_geno = false;
  // Distance score matrix based on Eq 8.1 from Gronau et al 2011 and Del Vecchyo et al. 2014
  pars->score[0][0] = pars->score[1][1] = pars->score[2][2] = 0; 
  pars->score[0][1] = pars->score[1][0] = pars->score[1][2] = pars->score[2][1] = 0.5; 
  pars->score[0][2] = pars->score[2][0] = 1; 
  pars->out_prefix = NULL;
  pars->n_threads = 1;
  pars->version = false;
  pars->verbose = 1;
  pars->seed = time(NULL);
}


// Parses command line args and stores them into struct params
int parse_cmd_args(int argc, char** argv, params* pars) {
  
  static struct option long_options[] =
    {
      {"geno", required_argument, NULL, 'g'},
      {"lkl", no_argument, NULL, 'l'},
      {"log", no_argument, NULL, 'L'},
      {"labels", required_argument, NULL, 'y'},
      {"n_ind", required_argument, NULL, 'n'},
      {"n_sites", required_argument, NULL, 's'},
      {"call_geno", no_argument, NULL, 'G'},
      {"nuc_diff", no_argument, NULL, 'd'},
      {"out", required_argument, NULL, 'o'},
      {"n_threads", required_argument, NULL, 'x'},
      {"version", no_argument, NULL, 'v'},
      {"verbose", required_argument, NULL, 'V'},
      {"seed", required_argument, NULL, 'S'},
      {0, 0, 0, 0}
    };
  
  int c = 0;
  while ( (c = getopt_long_only(argc, argv, "g:lLy:n:s:Gdo:x:vV:S:", long_options, NULL)) != -1 )
    switch (c) {
    case 'g':
      pars->in_geno = optarg;
      break;
    case 'l':
      pars->in_lkl = true;
      break;
    case 'L':
      pars->in_log = true;
      break;
    case 'y':
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
