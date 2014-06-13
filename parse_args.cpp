#include <getopt.h>
#include "ngsDist.hpp"


void init_pars(params *pars) {
  pars->in_geno = NULL;
  pars->in_bin = false;
  pars->n_ind = 0;
  pars->n_sites = 0;
  pars->in_labels = NULL;
  pars->in_probs = false;
  pars->in_logscale = false;
  pars->call_geno = false;
  // Distance score matrix based on Eq 8.1 from Gronau et al 2011 and Del Vecchyo et al. 2014
  pars->score[0][0] = pars->score[1][1] = pars->score[2][2] = 0; 
  pars->score[0][1] = pars->score[1][0] = pars->score[1][2] = pars->score[2][1] = 0.5; 
  pars->score[0][2] = pars->score[2][0] = 1; 
  pars->n_boot_rep = 0;
  pars->boot_block_size = 1;
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
      {"n_ind", required_argument, NULL, 'n'},
      {"n_sites", required_argument, NULL, 's'},
      {"labels", required_argument, NULL, 'L'},
      {"probs", no_argument, NULL, 'p'},
      {"log_scale", no_argument, NULL, 'l'},
      {"call_geno", no_argument, NULL, 'G'},
      {"alt_het_diff", no_argument, NULL, 'd'},
      {"n_boot_rep", required_argument, NULL, 'b'},
      {"boot_block_size", required_argument, NULL, 'B'},
      {"out_prefix", required_argument, NULL, 'o'},
      {"n_threads", required_argument, NULL, 'x'},
      {"version", no_argument, NULL, 'v'},
      {"verbose", required_argument, NULL, 'V'},
      {"seed", required_argument, NULL, 'S'},
      {0, 0, 0, 0}
    };
  
  int c = 0;
  while ( (c = getopt_long_only(argc, argv, "g:n:s:L:plGdb:B:o:x:vV:S:", long_options, NULL)) != -1 )
    switch (c) {
    case 'g':
      pars->in_geno = optarg;
      break;
    case 'n':
      pars->n_ind = atoi(optarg);
      break;
    case 's':
      pars->n_sites = atoi(optarg);
      break;
    case 'L':
      pars->in_labels = optarg;
      break;
    case 'p':
      pars->in_probs = true;
      break;
    case 'l':
      pars->in_logscale = true;
      break;
    case 'G':
      pars->call_geno = true;
      break;
    case 'd':
      pars->score[1][1] = 0.5;
      break;
    case 'b':
      pars->n_boot_rep = atoi(optarg);
      break;
    case 'B':
      pars->boot_block_size = atoi(optarg);
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
