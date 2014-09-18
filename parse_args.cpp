
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
  pars->N_thresh = 0;
  pars->call_thresh = 0;
  // Distance score matrix based on Eq 8.1 from Gronau et al 2011 and Del Vecchyo et al. 2014
  pars->score[0][0] = pars->score[1][1] = pars->score[2][2] = 0; 
  pars->score[0][1] = pars->score[1][0] = pars->score[1][2] = pars->score[2][1] = 0.5; 
  pars->score[0][2] = pars->score[2][0] = 1; 
  pars->indep_geno = false;
  pars->n_boot_rep = 0;
  pars->boot_block_size = 1;
  pars->out_prefix = NULL;
  pars->n_threads = 1;
  pars->version = false;
  pars->verbose = 1;
  pars->seed = time(NULL);
  pars->rnd_gen = NULL;
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
      {"call_geno", no_argument, NULL, 'c'},
      {"N_thresh", required_argument, NULL, 'N'},
      {"call_thresh", required_argument, NULL, 'C'},
      {"alt_het_diff", no_argument, NULL, 'd'},
      {"indep_geno", no_argument, NULL, 'I'},
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
  while ( (c = getopt_long_only(argc, argv, "g:n:s:L:plcN:C:dIb:B:o:x:vV:S:", long_options, NULL)) != -1 )
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
    case 'c':
      pars->call_geno = true;
      break;
    case 'N':
      pars->N_thresh = atof(optarg);
      pars->call_geno = true;
      break;
    case 'C':
      pars->call_thresh = atof(optarg);
      pars->call_geno = true;
      break;
    case 'd':
      pars->score[1][1] = 0.5;
      break;
    case 'I':
      pars->indep_geno = true;
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


  if(pars->verbose >= 1) {
    printf("==> Input Arguments:\n");
    printf("\tgeno: %s\n\tn_ind: %lu\n\tn_sites: %lu\n\tlabels: %s\n\tprobs: %s\n\tlog_scale: %s\n\tcall_geno: %s\n\tN_thresh: %f\n\tcall_thresh: %f\n\thet_dist: %f\n\tgeno_indep: %s\n\tn_boot_rep: %lu\n\tboot_block_size: %lu\n\tout_prefix: %s\n\tn_threads: %d\n\tversion: %s\n\tverbose: %d\n\tseed: %d\n\n",
           pars->in_geno,
           pars->n_ind,
           pars->n_sites,
           pars->in_labels,
           pars->in_probs ? "true":"false",
           pars->in_logscale ? "true":"false",
           pars->call_geno ? "true":"false",
           pars->N_thresh,
           pars->call_thresh,
           pars->score[1][1],
           pars->indep_geno ? "true":"false",
           pars->n_boot_rep,
           pars->boot_block_size,
           pars->out_prefix,
           pars->n_threads,
           pars->version ? "true":"false",
           pars->verbose,
           pars->seed);
  }
  if(pars->verbose > 4)
    printf("==> Verbose values greater than 4 for debugging purpose only. Expect large amounts of info on screen\n");



  /////////////////////
  // Check Arguments //
  /////////////////////
  if(pars->in_geno == NULL)
    error(__FUNCTION__, "Genotype input file (-geno) missing!");
  if(pars->n_ind == 0)
    error(__FUNCTION__, "Number of individuals (-n_ind) missing!");
  if(pars->n_sites == 0)
    error(__FUNCTION__, "Number of sites (-n_sites) missing!");
  if(pars->call_geno && !pars->in_probs)
    error(__FUNCTION__, "Can only call genotypes from probabilities!");
  if(pars->n_threads < 1)
    error(__FUNCTION__, "NUmber of threads cannot be less than 1!");


  return 0;
}
