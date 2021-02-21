
#include <getopt.h>
#include "ngsDist.hpp"


void init_pars(params *pars) {
  pars->in_geno = NULL;
  pars->in_bin = false;
  pars->in_probs = false;
  pars->in_logscale = false;
  pars->n_ind = 0;
  pars->n_sites = 0;
  pars->tot_sites = 0;
  pars->in_labels = NULL;
  pars->in_labels_header = false;
  pars->in_pos = NULL;
  pars->in_pos_header = false;
  pars->call_geno = false;
  pars->N_thresh = 0;
  pars->call_thresh = 0;
  pars->pairwise_del = false;
  // Distance score matrix based on:
  // Gronau I    et al. 2011 (eq 12)
  // Freedman AH et al. 2014 (eq 8.1)
  pars->score[0][0] = pars->score[1][1] = pars->score[2][2] = 0; 
  pars->score[0][1] = pars->score[1][0] = pars->score[1][2] = pars->score[2][1] = 0.5; 
  pars->score[0][2] = pars->score[2][0] = 1; 
  pars->evol_model = 1;
  pars->indep_geno = false;
  pars->n_boot_rep = 0;
  pars->boot_block_size = 1;
  pars->out = NULL;
  pars->n_threads = 1;
  pars->verbose = 1;
  pars->seed = time(NULL);
  pars->rnd_gen = NULL;
}


// Supported evolutionary models
const char* evol_model[] = {"Raw p-distance", // Raw p-distance
			    "Log transf. p-distance", // Logarithmic transformed p-distance
			    "JC69", 
			    "K80", // JC69 + (transiton rates != transversions rates)
			    "F81", // JC69 + (base freqs != 0.25)
			    "HKY85/F84", // K80 + F81 
			    "TN93", // HKY85/F84 + (A/G != T/C)
};


// Parses command line args and stores them into struct params
void parse_cmd_args(params* pars, int argc, char** argv) {
  
  static struct option long_options[] =
    {
      {"geno", required_argument, NULL, 'g'},
      {"probs", no_argument, NULL, 'p'},
      {"log_scale", no_argument, NULL, 'l'},
      {"n_ind", required_argument, NULL, 'n'},
      {"n_sites", required_argument, NULL, 's'},
      {"tot_sites", required_argument, NULL, 'S'},
      {"labels", required_argument, NULL, 'L'},
      {"labelsH", required_argument, NULL, 'H'},
      {"pos", required_argument, NULL, 'a'},
      {"posH", required_argument, NULL, 'A'},
      {"call_geno", no_argument, NULL, 'c'},
      {"N_thresh", required_argument, NULL, 'N'},
      {"call_thresh", required_argument, NULL, 'C'},
      {"pairwise_del", no_argument, NULL, 'D'},
      {"avg_nuc_dist", no_argument, NULL, 'd'},
      {"evol_model", required_argument, NULL, 'm'},
      {"indep_geno", no_argument, NULL, 'I'},
      {"n_boot_rep", required_argument, NULL, 'b'},
      {"boot_block_size", required_argument, NULL, 'B'},
      {"out", required_argument, NULL, 'o'},
      {"n_threads", required_argument, NULL, 'x'},
      {"verbose", required_argument, NULL, 'V'},
      {"seed", required_argument, NULL, 'r'},
      {0, 0, 0, 0}
    };
  
  int c = 0;
  while ( (c = getopt_long_only(argc, argv, "g:pln:s:S:a:A:L:H:cN:C:Ddm:Ib:B:o:x:V:r:", long_options, NULL)) != -1 )
    switch (c) {
    case 'g':
      pars->in_geno = optarg;
      break;
    case 'p':
      pars->in_probs = true;
      break;
    case 'l':
      pars->in_logscale = true;
      pars->in_probs = true;
      break;
    case 'n':
      pars->n_ind = atol(optarg);
      break;
    case 's':
      pars->n_sites = atol(optarg);
      break;
    case 'S':
      pars->tot_sites = atol(optarg);
      break;
    case 'L':
      pars->in_labels = optarg;
      pars->in_labels_header = false;
      break;
    case 'H':
      pars->in_labels = optarg;
      pars->in_labels_header = true;
      break;
    case 'a':
      pars->in_pos = optarg;
      pars->in_pos_header = false;
      break;
    case 'A':
      pars->in_pos = optarg;
      pars->in_pos_header = true;
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
    case 'D':
      pars->pairwise_del = true;
      break;
    case 'd':
      // Freedman AH et al. 2014 (eq 8.2)
      pars->score[1][1] = 0.5;
      break;
    case 'm':
      pars->evol_model = atol(optarg);
      break;
    case 'I':
      pars->indep_geno = true;
      break;
    case 'b':
      pars->n_boot_rep = atol(optarg);
      break;
    case 'B':
      pars->boot_block_size = atol(optarg);
      break;
    case 'o':
      pars->out = optarg;
      break;
    case 'x':
      pars->n_threads = atoi(optarg);
      break;
    case 'V':
      pars->verbose = atoi(optarg);
      break;
    case 'r':
      pars->seed = atoi(optarg);
      break;
    default:
      exit(-1);
    }


  if(pars->verbose >= 1) {
    fprintf(stderr, "==> Input Arguments:\n");
    fprintf(stderr, "\tgeno: %s\n\tprobs: %s\n\tlog_scale: %s\n\tn_ind: %lu\n\tn_sites: %lu\n\ttot_sites: %lu\n\tlabels: %s (%s header)\n\tpositions: %s (%s header)\n\tcall_geno: %s\n\tN_thresh: %f\n\tcall_thresh: %f\n\tpairwise_del: %s\n\tavg_nuc_dist: %s\n\tevol_model: %s\n\tgeno_indep: %s\n\tn_boot_rep: %lu\n\tboot_block_size: %lu\n\tout: %s\n\tn_threads: %d\n\tverbose: %d\n\tseed: %d\n\tversion: %s (%s @ %s)\n\n",
	    pars->in_geno,
	    pars->in_probs ? "true":"false",
	    pars->in_logscale ? "true":"false",
	    pars->n_ind,
	    pars->n_sites,
	    pars->tot_sites,
	    pars->in_labels,
	    pars->in_labels_header ? "WITH" : "WITHOUT",
	    pars->in_pos,
	    pars->in_pos_header ? "WITH" : "WITHOUT",
	    pars->call_geno ? "true":"false",
	    pars->N_thresh,
	    pars->call_thresh,
	    pars->pairwise_del ? "true":"false",
	    pars->score[1][1] == 0.5 ? "true":"false",
	    evol_model[pars->evol_model],
	    pars->indep_geno ? "true":"false",
	    pars->n_boot_rep,
	    pars->boot_block_size,
	    pars->out,
	    pars->n_threads,
	    pars->verbose,
	    pars->seed,
	    version, __DATE__, __TIME__);
  }
  if(pars->verbose > 4)
    fprintf(stderr, "==> Verbose values greater than 4 for debugging purpose only. Expect large amounts of info on screen\n");



  /////////////////////
  // Check Arguments //
  /////////////////////
  if(pars->in_geno == NULL)
    error(__FUNCTION__, "genotype input file (--geno) missing!");
  if(pars->n_ind == 0)
    error(__FUNCTION__, "number of individuals (--n_ind) missing!");
  if(pars->n_sites == 0)
    error(__FUNCTION__, "number of sites (--n_sites) missing!");
  if(pars->tot_sites > 0 && pars->pairwise_del)
    error(__FUNCTION__, "cannot specify total number of sites (--tot_sites) with pairwise deletion (--pairwise_del)!");
  if(pars->call_geno && !pars->in_probs)
    error(__FUNCTION__, "can only call genotypes from likelihoods/probabilities!");
  if(pars->evol_model < 0 || pars->evol_model > 6)
    error(__FUNCTION__, "invalid correction method specified!");
  if(pars->evol_model > 2 && pars->in_pos == NULL)
    error(__FUNCTION__, "use of more complex evolutionary models requires position information!");
  if(pars->out == NULL)
    error(__FUNCTION__, "output prefix (--out) missing!");
  if(pars->n_threads < 1)
    error(__FUNCTION__, "number of threads cannot be less than 1!");
}
