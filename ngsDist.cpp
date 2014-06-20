
/*
 *
 * ngsDist - NGS data individual inbreeding coefficients estimation.
 * Copyright (C) 2012  Filipe G. Vieira
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/

#include "ngsDist.hpp"
#include "read_data.cpp"
#include "emOptim2.cpp"
#include "threadpool.c"

char const* version = "0.0.1b";


int main (int argc, char** argv) {
  /////////////////////
  // Parse Arguments //
  /////////////////////
  params* pars = new params;
  init_pars(pars);
  parse_cmd_args(argc, argv, pars);

  if( pars->version ) {
    printf("ngsDist v%s\nCompiled on %s @ %s", version, __DATE__, __TIME__);
    exit(0);
  }

  if(pars->verbose >= 1) {
    printf("==> Input Arguments:\n");
    printf("\tgeno: %s\n\tn_ind: %lu\n\tn_sites: %lu\n\tlabels: %s\n\tprobs: %s\n\tlog_scale: %s\n\tcall_geno: %s\n\thet_dist: %f\n\tn_boot_rep: %lu\n\tboot_block_size: %lu\n\tout_prefix: %s\n\tn_threads: %d\n\tversion: %s\n\tverbose: %d\n\tseed: %d\n\n",
	   pars->in_geno,
	   pars->n_ind,
           pars->n_sites,
	   pars->in_labels,
	   pars->in_probs ? "true":"false",
	   pars->in_logscale ? "true":"false",
	   pars->call_geno ? "true":"false",
	   pars->score[1][1],
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
  
  
  
  ///////////////////////
  // Adjust Parameters //
  ///////////////////////
  // Calculate total number of combinations
  uint64_t n_comb = (pow(pars->n_ind, 2) - pars->n_ind) / 2;
  if(pars->verbose >= 1)
    printf("==> Analysis will be run in %lu combinations\n", n_comb);
  // Adjust thread number to chunks
  if(n_comb < pars->n_threads){
    if(pars->verbose >= 1)
      printf("==> Fewer combinations (%ld) than threads (%d). Reducing the number of threads...\n", n_comb, pars->n_threads);
    pars->n_threads = n_comb;
  }



  ///////////////////////
  // Check input files //
  ///////////////////////
  // Get file total size
  struct stat st;
  if( stat(pars->in_geno, &st) != 0 )
    error(__FUNCTION__, "cannot check file size!");
  if( strcmp(strrchr(pars->in_geno, '.'), ".gz") == 0 ){
    if(pars->verbose >= 1)
      printf("==> GZIP input file (never BINARY)\n");
    pars->in_bin = false;
  }else if( pars->n_sites == st.st_size/sizeof(double)/pars->n_ind/N_GENO ){
    if(pars->verbose >= 1)
      printf("==> BINARY input file\n");
    pars->in_bin = true;
  }else
    error(__FUNCTION__, "invalid/corrupt genotype input file!");
  


  ////////////////////////////
  // Prepare initial values //
  ////////////////////////////
  // Read labels file
  if(pars->in_labels){
    if(pars->verbose >= 1)
      printf("==> Reading labels\n");

    int64_t ret = read_file(pars->in_labels, &pars->ind_labels, 1000);
    if(pars->verbose >= 1)
      printf("> Found %ld labels in file\n", ret);

    if(ret == -1)
      error(__FUNCTION__, "cannot open labels file!");
    else if(ret != (int64_t) pars->n_ind)
      error(__FUNCTION__, "wrong number of labels provided!");
  }else{
    pars->ind_labels = init_char(pars->n_ind, BUFF_LEN, (const char*) "Ind_#");
    // Tweak initiation value (replace # by number)
    for(uint64_t i = 0; i < pars->n_ind; i++){
      char* pch = strchr(pars->ind_labels[i], '#');
      if(pch)
	sprintf(pch, "%lu", i);
    }
  }
  
  if(pars->verbose >= 4)
    for(uint64_t i = 0; i < pars->n_ind; i++)
      printf("%s\n", pars->ind_labels[i]);

  // Read from GENO file
  if(pars->verbose >= 1)
    printf("==> Reading genotype posterior probabilities\n");
  pars->in_geno_lkl = read_geno(pars->in_geno, pars->in_bin, pars->in_probs, pars->n_ind, pars->n_sites);
  
  // Initialize geno_lkl pointers
  pars->geno_lkl = new double**[pars->n_ind];
  for(uint64_t i = 0; i < pars->n_ind; i++)
    pars->geno_lkl[i] = new double*[pars->n_sites+1];

  for(uint64_t i = 0; i < pars->n_ind; i++)
    for(uint64_t s = 1; s <= pars->n_sites; s++){
      // Call genotypes
      if(pars->call_geno)
	call_geno(pars->in_geno_lkl[i][s], N_GENO, pars->in_logscale);
      // Convert space
      if(pars->in_logscale)
	conv_space(pars->in_geno_lkl[i][s], N_GENO, exp);
    }

  // Initialize random number generator
  if(pars->verbose >= 2)
    printf("==> Setting seed for random number generator\n");
  gsl_rng* rnd_gen = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rnd_gen, pars->seed);

 

  //////////////////////
  // Open output file //
  //////////////////////
  FILE* out_fh = fopen(strcat(pars->out_prefix, ".dist"), "w");
  if(out_fh == NULL)
    error(__FUNCTION__, "cannot open output file!");
 

  //////////////////
  // Analyze Data //
  //////////////////
  fflush(stdout);
  uint64_t comb_id = 0;
  double** dist_matrix = init_double(pars->n_ind, pars->n_ind, 0);

  // Create pthread structure array
  pth_struct* pth = new pth_struct[n_comb];
  // Initialize pars and distance matrix pointers
  for(comb_id = 0; comb_id < n_comb; comb_id++){
    pth[comb_id].pars = pars;
    pth[comb_id].dist_matrix = dist_matrix;
  }

  // Loop for bootstrap analyses
  for(uint64_t rep = 0; rep <= pars->n_boot_rep; rep++){
    comb_id = 0;

    if(pars->verbose >= 1){
      if(rep == 0)
	// Full dataset analyses
	printf("==> Analyzing full dataset...\n");
      else
	// Bootstrap analyses
	printf("==> Bootstrap replicate # %lu ...\n", rep);
    }


    /* Create a threadpool of thread workers. */
    struct threadpool *thread_pool;
    if ((thread_pool = threadpool_init(pars->n_threads, n_comb)) == NULL)
      error(__FUNCTION__, "failed to create thread pool!");


    // Map from in_geno_lkl data to geno_lkl
    if(pars->n_sites % pars->boot_block_size != 0)
      error(__FUNCTION__, "bootstrap block size must be a multiple of the total length!");
    uint64_t n_blocks = (uint64_t) ceil(pars->n_sites/pars->boot_block_size);
    
    for(uint64_t b = 0; b < n_blocks; b++){
      uint64_t rnd_b = b;
      if(rep > 0)
        rnd_b = (uint64_t) floor( draw_rnd(rnd_gen, 0, n_blocks) );

      for(uint64_t s = 1; s <= pars->boot_block_size; s++){
	uint64_t orig =     b * pars->boot_block_size + s;
	uint64_t rnd  = rnd_b * pars->boot_block_size + s;

	if(pars->verbose > 5)
	  printf("block: %lu\torig_site: %lu\trand_block:%lu\trand_site: %lu\n", b, orig, rnd_b, rnd);

	for(uint64_t i = 0; i < pars->n_ind; i++)
	  pars->geno_lkl[i][orig] = pars->in_geno_lkl[i][rnd];
      }
    }


    if(pars->verbose >= 2)
      printf("> Calculating pairwise genetic distances\n");

    for(uint64_t i1 = 0; i1 < pars->n_ind; i1++)
      for(uint64_t i2 = i1+1; i2 < pars->n_ind; i2++){
	// Set which individuals to analyze
	pth[comb_id].i1 = i1;
	pth[comb_id].i2 = i2;

	// Add task to thread pool
	int ret = threadpool_add_task(thread_pool, gen_dist_slave, (void*) &pth[comb_id], 0);
	if (ret == -1)
	  error(__FUNCTION__, "error while adding task to thread pool!");
	else if (ret == -2)
	  error(__FUNCTION__, "task did not execute since the pool was overloaded!");

	comb_id++;

	if(pars->verbose >= 5)
	  printf("> Launched thread for individuals %lu and %lu (comb# %lu).\n", i1, i2, comb_id);
      }


  
    //////////////////////////
    // Wait for all threads //
    //////////////////////////
    threadpool_free(thread_pool,1);


    if(n_comb != comb_id)
      error(__FUNCTION__, "some combinations are missing!");


    ///////////////////////////
    // Print Distance Matrix //
    ///////////////////////////
    if(pars->verbose >= 2)
      printf("> Printing distance matrix\n");

    fprintf(out_fh, "\n  %lu\n", pars->n_ind);
    for(uint64_t i = 0; i < pars->n_ind; i++){
      char* buf = join(dist_matrix[i], pars->n_ind, "\t");
      fprintf(out_fh, "%s\t%s\n", pars->ind_labels[i], buf);
      delete [] buf;
    }

  }
  fclose(out_fh);



  /////////////////
  // Free Memory //
  /////////////////
  if(pars->verbose >= 1)
    printf("====> Freeing memory...\n");

  free_ptr((void**) dist_matrix, pars->n_ind);
  delete [] pth;
  // pars struct
  free_ptr((void***) pars->in_geno_lkl, pars->n_ind, pars->n_sites+1);
  free_ptr((void**) pars->ind_labels, pars->n_ind);
  //free_ptr((void*) pars->in_geno);

  if(pars->verbose >= 1)
    printf("Done!\n");
  delete pars;

  gsl_rng_free(rnd_gen);

  return 0;
}




double gen_dist(params* p, uint64_t i1, uint64_t i2){
  double dist = 0;

  Matrix<double> GL1 = alloc(1,N_GENO);
  Matrix<double> GL2 = alloc(1,N_GENO);
  int dim = GL1.y*GL2.y;

  for(uint64_t s = 1; s <= p->n_sites; s++){
    double* sfs = init_double(dim, (double) 1/dim);
    GL1.mat[0][0] = p->geno_lkl[i1][s][0];
    GL1.mat[0][1] = p->geno_lkl[i1][s][1];
    GL1.mat[0][2] = p->geno_lkl[i1][s][2];
    GL2.mat[0][0] = p->geno_lkl[i2][s][0];
    GL2.mat[0][1] = p->geno_lkl[i2][s][1];
    GL2.mat[0][2] = p->geno_lkl[i2][s][2];
    em2(sfs, &GL1, &GL2, 0.001, 50, dim);

    for(uint64_t g1 = 0; g1 < N_GENO; g1++)
      for(uint64_t g2 = 0; g2 < N_GENO; g2++)
	dist += p->score[g1][g2] * sfs[3*g1+g2];

    free_ptr(sfs);
  }

  dalloc(GL1, 1);
  dalloc(GL2, 1);
  return dist/p->n_sites;
}



void gen_dist_slave(void* pth){
  pth_struct* p = (pth_struct*) pth;

  p->dist_matrix[p->i1][p->i2] = p->dist_matrix[p->i2][p->i1] = gen_dist(p->pars, p->i1, p->i2);
}
