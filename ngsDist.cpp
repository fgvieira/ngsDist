
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
#include "emOptim2.cpp"

char const* version = "1.0.0";

void rnd_map_data(params *pars, uint64_t n_blocks);

int main (int argc, char** argv) {
  /////////////////////
  // Parse Arguments //
  /////////////////////
  params* pars = new params;
  init_pars(pars);
  parse_cmd_args(pars, argc, argv);



  ///////////////////////
  // Adjust Parameters //
  ///////////////////////
  // Calculate total number of combinations
  uint64_t n_comb = (pow(pars->n_ind, 2) - pars->n_ind) / 2;
  if(pars->verbose >= 1)
    fprintf(stderr, "==> Analysis will be run in %lu combinations\n", n_comb);
  // Adjust thread number to combinations
  if(n_comb < pars->n_threads){
    if(pars->verbose >= 1)
      fprintf(stderr, "==> Fewer combinations (%ld) than threads (%d). Reducing the number of threads...\n", n_comb, pars->n_threads);
    pars->n_threads = n_comb;
  }
  // If input are genotypes (either called ot not) assume independence between genotypes (faster)
  if(pars->call_geno || !pars->in_probs)
    pars->indep_geno = true;

  if(pars->indep_geno)
    if(pars->verbose >= 1)
      fprintf(stderr, "==> Using faster algorithm (assuming independence of genotypes)!\n");



  ///////////////////////
  // Check input files //
  ///////////////////////
  // Get file total size
  if( strcmp(pars->in_geno, "-") == 0 ){
    if(pars->verbose >= 1)
      fprintf(stderr, "==> Reading from STDIN (BINARY)\n");
    pars->in_bin = true;
  }else{
    struct stat st;
    if( stat(pars->in_geno, &st) != 0 )
      error(__FUNCTION__, "cannot check GENO file size!");

    if( strcmp(strrchr(pars->in_geno, '.'), ".gz") == 0 ){
      if(pars->verbose >= 1)
	fprintf(stderr, "==> GZIP input file (never BINARY)\n");
      pars->in_bin = false;
    }else{
      if(pars->verbose >= 1)
	fprintf(stderr, "==> BINARY input file\n");
      pars->in_bin = true;
      pars->in_probs = true;

      if( pars->n_sites != st.st_size/sizeof(double)/pars->n_ind/N_GENO )
	error(__FUNCTION__, "invalid/corrupt genotype input file!");
    }
  }
  


  ////////////////////////////
  // Prepare initial values //
  ////////////////////////////
  // Read labels file
  if(pars->in_labels){
    if(pars->verbose >= 1)
      fprintf(stderr, "==> Reading labels\n");

    int64_t ret = read_file(pars->in_labels, &pars->ind_labels, 1000);
    if(ret == -1)
      error(__FUNCTION__, "cannot open labels file!");

    if(ret == (int64_t) pars->n_ind +1){
      fprintf(stderr, "> Assuming label file has a header\n");
      ret--;
      char **tmp = pars->ind_labels;
      pars->ind_labels = init_ptr(ret, 0, (const char*) '\0');
      memcpy(pars->ind_labels, tmp+1, ret*sizeof(char*));
      free_ptr((void*) *tmp); // Free first line (the header)
      free_ptr((void**) tmp); // Free remaining pointers
    }
    if(pars->verbose >= 1)
      fprintf(stderr, "> Found %ld labels in file\n", ret);
    if(ret != (int64_t) pars->n_ind)
      error(__FUNCTION__, "wrong number of labels provided!");

    // Fix labels...
    char* ptr;
    for(int64_t i = 0; i < ret; i++){
      ptr = strchr(pars->ind_labels[i], '\t');
      if(ptr != NULL)
	*ptr = '\0';
    }
  }else{
    pars->ind_labels = init_ptr(pars->n_ind, BUFF_LEN, (const char*) "Ind_#");
    // Tweak initiation value (replace # by number)
    for(uint64_t i = 0; i < pars->n_ind; i++){
      char* pch = strchr(pars->ind_labels[i], '#');
      if(pch)
	sprintf(pch, "%lu", i);
    }
  }
  
  if(pars->verbose >= 4)
    for(uint64_t i = 0; i < pars->n_ind; i++)
      fprintf(stderr, "%s\n", pars->ind_labels[i]);

  // Read from GENO file
  if(pars->verbose >= 1)
    fprintf(stderr, "==> Reading genotype data\n");
  pars->in_geno_lkl = read_geno(pars->in_geno, pars->in_bin, pars->in_probs, pars->in_logscale, pars->n_ind, pars->n_sites);
  // Read_geno always returns genos in logscale
  pars->in_logscale = true;

  // Make copy of GLs in case of bootstrap
  pars->geno_lkl = init_ptr(pars->n_ind, pars->n_sites+1, 0, -INF);
  for(uint64_t i = 0; i < pars->n_ind; i++)
    memcpy(pars->geno_lkl[i], pars->in_geno_lkl[i], (pars->n_sites+1)*sizeof(double*));

  for(uint64_t i = 0; i < pars->n_ind; i++)
    for(uint64_t s = 1; s <= pars->n_sites; s++){
      // Call genotypes
      if(pars->call_geno)
	call_geno(pars->in_geno_lkl[i][s], N_GENO, pars->in_logscale, pars->N_thresh, pars->call_thresh, 0);

      // Convert space - from now on, all in NORMAL space!
      if(pars->in_logscale)
	conv_space(pars->in_geno_lkl[i][s], N_GENO, exp);
    }

  // Initialize random number generator
  if(pars->verbose >= 2)
    fprintf(stderr, "==> Setting seed for random number generator\n");
  pars->rnd_gen = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(pars->rnd_gen, pars->seed);

 

  //////////////////////
  // Open output file //
  //////////////////////
  FILE* out_fh = fopen(pars->out, "w");
  if(out_fh == NULL)
    error(__FUNCTION__, "cannot open output file!");
 


  /////////////////////
  // Prepare threads //
  /////////////////////
  // Create pthread structure array
  pth_struct* pth = new pth_struct[n_comb];
  // Initialize pars and distance matrix pointers
  uint64_t comb_id = 0;
  double** dist_matrix = init_ptr(pars->n_ind, pars->n_ind, 0.0);
  for(comb_id = 0; comb_id < n_comb; comb_id++){
    pth[comb_id].pars = pars;
    pth[comb_id].dist_matrix = dist_matrix;
  }

  // Create threadpool
  if( (pars->thread_pool = threadpool_create(pars->n_threads, n_comb, 0)) == NULL )
    error(__FUNCTION__, "failed to create thread pool!");



  //////////////////
  // Analyze Data //
  //////////////////
  fflush(stdout);
  // Loop for bootstrap analyses
  for(uint64_t rep = 0; rep <= pars->n_boot_rep; rep++){
    comb_id = 0;

    if(pars->verbose >= 1){
      if(rep == 0)
	// Full dataset analyses
	fprintf(stderr, "==> Analyzing full dataset...\n");
      else
	// Bootstrap analyses
	fprintf(stderr, "==> Bootstrap replicate # %lu ...\n", rep);
    }


    // Map from in_geno_lkl data to geno_lkl
    if(pars->verbose >= 2)
      fprintf(stderr, "> Mapping positions...\n");

    // If not rep 0, then adjust number of sites and random sample blocks from original dataset
    if(rep > 0){
      pars->n_sites -= pars->n_sites % pars->boot_block_size;
      rnd_map_data(pars, pars->n_sites/pars->boot_block_size);
    }

    // Calculate pairwise genetic distances
    if(pars->verbose >= 2)
      fprintf(stderr, "> Calculating pairwise genetic distances...\n");

    for(uint64_t i1 = 0; i1 < pars->n_ind; i1++)
      for(uint64_t i2 = i1+1; i2 < pars->n_ind; i2++){
	// Set which individuals to analyze
	pth[comb_id].i1 = i1;
	pth[comb_id].i2 = i2;

	// Add task to thread pool
	int ret = threadpool_add(pars->thread_pool, gen_dist_slave, (void*) &pth[comb_id++], 0);
	if(ret == -1)
          error(__FUNCTION__, "invalid thread pool!");
	else if(ret == -2)
          error(__FUNCTION__, "thread pool lock failure!");
	else if(ret == -3)
          error(__FUNCTION__, "queue full!");
	else if(ret == -4)
          error(__FUNCTION__, "thread pool is shutting down!");
	else if(ret == -5)
          error(__FUNCTION__, "thread failure!");
      }


  
    //////////////////////////
    // Wait for all threads //
    //////////////////////////
    threadpool_wait(pars->thread_pool, 0.1);

    if(n_comb != comb_id)
      error(__FUNCTION__, "some combinations are missing!");



    ///////////////////////////
    // Print Distance Matrix //
    ///////////////////////////
    if(pars->verbose >= 2)
      fprintf(stderr, "> Printing distance matrix\n");

    fprintf(out_fh, "\n%lu\n", pars->n_ind);
    for(uint64_t i = 0; i < pars->n_ind; i++){
      char* buf = join(dist_matrix[i], pars->n_ind, "\t");
      fprintf(out_fh, "%s\t%s\n", pars->ind_labels[i], buf);
      delete [] buf;
    }

  }

  threadpool_wait(pars->thread_pool, 0.1);
  if(threadpool_destroy(pars->thread_pool, threadpool_graceful) != 0)
    error(__FUNCTION__, "cannot free thread pool!");

  fclose(out_fh);



  /////////////////
  // Free Memory //
  /////////////////
  if(pars->verbose >= 1)
    fprintf(stderr, "==> Freeing memory...\n");

  free_ptr((void**) dist_matrix, pars->n_ind);
  delete [] pth;
  // pars struct
  free_ptr((void***) pars->in_geno_lkl, pars->n_ind, pars->n_sites+1);
  free_ptr((void**) pars->geno_lkl, pars->n_ind);
  free_ptr((void**) pars->ind_labels, pars->n_ind);
  //free_ptr((void*) pars->in_geno);
  gsl_rng_free(pars->rnd_gen);

  if(pars->verbose >= 1)
    fprintf(stderr, "Done!\n");

  delete pars;

  return 0;
}




double gen_dist(params *p, uint64_t i1, uint64_t i2){
  uint64_t cnt = 0;
  double dist = 0;

  Matrix<double> GL1 = alloc(1,N_GENO);
  Matrix<double> GL2 = alloc(1,N_GENO);
  int dim = GL1.y*GL2.y;

  for(uint64_t s = 1; s <= p->n_sites; s++){
    // Skip missing data
    if( p->pairwise_del && 
	(miss_data(p->geno_lkl[i1][s]) ||
	 miss_data(p->geno_lkl[i2][s])) )
      continue;

    double* sfs = init_ptr(dim, (double) 1/dim);
    GL1.mat[0][0] = p->geno_lkl[i1][s][0];
    GL1.mat[0][1] = p->geno_lkl[i1][s][1];
    GL1.mat[0][2] = p->geno_lkl[i1][s][2];
    GL2.mat[0][0] = p->geno_lkl[i2][s][0];
    GL2.mat[0][1] = p->geno_lkl[i2][s][1];
    GL2.mat[0][2] = p->geno_lkl[i2][s][2];

    if(!p->indep_geno)
      em2(sfs, &GL1, &GL2, 0.001, 50, dim);

    for(uint64_t g1 = 0; g1 < N_GENO; g1++)
      for(uint64_t g2 = 0; g2 < N_GENO; g2++){
	dist += p->score[g1][g2] * (p->indep_geno ? p->geno_lkl[i1][s][g1]*p->geno_lkl[i2][s][g2] : sfs[3*g1+g2]);

	if(p->verbose >= 9)
	  fprintf(stderr, "%lu\t%lu <-> %lu\t%lu - %lu\t%f\t%f\n", s, i1, i2, g1, g2, p->geno_lkl[i1][s][g1]*p->geno_lkl[i2][s][g2], sfs[3*g1+g2]);
      }

    if(p->verbose >= 8)
      fprintf(stderr, "Cumulative distance between indiv %lu and %lu at site %lu: %f\n", i1, i2, s, dist);

    cnt++;
    free_ptr(sfs);
  }

  if(p->verbose >=3)
    fprintf(stderr, "\tFound %lu sites valid between Ind %lu and Ind %lu!\n", cnt, i1, i2);

  dalloc(GL1, 1);
  dalloc(GL2, 1);

  // Calculates raw distance
  dist /= cnt;
  // Logarithmic transformation (log(1-d)) to make distance additive (assuming constant Ne) and avoid violating minimum evolution and NJ assumptions.
  dist = -log(1-dist);

  return dist;
}



void gen_dist_slave(void *pth){
  pth_struct* p = (pth_struct*) pth;

  p->dist_matrix[p->i1][p->i2] = p->dist_matrix[p->i2][p->i1] = gen_dist(p->pars, p->i1, p->i2);
}



void rnd_map_data(params *pars, uint64_t n_blocks){
  uint64_t block, rnd_block;
  uint64_t block_s, rnd_block_s;

  // For each of the original blocks
  for(block = 0; block < n_blocks; block++){
    // Pick a random one to replace it
    rnd_block = (uint64_t) floor( draw_rnd(pars->rnd_gen, 0, n_blocks) );

    // And copy its content
    for(uint64_t s = 1; s <= pars->boot_block_size; s++){
      block_s = block * pars->boot_block_size + s;
      rnd_block_s = rnd_block * pars->boot_block_size + s;

      if(pars->verbose >= 5)
	fprintf(stderr, "block: %lu\torig_site: %lu\trand_block:%lu\trand_site: %lu\n", block, block_s, rnd_block, rnd_block_s);

      for(uint64_t i = 0; i < pars->n_ind; i++)
	pars->geno_lkl[i][block_s] = pars->in_geno_lkl[i][rnd_block_s];
    }
  }
}
