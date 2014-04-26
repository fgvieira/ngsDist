
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
    printf("\tgeno file: %s\n\tlog-scale: %s\n\tlabels file: %s\n\tn_ind: %lu\n\tn_sites: %lu\n\tcall_geno: %s\n\tout prefix: %s\n\tthreads: %d\n\tversion: %s\n\tverbose: %d\n\tseed: %d\n\n",
	   pars->in_geno, pars->in_log ? "true":"false", pars->in_labels, pars->n_ind, pars->n_sites, pars->call_geno ? "true":"false", pars->out_prefix, pars->n_threads, pars->version ? "true":"false", pars->verbose, pars->seed);
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
      printf("==> BINARY input file (never log)\n");
    pars->in_bin = true;
    pars->in_log = false;
  }else
    error(__FUNCTION__, "invalid/corrupt genotype input file!");
  


  ////////////////////////////
  // Prepare initial values //
  ////////////////////////////
  // Read labels files
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
  pars->geno_lkl = read_geno(pars->in_geno, pars->in_bin, pars->in_lkl, pars->n_ind, pars->n_sites);
  
  for(uint64_t i = 0; i < pars->n_ind; i++)
    for(uint64_t s = 1; s <= pars->n_sites; s++){
      // Call genotypes
      if(pars->call_geno)
	call_geno(pars->geno_lkl[i][s], N_GENO, pars->in_log);
      // Convert space
      if(pars->in_log)
	conv_space(pars->geno_lkl[i][s], N_GENO, exp);
    }

  

  //////////////////
  // Analyze Data //
  //////////////////
  uint64_t comb_id = 0;
  double** dist_matrix = init_double(pars->n_ind, pars->n_ind, 0);
  // Create pthread structure
  pth_struct* pth = new pth_struct[n_comb];

  // Initialize semaphore
  if( sem_init(&pars->pth_sem, 0, pars->n_threads) )
    error(__FUNCTION__, "cannot initialise pthread sempahore!");

  if(pars->verbose >= 1)
    printf("==> Calculating pairwise genetic distances\n");

  for(uint64_t i1 = 0; i1 < pars->n_ind; i1++)
    for(uint64_t i2 = i1+1; i2 < pars->n_ind; i2++){
      sem_wait(&pars->pth_sem);
      // Initialize and set pthread detached attribute
      pthread_attr_init(&pth[comb_id].attr);
      pthread_attr_setdetachstate(&pth[comb_id].attr, PTHREAD_CREATE_DETACHED);
      // Initialize pthread structure
      pth[comb_id].pars = pars;
      pth[comb_id].dist_matrix = dist_matrix;
      pth[comb_id].i1 = i1;
      pth[comb_id].i2 = i2;
      // Launch pthread
      if( pthread_create(&pth[comb_id].id, &pth[comb_id].attr, gen_dist_slave, (void*) &pth[comb_id]) )
	error(__FUNCTION__, "cannot create thread!");
      comb_id++;

      if(pars->verbose >= 5)
	printf("> Launched thread for individuals %lu and %lu (comb# %lu).\n", i1, i2, comb_id);
    }

  if(n_comb != comb_id)
    error(__FUNCTION__, "missing combinations!");


  
  //////////////////////////
  // Wait for all threads //
  //////////////////////////
  int n_running_pthreads = 0;
  while(n_running_pthreads - pars->n_threads){
    sem_getvalue(&pars->pth_sem, &n_running_pthreads);
    sleep(1);
  }



  ///////////////////////////
  // Print Distance Matrix //
  ///////////////////////////
  if(pars->verbose >= 1)
    printf("==> Printing distance matrix\n");

  FILE* out_fh = fopen(strcat(pars->out_prefix, ".dist"), "w");
  if(out_fh == NULL)
    error(__FUNCTION__, "cannot open output file!");

  fprintf(out_fh, " %lu\n", pars->n_ind);
  for(uint64_t i = 0; i < pars->n_ind; i++){
    char* buf = join(dist_matrix[i], pars->n_ind, "\t");
    fprintf(out_fh, "%s\t%s\n", pars->ind_labels[i], buf);
    delete [] buf;
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
  free_ptr((void***) pars->geno_lkl, pars->n_ind, pars->n_sites+1);
  free_ptr((void**) pars->ind_labels, pars->n_ind);
  //free_ptr((void*) pars->in_geno);

  if(pars->verbose >= 1)
    printf("Done!\n");
  delete pars;

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
        //dist += p->geno_lkl[i1][s][g1] * p->geno_lkl[i2][s][g2] * p->score[g1][g2] * sfs[3*g1+g2];

    free_ptr(sfs);
  }

  dalloc(GL1, 1);
  dalloc(GL2, 1);
  return dist/p->n_sites;
}



void* gen_dist_slave(void* pth){
  pth_struct* p = (pth_struct*) pth;

  p->dist_matrix[p->i1][p->i2] = p->dist_matrix[p->i2][p->i1] = gen_dist(p->pars, p->i1, p->i2);

  // Free one slot for another thread
  if(sem_post(&p->pars->pth_sem))
    printf("WARN: semaphore post failed!\n");

  // Debug
  if(p->pars->verbose >= 6){
    int n_free_threads = 0;
    sem_getvalue(&p->pars->pth_sem, &n_free_threads);
    printf("Thread finished! Still running: %d\n", p->pars->n_threads - n_free_threads);
  }

  pthread_exit(NULL);
}
