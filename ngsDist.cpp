
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

#include <sys/stat.h>
#include "ngsDist.hpp"


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
    printf("\tgeno file: %s\n\tlog-scale: %s\n\tlabels file: %s\n\tn_ind: %lu\n\tn_sites: %lu\n\tcall_geno: %s\n\tout prefix: %s\n\tthreads: %d\n\tchunk size: %d\n\tversion: %s\n\tverbose: %d\n\tseed: %d\n\n",
	   pars->in_pp, pars->in_log ? "true":"false", pars->in_labels, pars->n_ind, pars->n_sites, pars->call_geno ? "true":"false", pars->out_prefix, pars->n_threads, pars->max_chunk_size, pars->version ? "true":"false", pars->verbose, pars->seed);
  }
  if(pars->verbose > 4) printf("==> Verbose values greater than 4 for debugging purpose only. Expect large amounts of info on screen\n");



  /////////////////////
  // Check Arguments //
  /////////////////////
  if(pars->in_pp == NULL)
    error("Genotype input file (-geno) missing!");
  if(pars->n_ind == 0)
    error("Number of individuals (-n_ind) missing!");
  if(pars->n_sites == 0)
    error("Number of sites (-n_sites) missing!");
  
  
  
  ///////////////////////
  // Adjust Parameters //
  ///////////////////////
  // Adjust max_chunk_size in case of fewer sites
  if(pars->max_chunk_size > pars->n_sites)
    pars->max_chunk_size = pars->n_sites;

  // Calculate total number of chunks
  pars->n_chunks = ceil( (double) pars->n_sites / (double) pars->max_chunk_size );
  if( pars->verbose >= 1 ) printf("==> Analysis will be run in %d chunk(s)\n", pars->n_chunks);

  // Adjust thread number to chunks
  if(pars->n_chunks < pars->n_threads)
    pars->n_threads = pars->n_chunks;



  ///////////////////////
  // Check input files //
  ///////////////////////
  // Get file total size
  struct stat st;
  if( stat(pars->in_pp, &st) != 0 )
    error("cannot check file size!");
  if( strcmp(strrchr(pars->in_pp, '.'), ".gz") == 0 ){
    printf("==> GZIP input file (never BINARY)\n");
    pars->in_bin = false;
  }else if( pars->n_sites == st.st_size/sizeof(double)/pars->n_ind/N_GENO ){
    printf("==> BINARY input file (never log)\n");
    pars->in_bin = true;
    pars->in_log = false;
  }else
    error("invalid/corrupt genotype input file!");
  


  ////////////////////////////
  // Prepare initial values //
  ////////////////////////////
  // Read labels files
  pars->ind_labels = init_char(pars->n_ind, BUFF_LEN, (const char*) "Ind_#");
  if(pars->in_labels){
    if(pars->verbose >= 1)
      printf("==> Reading labels\n");
    read_labels(pars);
  }

  // Read from GENO file
  if(pars->verbose >= 1)
    printf("==> Reading genotype posterior probabilities\n");
  read_geno(pars);
  
  // Call genotypes
  if(pars->call_geno)
    for(uint64_t i = 0; i < pars->n_ind; i++)
      for(uint64_t s = 1; s <= pars->n_sites; s++)
	call_geno(pars->post_prob[i][s], N_GENO);

  

  //////////////////
  // Analyze Data //
  //////////////////
  double dist_matrix[pars->n_ind][pars->n_ind];

  if(pars->verbose >= 1)
    printf("==> Calculating pairwise genetic distances\n");

  for(uint64_t i1 = 0; i1 < pars->n_ind; i1++)
    for(uint64_t i2 = 0; i2 < pars->n_ind; i2++)
      if(i1 == i2)
	dist_matrix[i1][i2] = 0;
      else
	dist_matrix[i1][i2] = gen_dist(pars, i1, i2);



  ///////////////////////////
  // Print Distance Matrix //
  ///////////////////////////
  if(pars->verbose >= 1)
    printf("==> Printing distance matrix\n");

  FILE* out_fh = fopen(strcat(pars->out_prefix, ".dist"), "w");
  if(out_fh == NULL)
    error("cannot open output file!");

  fprintf(out_fh, " %lu\n", pars->n_ind);
  for(uint64_t i = 0; i < pars->n_ind; i++){
    char* buf = merge(dist_matrix[i], pars->n_ind, "\t");
    fprintf(out_fh, "%s\t%s\n", pars->ind_labels[i], buf);
    delete [] buf;
  }

  fclose(out_fh);



  /////////////////
  // Free Memory //
  /////////////////
  if(pars->verbose >= 1)
    printf("Freeing memory...\n");
  // pars struct
  free_ptr((void***) pars->post_prob, pars->n_ind, pars->n_sites+1);
  free_ptr((void**) pars->ind_labels, pars->n_ind);
  //free_ptr((void*) pars->in_pp);

  if(pars->verbose >= 1)
    printf("Done!\n");
  delete pars;

  return 0;
}




double gen_dist(params* p, uint64_t i1, uint64_t i2){
  double dist = 0;

  for(uint64_t s = 1; s <= p->n_sites; s++)
    for(uint64_t g1 = 0; g1 < N_GENO; g1++)
      for(uint64_t g2 = 0; g2 < N_GENO; g2++)
	dist += p->post_prob[i1][s][g1] * p->post_prob[i2][s][g2] * p->score[g1][g2];
	  
  return dist/p->n_sites;
}
