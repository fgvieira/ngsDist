#include "ngsDist.hpp"



int read_geno(params* pars){
  char* buf = new char[BUFF_LEN];

  // Allocate memory
  pars->post_prob = init_double(pars->n_ind, pars->n_sites+1, N_GENO, -INFINITY);
  
  // Open GENO file
  gzFile in_geno_fh;
  if( (in_geno_fh = gzopen(pars->in_geno, pars->in_bin ? "rb" : "r")) == NULL )
    error("cannot open genotype file!");

  for(uint64_t s = 1; s <= pars->n_sites; s++){
    if(pars->in_bin){
      for(uint64_t i = 0; i < pars->n_ind; i++){
	if( gzread(in_geno_fh, pars->post_prob[i][s], N_GENO * sizeof(double)) != N_GENO * sizeof(double) )
	  error("cannot read GENO file!");
      }
    }
    else{
      if( gzgets(in_geno_fh, buf, BUFF_LEN) == NULL)
	error("cannot read GENO file!");
      
      double* t = NULL;
      if( split(buf, (const char*) " \t\r\n", &t) != pars->n_ind * N_GENO )
	error("wrong GENO file format!");
      
      for(uint64_t i = 0; i < pars->n_ind; i++)
	for(uint64_t g = 0; g < N_GENO; g++)
	  pars->post_prob[i][s][g] = pars->in_log ? exp(t[i*N_GENO+g]) : t[i*N_GENO+g];
      
      delete [] t;
    }
  }
  
  gzclose(in_geno_fh);
  delete [] buf;
  return 0;
}




uint64_t read_chunk(double** chunk_data, params* pars, uint64_t chunk) {
  uint64_t total_elems_read = 0;
  
  if(chunk >= pars->n_chunks)
    error("invalid chunk number!");
  
  // Define chunk start and end positions
  uint64_t start_pos = chunk * pars->max_chunk_size;
  uint64_t end_pos = start_pos + pars->max_chunk_size - 1;
  if(end_pos >= pars->n_sites)	end_pos = pars->n_sites - 1;
  uint64_t chunk_size = end_pos - start_pos + 1;
  if( pars->verbose >= 6 ) printf("\tReading chunk %lu from position %lu to %lu (%lu)\n", chunk+1, start_pos, end_pos, chunk_size);
  
  // Read data from file
  for(uint64_t c = 0; c < chunk_size; c++) {
    //    chunk_data[c] = pars->post_prob[start_pos+c];
    uint64_t elems_read = pars->n_ind * 3;

    if( elems_read != pars->n_ind * 3 )
      error("cannot read GLF file!");
    total_elems_read += elems_read;
  }

  return( total_elems_read/(pars->n_ind * 3) );
}
