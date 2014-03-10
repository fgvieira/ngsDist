#include "ngsDist.hpp"



int read_geno(params* pars){
  char* buf = new char[BUFF_LEN];

  // Allocate memory
  pars->post_prob = init_double(pars->n_ind, pars->n_sites+1, N_GENO, 0);
  
  // Open GENO file
  gzFile in_pp_fh = gzopen(pars->in_pp, pars->in_bin ? "rb" : "r");
  if(in_pp_fh == NULL)
    error("cannot open genotype file!");

  for(uint64_t s = 1; s <= pars->n_sites; s++){
    if(pars->in_bin){
      for(uint64_t i = 0; i < pars->n_ind; i++){
	if( gzread(in_pp_fh, pars->post_prob[i][s], N_GENO * sizeof(double)) != N_GENO * sizeof(double) )
	  error("cannot read GENO file!");
      }
    }
    else{
      if( gzgets(in_pp_fh, buf, BUFF_LEN) == NULL)
	error("cannot read GENO file!");
      
      double* t = NULL;
      double* ptr;
      uint64_t n_fields = split(buf, (const char*) " \t\r\n", &t);

      if(n_fields == pars->n_ind * N_GENO)
	ptr = t;
      else if(n_fields == pars->n_ind * N_GENO + 2)
	// If ANGSD format (chr, pos, ...), skip first 2 columns
	ptr = t + 2;
      else if(n_fields == pars->n_ind * N_GENO + 3)
	// If BEAGLE format (id, allele1, allele2, ...), skip first 3 columns
	ptr = t + 3;
      else
	error("wrong GENO file format!");

      for(uint64_t i = 0; i < pars->n_ind; i++)
	for(uint64_t g = 0; g < N_GENO; g++)
	  pars->post_prob[i][s][g] = pars->in_log ? exp(ptr[i*N_GENO+g]) : ptr[i*N_GENO+g];
      
      delete [] t;
    }
  }
  
  gzclose(in_pp_fh);
  delete [] buf;
  return 0;
}




int read_labels(params* pars){
  FILE* in_labels_fh = fopen(pars->in_labels, "r");
  if(in_labels_fh == NULL)
    error("cannot open labels file!");

  for(uint64_t i = 0; i < pars->n_ind; i++){
    fgets(pars->ind_labels[i], BUFF_LEN, in_labels_fh);
    // Remove trailing newline
    pars->ind_labels[i][strlen(pars->ind_labels[i])-1] = '\0';
  }

  fclose(in_labels_fh);
  return 0;
}
