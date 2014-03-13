#include "ngsDist.hpp"



int read_geno(params* pars){
  uint64_t n_fields;
  // Depending on input we will have either 1 or 3 genot
  uint64_t n_geno = (pars->in_lkl ? N_GENO : 1);
  double* t;
  double* ptr;
  char* buf = new char[BUFF_LEN];

  // Allocate memory
  pars->geno_pp = init_double(pars->n_ind, pars->n_sites+1, N_GENO, -INFINITY);
  
  // Open GENO file
  gzFile in_geno_fh;
  if( (in_geno_fh = gzopen(pars->in_geno, pars->in_bin ? "rb" : "r")) == NULL )
    error("cannot open genotype file!");

  for(uint64_t s = 1; s <= pars->n_sites; s++){
    if(pars->in_bin){
      for(uint64_t i = 0; i < pars->n_ind; i++)
        if( gzread(in_geno_fh, pars->geno_pp[i][s], N_GENO * sizeof(double)) != N_GENO * sizeof(double) )
          error("cannot read GENO file!");
    }
    else{
      if( gzgets(in_geno_fh, buf, BUFF_LEN) == NULL)
        error("cannot read GENO file!");
      
      // Parse input line into array
      n_fields = split(buf, (const char*) " \t\r\n", &t);

      // Check if header and skip
      if(!n_fields){
        s--;
        continue;
      }

      if(n_fields < pars->n_ind * n_geno)
        error("wrong GENO file format!");
      
      // Use last "n_ind * n_geno" columns
      ptr = t + (n_fields - pars->n_ind * n_geno);
      
      if(pars->in_lkl)
        for(uint64_t i = 0; i < pars->n_ind; i++)
          for(uint64_t g = 0; g < N_GENO; g++)
            pars->geno_pp[i][s][g] = pars->in_log ? exp(ptr[i*N_GENO+g]) : ptr[i*N_GENO+g];
      else
        for(uint64_t i = 0; i < pars->n_ind; i++){
          int g = (int) ptr[i];
          pars->geno_pp[i][s][g] = log(1);
        }

      delete [] t;
    }
  }

  gzclose(in_geno_fh);
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
