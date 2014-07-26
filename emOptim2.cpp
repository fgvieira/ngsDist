/*
  The functionality of this file, has replaced the old emOptim and testfolded.c programs.

  part of ANGSD

  GNU license or whetever its called

  thorfinn@binf.ku.dk

  fixme: minor leaks in structures related to the thread structs, and the append function.
  
  Its july 13 2013, it is hot outside

 */


template <typename T>
struct Matrix{
  size_t x;
  size_t y;
  T** mat;
};



size_t fsize(const char* fname){
  struct stat st ;
  stat(fname,&st);
  return st.st_size;
}


Matrix<double> alloc(size_t x,size_t y){
  //fprintf(stderr,"def=%f\n",def);
  Matrix<double> ret;
  ret.x=x;
  ret.y=y;
  ret.mat= new double*[x];
  for(size_t i=0;i<x;i++)
    ret.mat[i]=new double[y];
  return ret;
}
void dalloc(Matrix<double> &ret,size_t x){
  for(size_t i=0;i<x;i++)
    delete [] ret.mat[i];
  delete [] ret.mat;
}



typedef struct emPars_t{
  int threadId; //size_t is the largest primitive datatype.
  double *inparameters;
  double *outparameters;
  Matrix<double> *GL1;
  Matrix<double> *GL2;
  int from;
  int to;
  int dim;
  double lik;
  double *sfs;//shared for all threads
  double *post;//allocated for every thread
  double *grad;//used for bfgs, length is dim-1; 
} emPars;

emPars *emp = NULL;


void normalize(double *tmp,int len){
  double s=0;
  for(int i=0;i<len;i++)
    s += tmp[i];
  for(int i=0;i<len;i++)
    tmp[i] /=s;
}

double lik2(double *sfs,Matrix<double> *GL1,Matrix<double> *GL2,size_t start,size_t stop){
  double res =0;
  for(size_t s=start; s<stop; s++){
    //    fprintf(stderr,"s=%d\n",s);
    double tmp =0;
    int inc =0;
    for(size_t x=0;x<GL1->y;x++)
      for(size_t y=0;y<GL2->y;y++)
	tmp += sfs[inc++]* GL1->mat[s][x]*GL2->mat[s][y];
    res +=log(tmp);
  }
  return res;
}

void emStep2(double *pre,Matrix<double> *GL1,Matrix<double> *GL2,double *post,int start,int stop,int dim){
  double inner[dim];
  for(int x=0;x<dim;x++)
    post[x] =0.0;
    
  for(int s=start; s<stop; s++){
    int inc=0;
    for(size_t x=0;x<GL1->y;x++)
      for(size_t y=0;y<GL2->y;y++){
	inner[inc] = pre[inc]*GL1->mat[s][x]*GL2->mat[s][y];
	inc++;
      }
   normalize(inner,dim);
   for(int x=0;x<dim;x++)
     post[x] += inner[x];
  }
  normalize(post,dim);
 
}


void em2(double *sfs,Matrix<double> *GL1,Matrix<double> *GL2,double tole=0.01,int maxIter=10,int dim=0){
  double oldLik,lik;
  oldLik = lik2(sfs,GL1,GL2,0,GL1->x);
  //printf("startlik=%f\n",oldLik);

  double *tmp = new double[dim];//<- wont be cleaned up, but is only allocated once
  for(int it=0; it<maxIter; it++) {
    emStep2(sfs,GL1,GL2,tmp,0,GL1->x,dim);

    for(int i=0;i<dim;i++)
      sfs[i]= tmp[i];
    
    lik = lik2(sfs,GL1,GL2,0,GL1->x);
    
    //printf("[%d] lik=%f diff=%g\n",it,lik,fabs(lik-oldLik));

    if(fabs(lik-oldLik)<tole){
      oldLik=lik;
      break;
    }
    oldLik=lik;
  }
  delete [] tmp;
}


void print(Matrix<double> &mat,FILE *fp){
  for(size_t x=0;x<mat.x;x++){
    //    fprintf(stderr,"x=%d\n",x);
    for(size_t y=0;y<mat.y;y++)
      fprintf(fp,"%f ",mat.mat[x][y]);
    fprintf(fp,"\n");
  }
}
