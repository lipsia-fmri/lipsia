
#include <viaio/Vlib.h>

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>


#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/


extern double SIMPLS1(gsl_matrix *,gsl_vector *,gsl_vector *,size_t A);
extern void KNorm(gsl_matrix *X);


void MatFill(gsl_matrix_float *Xcorr,gsl_matrix *Xtrain,gsl_matrix *Xtest,
	     int *train,int *test,long *select)
{
  long i,j,k;
  size_t ntrain = Xtrain->size1;
  size_t dimX = Xtrain->size2;
  double u=0;

  for (i=0; i<dimX; i++) {
    k = select[i];
    for (j=0; j<ntrain; j++) {
      u = (double)gsl_matrix_float_get(Xcorr,train[j],k);
      gsl_matrix_set(Xtrain,j,i,u);
    }
    for (j=0; j<Xtest->size1; j++) {
      u = (double)gsl_matrix_float_get(Xcorr,test[j],k);
      gsl_matrix_set(Xtest,j,i,u);
    }
  }
}



double Loop(gsl_matrix_float *Xcorr,
	    gsl_vector *ytrain,gsl_vector *ytest,int *train,int *test,
	    gsl_matrix_long *S,size_t npls,gsl_vector *result)
{
  size_t iter;
  size_t niter = S->size1;
  size_t dimX  = S->size2;
  size_t ntest = result->size;
  size_t ntrain = ytrain->size;
  double kx=0;

  gsl_vector_set_zero(result);


#pragma omp parallel for schedule(guided)
  for (iter=0; iter<niter; iter++) {
    gsl_matrix *Xtrain = gsl_matrix_calloc(ntrain,dimX);
    gsl_matrix *Xtest  = gsl_matrix_calloc(ntest,dimX);
    gsl_vector *beta   = gsl_vector_calloc (dimX);
    gsl_vector *ytmp   = gsl_vector_calloc (ntest);
    long *pp = gsl_matrix_long_ptr(S,iter,0);
    
    MatFill(Xcorr,Xtrain,Xtest,train,test,pp);
    KNorm(Xtrain);
    KNorm(Xtest);

   
    SIMPLS1(Xtrain,ytrain,beta,npls);
    gsl_blas_dgemv(CblasNoTrans,1,Xtest,beta,0,ytmp);
 
#pragma omp critical
    {
      gsl_vector_add(result,ytmp);
      kx++;
    }
    gsl_matrix_free(Xtrain);
    gsl_matrix_free(Xtest);
    gsl_vector_free(beta);
    gsl_vector_free(ytmp);
  }
  gsl_vector_scale(result,1.0/kx);
  return 0;
}
