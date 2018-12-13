/*
** Test for irreducibility of a symmetric matrix R = X^T X,
** i.e. there is at least one positive off-diagonal entry in each row and column.
** Since R is symmetric, it suffices to show that each row has one non-negative entry.
** i.e. (X^TX v - v)_i > 0 for all i=0,...,n-1 
** exclude main diagonal
**
** G.Lohmann, MPI-KYB, Nov 2018
*/

#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern void SignFlip(gsl_matrix_float *X);
  
/* 
** if irreducible return 1, else return 0 
*/
int Irreducible(gsl_matrix_float *X,int type)
{
  size_t i;
  size_t nvox = X->size1;
  size_t nt = X->size2;
  float u=0,tiny=1.0e-6;

  /* ini */
  gsl_vector_float *x = gsl_vector_float_calloc(nvox);
  gsl_vector_float *y = gsl_vector_float_calloc(nvox);
  gsl_vector_float *z = gsl_vector_float_calloc(nt);
  for (i=0; i<nvox; i++) x->data[i] = y->data[i] = 1.0;

  
  /* R = X^t X */
  gsl_blas_sgemv (CblasTrans,1.0,X,x,0.0,z);
  if (type == 7) SignFlip(X);
  gsl_blas_sgemv (CblasNoTrans,1.0,X,z,0.0,y);
  if (type == 7) SignFlip(X);
  
  
  /* count off-diagonal entries with value '0' */
  size_t n=0;
  for (i=0; i<nvox; i++) {
    u = y->data[i] - x->data[i];  /* subtract main diagonal */
    if (fabs(u) < tiny) n++;
  }
  gsl_vector_float_free(x);
  gsl_vector_float_free(y);
  gsl_vector_float_free(z);

  
  /* return code */
  int rtcode = 0;
  if (n > 0) rtcode = 0;
  else rtcode = 1;
  if (rtcode == 0) {
    VWarning(" Matrix is not irreducible,   n=%lu,  nvox= %lu\n",n,nvox);
  }
  return rtcode;
}
