/*
**
** G.Lohmann, Apr 2020
*/
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>


extern void RNorm(gsl_vector *vec);
extern void LNorm(gsl_vector *vec);
extern void ZNorm(gsl_vector *vec);
extern void QNorm(gsl_vector *vec);
extern void MNorm(gsl_vector *vec);

extern gsl_vector *ReadRegressor(VString filename);

extern void printmat(gsl_matrix *A);
extern void YCheckImage(VImage src);


/* flag=1: training data,  flag=0: test data */
void Confounds(gsl_matrix *C,gsl_vector *y,gsl_vector *bx,double offset,int *table,int flag)
{
  size_t i,j,k,n=y->size;
  double u,chisq=0;

  
  /* fill matrix of confounds */
  gsl_matrix *X = gsl_matrix_calloc(n,C->size2);
  gsl_matrix_set_all(X,1.0);
  for (i=0; i<X->size1; i++) {
    k = table[i];
    for (j=0; j<X->size2; j++) {
      u = gsl_matrix_get(C,k,j);
      gsl_matrix_set(X,i,j,u);
    }
  }

  /* OLS fit only for training data, not for test data */
  if (flag > 0) {
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(X->size1,X->size2);
    gsl_matrix *cov = gsl_matrix_alloc (X->size2,X->size2);
    gsl_multifit_linear(X,y,bx,cov,&chisq,work);
    gsl_matrix_free(cov);
    gsl_multifit_linear_free(work);
  }

  /* apply centering for test data using offset from training data */
  if (flag == 0) {
    for (i=0; i<n; i++) y->data[i] -= offset;
  }

  /* get residuals */
  gsl_vector *residual = gsl_vector_calloc(n);
  gsl_multifit_linear_residuals(X,y,bx,residual);
  for (i=0; i<n; i++) y->data[i] = residual->data[i];

  /* reverse centering */
  if (flag == 0) {
    for (i=0; i<n; i++) y->data[i] += offset;
  }
  
  gsl_vector_free(residual);
  gsl_matrix_free(X);
}



void XConfound(gsl_matrix *C,gsl_vector *y)
{
  int i;
  int *table = (int *) VCalloc(y->size,sizeof(int));
  for (i=0; i<y->size; i++) table[i] = i;
  gsl_vector *bx = gsl_vector_calloc(C->size2);
  
  double mean0 = gsl_stats_mean(y->data,1,y->size);
  gsl_vector_add_constant(y,-mean0);
  
  Confounds(C,y,bx,mean0,table,1);

  gsl_vector_add_constant(y,mean0);
  
  VFree(table);
  gsl_vector_free(bx);
}
