
/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/Vlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

/* From the standard C library: */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../cylutils/cyl.h"

double LaminarFit(gsl_vector *y,gsl_vector *mvec,gsl_vector *zval)
{
  size_t i,p;
  size_t n = y->size;
  double u=0,r2=-1;
  double chisq=0;

  gsl_vector_set_zero(zval);
  double tss = gsl_stats_tss(y->data,1,y->size);
  if (tss < TINY) return -1;


  /* set up GLM */
  for (p=1; p<7; p++) {
    size_t pa=p+1;
    gsl_vector *beta = gsl_vector_calloc(pa);
    gsl_matrix *cov = gsl_matrix_calloc(pa,pa);
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n,pa);
    gsl_matrix *X = gsl_matrix_calloc(n,pa);

    for (i=0; i<n; i++) {
      u = mvec->data[i];
      u = 2.0*u - 1.0;    
      gsl_matrix_set(X,i,0,1.0);
      if (p > 0) gsl_matrix_set(X,i,1,P1(u));
      if (p > 1) gsl_matrix_set(X,i,2,P2(u));
      if (p > 2) gsl_matrix_set(X,i,3,P3(u));
      if (p > 3) gsl_matrix_set(X,i,4,P4(u));
      if (p > 4) gsl_matrix_set(X,i,5,P5(u));
    }
  
    /* GLM */
    gsl_multifit_linear(X,y,beta,cov,&chisq,work);


    /* goodness of fit, R^2 */
    r2 = 1.0 - chisq/tss;
    zval->data[p-1] = r2;
    
    gsl_vector_free(beta);
    gsl_multifit_linear_free(work);
    gsl_matrix_free(cov);
    gsl_matrix_free(X);
  }
  return r2;
}
