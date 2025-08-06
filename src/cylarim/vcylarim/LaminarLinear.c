
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
#include <gsl/gsl_chebyshev.h>

/* From the standard C library: */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../cylutils/cyl.h"


double LaminarLinear(gsl_vector *y,gsl_vector *mvec,gsl_vector *beta)
{
  size_t i,p=2;
  size_t n=y->size;
  double u=0,r2=-1;
  gsl_matrix *X = NULL;
  gsl_multifit_linear_workspace *work = NULL;

  gsl_vector_set_zero(beta);
  double tss = gsl_stats_tss(y->data,1,y->size);
  if (tss < TINY) return -1;

  /* set up GLM */
  double chisq=0;
  gsl_vector *zbeta = gsl_vector_calloc(p);
  gsl_matrix *cov = gsl_matrix_calloc(p,p);
  work = gsl_multifit_linear_alloc(n,p);
  X = gsl_matrix_calloc(n,p);
  
  for (i=0; i<n; i++) {
    u = mvec->data[i];
    gsl_matrix_set(X,i,0,1);
    gsl_matrix_set(X,i,1,u);
  }

  /* GLM */
  gsl_multifit_linear(X,y,zbeta,cov,&chisq,work);
  r2 = 1.0 - chisq/tss;

  beta->data[0] = zbeta->data[1];

  gsl_vector_free(zbeta);
  gsl_multifit_linear_free(work);
  gsl_matrix_free(cov);
  gsl_matrix_free(X);

  return r2;
}
