
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


double TFit4(double x,gsl_vector *beta)
{
  double s = beta->data[0] + beta->data[1]*P1(x) + beta->data[2]*P2(x)
    + beta->data[3]*P3(x) + beta->data[4]*P4(x);
  return s;
}

double SFit4(double x,gsl_vector *beta)
{
  double s = beta->data[0] + beta->data[1]*S1(x) + beta->data[2]*S2(x)
    + beta->data[3]*S3(x) + beta->data[4]*S4(x);
  return s;
}


double GetMax(gsl_vector *mvec,gsl_vector *z)
{
  size_t i;
  double zmax=-99999,umax=0;
  for (i=0; i<z->size; i++) {
    if (z->data[i] > zmax) {
      zmax = z->data[i];
      umax = mvec->data[i];
    }
  }
  if (umax < 0 || umax > 1) umax = 0;
  return umax;
}

double GetMin(gsl_vector *mvec,gsl_vector *z)
{
  size_t i;
  double zmin=99999,umin=0;
  for (i=0; i<z->size; i++) {
    if (z->data[i] < zmin) {
      zmin = z->data[i];
      umin = mvec->data[i];
    }
  }
  if (umin < 0 || umin > 1) umin = 0;
  return umin;
}

double LaminarPeak(gsl_vector *y0,gsl_vector *mvec0,gsl_vector *result,int verbose)
{
  size_t i,n,p=5;
  double u=0,r2=-1;
  double xmin=0.1,xmax=0.9;

  gsl_vector_set_zero(result);

  n=0;
  for (i=0; i<y0->size; i++) {
    if (mvec0->data[i] < xmin || mvec0->data[i] > xmax) continue;
    n++;
  }
  if (n < 20) return 0;

  gsl_vector *y = gsl_vector_calloc(n);
  gsl_vector *mvec = gsl_vector_calloc(n);
  n=0;
  for (i=0; i<y0->size; i++) {
    if (mvec0->data[i] < xmin || mvec0->data[i] > xmax) continue;
    mvec->data[n] = mvec0->data[i];
    y->data[n] = y0->data[i];
    n++;
  }
  double tss = gsl_stats_tss(y->data,1,y->size);
  if (tss < TINY) return -1;

  
  /* set up GLM */
  double chisq=0;
  gsl_vector *beta = gsl_vector_calloc(p);
  gsl_matrix *cov = gsl_matrix_calloc(p,p);
  gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n,p);
  gsl_matrix *X = gsl_matrix_calloc(n,p);
  
  for (i=0; i<n; i++) {
    u = mvec->data[i];
    gsl_matrix_set(X,i,0,1);
    gsl_matrix_set(X,i,1,P1(u));
    gsl_matrix_set(X,i,2,P2(u));
    gsl_matrix_set(X,i,3,P3(u));
    gsl_matrix_set(X,i,4,P4(u));
  }

  /* GLM */
  gsl_multifit_linear(X,y,beta,cov,&chisq,work);
  r2 = 1.0 - chisq/tss;

  /* get trough,peak */
  gsl_vector *z = gsl_vector_calloc(y->size);
  gsl_blas_dgemv(CblasNoTrans,1,X,beta,0,z);
  result->data[0] = GetMin(mvec,z);
  result->data[1] = GetMax(mvec,z);

  /*
  if (verbose) {
    fprintf(stderr,"\n valley: %f\n",result->data[0]);
    FILE *fp = fopen("l.txt","w");
    for (i=0; i<n; i++) {
      u = mvec->data[i];
      fprintf(fp," %f %f\n",u,y->data[i]);
    }
    fprintf(fp,"\n\n\n");

    for (u=0; u<=1; u+=0.01) {
      fprintf(fp," %f %f\n", u,TFit4(u,beta));
    }
    fclose(fp);
    exit(0);
  }
  */
  
  gsl_vector_free(mvec);
  gsl_vector_free(y);
  gsl_vector_free(z);
  gsl_vector_free(beta);
  gsl_multifit_linear_free(work);
  gsl_matrix_free(cov);
  gsl_matrix_free(X);
  return r2;
}
