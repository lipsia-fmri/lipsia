
/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/Vlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>

/* From the standard C library: */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../cylutils/cyl.h"

/* 2nd derivatives of chebychev polynomials */
#define D2(x) (4.0)
#define D3(x) (24.0*(x))
#define D4(x) ((96.0*(x)*(x) - 16.0))

double TDeriv2_order4(double x,gsl_vector *beta)
{
  double z = beta->data[2]*D2(x) + beta->data[3]*D3(x) + beta->data[4]*D4(x);
  return z;
}


/*  2nd deriv, Chebychev order 3 */
double TDeriv2_order3(double x,gsl_vector *beta)
{
  double z = 4.0*beta->data[2] + 24.0*beta->data[3]*x; 
  return z;
}

double TFit3(double x,gsl_vector *beta)
{
  double s = beta->data[0] + beta->data[1]*P1(x) + beta->data[2]*P2(x) + beta->data[3]*P3(x);
  return s;
}


/*  2nd deriv, standard polynomial, order 4 */
double SDeriv2(double x,gsl_vector *beta)
{
  double z = 2.0*beta->data[2] + 6.0*beta->data[3]*x + 12.0*beta->data[4]*x*x;
  return z;
}

double SFit3(double x,gsl_vector *beta)
{
  double s = beta->data[0] + beta->data[1]*x + beta->data[2]*S2(x) + beta->data[3]*S3(x);
  return s;
}



void YNorm(gsl_vector *x)
{
  double mean = gsl_stats_mean(x->data,1,x->size);
  gsl_vector_add_constant(x,-mean);
}


/* fraction of points in [0,1] for which f''(x) > 0 */
double Fraction(gsl_vector *beta)
{
  double a = beta->data[2];
  double b = beta->data[3];

  if (fabs(b) < TINY) return 0;
  if (fabs(b) < TINY && a > 0) return 1;
  if (fabs(b) < TINY && a <= 0) return 0;

  /* zero crossing */
  /* double x = -a/(3.0*b);  */  /* standard polynomial */
  double x = -a/(6.0*b);  /* chebychev polynomial */
  
  double z = 0;
  if (b > 0) z = 1.0 - x;
  if (b < 0) z = x;
  if (z < 0) z = 0;
  if (z > 1) z = 1;
  return z;
}


double LaminarConvex(gsl_vector *y0,gsl_vector *mvec0,gsl_vector *result,int verbose)
{
  size_t i,n,p=4;
  double u=0,v=0,z=0,r2=-1;
  double xmin=0.2,xmax=0.8;
  gsl_matrix *X = NULL;
  gsl_multifit_linear_workspace *work = NULL;

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
  work = gsl_multifit_linear_alloc(n,p);
  X = gsl_matrix_calloc(n,p);
  
  for (i=0; i<n; i++) {
    u = mvec->data[i];
    gsl_matrix_set(X,i,0,1);
    gsl_matrix_set(X,i,1,P1(u));
    gsl_matrix_set(X,i,2,P2(u));
    gsl_matrix_set(X,i,3,P3(u));
  }

  /* GLM */
  gsl_multifit_linear(X,y,beta,cov,&chisq,work);
  r2 = 1.0 - chisq/tss;

  /* normalize */
  /* if 2nd deriv > 0: concave down */
  /* if 2nd deriv < 0: concave up */
  /* if 1st deriv = 0: stationary point (min,max,inflection) */
  /* 1st deriv==0, and 2nd deriv > 0: local min */
  /* 1st deriv==0, and 2nd deriv < 0: local max */

  /* fraction of points in [0,1] for which f''(x) > 0 */
  u = Fraction(beta);
  v = 1-u;
  z = u;
  if (v > u) z = -v;
  result->data[0] = z;

  /*
  if (verbose) {
    FILE *fp = fopen("l.txt","w");
    for (i=0; i<n; i++) {
      u = mvec->data[i];
      fprintf(fp," %f %f\n",u,y->data[i]);
    }
    fprintf(fp,"\n\n\n");

    for (u=0; u<=1; u+=0.01) {
      fprintf(fp," %f %f\n", u,TFit3(u,beta));
    }
    fclose(fp);
    exit(0);
  }
  */
  
  /* free mem */
  gsl_vector_free(mvec);
  gsl_vector_free(y);
  gsl_vector_free(beta);
  gsl_multifit_linear_free(work);
  gsl_matrix_free(cov);
  gsl_matrix_free(X);

  return r2;
}
