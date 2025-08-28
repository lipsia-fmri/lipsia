
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

extern double Xt2z(double t,double df);
extern void YNorm(gsl_vector *x);


double GetMax(gsl_vector *mvec,gsl_vector *z)
{
  size_t i;
  double zmax=-99999,umax=0;
  for (i=0; i<z->size; i++) {
    if (z->data[i] > zmax && mvec->data[i] > 0.1 && mvec->data[i] < 0.9) {
      zmax = z->data[i];
      umax = mvec->data[i];
    }
  }
  return umax;
}

double GetMin(gsl_vector *mvec,gsl_vector *z)
{
  size_t i;
  double zmin=99999,umin=0;
  for (i=0; i<z->size; i++) {
    if (z->data[i] < zmin && mvec->data[i] > 0.1 && mvec->data[i] < 0.9) {
      zmin = z->data[i];
      umin = mvec->data[i];
    }
  }
  return umin;
}

void PrintPeak(VImage zmap,VImage metric,VImage dest,Cylinders *cyl,size_t cid,char *filename)
{
  int b,r,c;
  size_t n0 = cyl->addr[cid]->size;
  size_t i,j,k;
  size_t p = 5;
  double u=0,r2=0;
  gsl_matrix *X = NULL;
  gsl_multifit_linear_workspace *work = NULL;

  gsl_vector *y0 = gsl_vector_calloc(n0);
  gsl_vector *mvec0 = gsl_vector_calloc(n0);
  for (i=0; i<n0; i++) {
    k = cyl->addr[cid]->data[i];
    b = gsl_matrix_int_get(cyl->xmap,k,0);
    r = gsl_matrix_int_get(cyl->xmap,k,1);
    c = gsl_matrix_int_get(cyl->xmap,k,2);
    mvec0->data[i] = VPixel(metric,b,r,c,VFloat);
    y0->data[i] = VPixel(zmap,b,r,c,VFloat);

    if (mvec0->data[i] > 0.05 && mvec0->data[i] < 0.95)
      VPixel(dest,b,r,c,VFloat) = y0->data[i];
  }
 
  
  size_t n = 0;
  for (i=0; i<mvec0->size; i++) {
    if (mvec0->data[i] > 0.05 && mvec0->data[i] < 0.95) n++;
  }
  gsl_vector *y = gsl_vector_calloc(n);
  gsl_vector *mvec = gsl_vector_calloc(n);

  j=0;
  for (i=0; i<mvec0->size; i++) {
    if (mvec0->data[i] > 0.05 && mvec0->data[i] < 0.95) {
      y->data[j] = y0->data[i];
      mvec->data[j] = mvec0->data[i];
      j++;
    }
  }

  /* set up GLM */
  double chisq=0;
  gsl_vector *beta = gsl_vector_calloc(p);
  gsl_matrix *cov = gsl_matrix_calloc(p,p);
  work = gsl_multifit_linear_alloc(n,p);
  X = gsl_matrix_calloc(n,p);
  
  for (i=0; i<n; i++) {
    u = mvec->data[i];
    u = 2.0*u - 1.0;    
    gsl_matrix_set(X,i,0,1);
    gsl_matrix_set(X,i,1,T1(u));
    gsl_matrix_set(X,i,2,T2(u));
    gsl_matrix_set(X,i,3,T3(u));
    gsl_matrix_set(X,i,4,T4(u));
  }

  /* GLM */
  gsl_multifit_linear(X,y,beta,cov,&chisq,work);

  gsl_vector *z = gsl_vector_calloc(y->size);
  gsl_blas_dgemv(CblasNoTrans,1,X,beta,0,z);
  
  /* goodness of fit, R^2 */
  double tss = gsl_stats_tss(y->data,1,y->size);
  r2 = 1.0 - chisq/tss;
  
  double xmin = GetMin(mvec,z);
  double xmax = GetMax(mvec,z);

  fprintf(stderr," fit, R2= %.3f\n",r2);
  fprintf(stderr,"  trough: %.3f,   peak %.3f\n",xmin,xmax);
  FILE *fp = fopen(filename,"w");
  for (i=0; i<n; i++) {
    fprintf(fp," %f %f %f\n",mvec->data[i],y->data[i],z->data[i]);
  }
  fclose(fp);
  
  gsl_vector_free(y);
  gsl_vector_free(mvec);
  gsl_multifit_linear_free(work);
  gsl_matrix_free(cov);
  gsl_matrix_free(X);
}
