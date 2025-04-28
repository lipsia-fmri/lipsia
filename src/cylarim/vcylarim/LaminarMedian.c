
/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/Vlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>


/* From the standard C library: */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "../cylutils/cyl.h"

extern double Xt2z(double t,double df);
extern void AssignLabels(gsl_vector *y,int *label,double *y0,double *y1,double *y2,size_t n);
extern int LayerLabels(gsl_vector *mvec,int *labels,size_t *n0,size_t *n1,size_t *n2);

extern double kth_smallest(double *a, size_t n, size_t k);
#define Median(a,n) kth_smallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))


/* median  */
double XMedian(double *data,double *tmp,size_t n)
{
  size_t i;
  double median=0;
  for (i=0; i<n; i++) tmp[i] = data[i];
  if (n > 10) median = Median(tmp,n);
  else median = gsl_stats_median(tmp,1,n);
  return median;
}

double MedianTest(double *x1,double *x2,size_t n1,size_t n2,int computeZ)
{
  size_t i,n=n1+n2;
  double *combined = (double *)VCalloc(n,sizeof(double));
  double *tmp = (double *)VCalloc(n,sizeof(double));
  for (i=0; i<n1; i++) combined[i] = x1[i];
  for (i=0; i<n2; i++) combined[i+n1] = x2[i];

  double xmedian = XMedian(combined,tmp,n);
  VFree(tmp);

  double a=0,b=0,c=0,d=0;
  for (i=0; i<n1; i++) {
    if (x1[i] > xmedian) a++;
    if (x1[i] < xmedian) c++;
  }
  for (i=0; i<n2; i++) {
    if (x2[i] > xmedian) b++;
    if (x2[i] < xmedian) d++;
  }
  VFree(combined);

  double nx = a+b+c+d;
  double s = (a+b)*(c+d)*(a+c)*(b+d);
  double chisq = (nx*(a*d-b*c)*(a*d-b*c)) / s;
  double z = sqrt(chisq);
  if (a < b) z = -z;

  if (computeZ == 1) return z;
  else {
    double pval = gsl_cdf_ugaussian_P(z);
    return pval;
  }
}


int LaminarMedian(gsl_vector *y,gsl_vector *mvec,
		  size_t numperm,long seed,gsl_vector *beta,gsl_vector *zval)
{
  size_t i;
  size_t perm;
  size_t n = y->size;
  double u,v,w,u0=0,v0=0,w0=0;

  gsl_vector_set_zero(beta);
  gsl_vector_set_zero(zval);

  /* get layer labels */
  size_t n0=0,n1=0,n2=0;
  int *label = (int *)VCalloc(n,sizeof(int));
  if (LayerLabels(mvec,label,&n0,&n1,&n2) < 0) {
    VFree(label);
    return -1;
  }
  
 
  /* no-permutation */
  double *y0 = (double *)VCalloc(n0,sizeof(double));
  double *y1 = (double *)VCalloc(n1,sizeof(double));
  double *y2 = (double *)VCalloc(n2,sizeof(double));

  double *tmp0 = (double *)VCalloc(n0,sizeof(double));
  double *tmp1 = (double *)VCalloc(n1,sizeof(double));
  double *tmp2 = (double *)VCalloc(n2,sizeof(double));
  AssignLabels(y,label,y0,y1,y2,n);

  beta->data[0] = XMedian(y0,tmp0,n0);
  beta->data[1] = XMedian(y1,tmp1,n1);
  beta->data[2] = XMedian(y2,tmp2,n2);

  /* skip permutations */
  if (numperm < 1) {
    int computeZ = 1;
    zval->data[0] = MedianTest(y0,y1,n0,n1,computeZ);
    zval->data[1] = MedianTest(y0,y2,n0,n2,computeZ);
    zval->data[2] = MedianTest(y1,y2,n1,n2,computeZ);
    goto noperm;
  }
  else {
    u0 = beta->data[0] - beta->data[1];
    v0 = beta->data[0] - beta->data[2];
    w0 = beta->data[1] - beta->data[2];
  }

  /* ini random generator */
  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *rx = gsl_rng_alloc(T);
  gsl_rng_set(rx,(unsigned long int) seed);

  int *tmplabel = (int *)VCalloc(n,sizeof(int));
  
  /* permute depth labels */
  double kx[3],xm0=0,xm1=0,xm2=0;
  kx[0]=kx[1]=kx[2]=0;
  for (perm=0; perm<numperm; perm++) {

    memcpy(tmplabel,label,n*sizeof(int));
    gsl_ran_shuffle(rx,tmplabel,n,sizeof(int));
    AssignLabels(y,tmplabel,y0,y1,y2,n);

    xm0 = XMedian(y0,tmp0,n0);
    xm1 = XMedian(y1,tmp1,n1);
    xm2 = XMedian(y2,tmp2,n2); 
    u = xm0-xm1;
    v = xm0-xm2;
    w = xm1-xm2;
    if (u > u0) kx[0]++;
    if (v > v0) kx[1]++;
    if (w > w0) kx[2]++;
  }
  gsl_rng_free(rx);
  VFree(tmp0);
  VFree(tmp1);
  VFree(tmp2);
  VFree(tmplabel);
  
  /* compute z-stats */
  double pmin = DBL_EPSILON;
  double pmax = 1.0-2.0*pmin;
  double nx = (double)numperm;
  double pval = 0;
  for (i=0; i<3; i++) {
    pval = kx[i]/nx;
    if (pval < pmin) pval = pmin;
    if (pval > pmax) pval = pmax;    
    zval->data[i] = gsl_cdf_ugaussian_Qinv(pval);
  }

 noperm: ;
  VFree(y0);
  VFree(y1);
  VFree(y2);
  VFree(label);
  return 1;
}

