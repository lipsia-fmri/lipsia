
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
#include <string.h>
#include <math.h>

#include "../cylutils/cyl.h"


extern double Xt2z(double t,double df);

/* two-sample Welsh test */
double TwosampleTest(double *x1,double *x2,size_t n1,size_t n2,int computeZ)
{
  double mean1 = gsl_stats_mean(x1,1,n1);
  double mean2 = gsl_stats_mean(x2,1,n2);
  
  double var1 = gsl_stats_variance(x1,1,n1);
  double var2 = gsl_stats_variance(x2,1,n2);

  double nx1 = (double)n1;
  double nx2 = (double)n2;

  double var = (var1/nx1 + var2/nx2);  
  double t = (mean1 - mean2)/sqrt(var);
  if (computeZ == 0) return t;

  double s1 = (var1*var1)/(nx1*nx1*(nx1-1.0));
  double s2 = (var2*var2)/(nx2*nx2*(nx2-1.0));
  double df = (var*var)/(s1 + s2);
  double z = Xt2z(t,df);
  return z;
}

void AssignLabels(gsl_vector *y,int *label,double *y0,double *y1,double *y2,size_t n)
{
  size_t i,k0=0,k1=0,k2=0;
  for (i=0; i<n; i++) {
    if (label[i] == 0) y0[k0++] = y->data[i];
    if (label[i] == 1) y1[k1++] = y->data[i];
    if (label[i] == 2) y2[k2++] = y->data[i];
  }
}


  /* get layer labels */
int LayerLabels(gsl_vector *mvec,int *label,size_t *n0,size_t *n1,size_t *n2)
{
  (*n0) = (*n1) = (*n2) = 0;
  size_t i,n=mvec->size;
  double t1 = 1.0/3.0;
  double t2 = 2.0/3.0;
  size_t nx0=0,nx1=0,nx2=0;
  for (i=0; i<n; i++) {
    if (mvec->data[i] < t1) { label[i] = 0; nx0++; }
    if (mvec->data[i] >= t1 && mvec->data[i] <= t2) { label[i] = 1; nx1++; }
    if (mvec->data[i] > t2) { label[i] = 2; nx2++; }
  }
  if (nx0 < 3 || nx1 < 3 || nx2 < 3) return -1;
  (*n0) = nx0;
  (*n1) = nx1;
  (*n2) = nx2;
  return 1;
}


int LaminarMean(gsl_vector *y,gsl_vector *mvec,
		size_t numperm,long seed,gsl_vector *beta,gsl_vector *zval)
{
  size_t i;
  size_t perm;
  size_t n = y->size;
  double u,v,w,u0,v0,w0;

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
  AssignLabels(y,label,y0,y1,y2,n);

  beta->data[0] = gsl_stats_mean(y0,1,n0);
  beta->data[1] = gsl_stats_mean(y1,1,n1);
  beta->data[2] = gsl_stats_mean(y2,1,n2);
 

  int computeZ=0;
  if (numperm < 1) computeZ=1;

  u0 = TwosampleTest(y0,y1,n0,n1,computeZ);
  v0 = TwosampleTest(y0,y2,n0,n2,computeZ);
  w0 = TwosampleTest(y1,y2,n1,n2,computeZ);

  /* skip permutations */
  if (numperm < 1) {
    zval->data[0] = u0;
    zval->data[1] = v0;
    zval->data[2] = w0;
    goto noperm;
  }


  /* ini random generator */
  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *rx = gsl_rng_alloc(T);
  gsl_rng_set(rx,(unsigned long int)seed);
  int *tmplabel = (int *)VCalloc(n,sizeof(int));
 
  /* permute layers labels */
  double kx[3];
  kx[0]=kx[1]=kx[2]=0;
  for (perm=0; perm<numperm; perm++) {

    memcpy(tmplabel,label,n*sizeof(int));
    gsl_ran_shuffle(rx,tmplabel,n,sizeof(int));
    AssignLabels(y,tmplabel,y0,y1,y2,n);

    u = TwosampleTest(y0,y1,n0,n1,computeZ);
    v = TwosampleTest(y0,y2,n0,n2,computeZ);
    w = TwosampleTest(y1,y2,n1,n2,computeZ);

    if (u > u0) kx[0]++;
    if (v > v0) kx[1]++;
    if (w > w0) kx[2]++;
  }

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
  gsl_rng_free(rx);
  VFree(tmplabel);

 noperm: ;
  VFree(y0);
  VFree(y1);
  VFree(y2);
  VFree(label);
  return 1;
}

