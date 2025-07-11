
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

/* 
** approximation: convert t to z values 
*/
double Xt2z(double t,double df)
{
  double z=0,u;
  if (df < 0) return t;

  u = df*log(1.0+t*t/df)*(1.0-0.5/df);
  if (u <= 0) return 0;
  z = sqrt(u);
  if (t < 0) z = -z;
  return z;
}

extern void AssignLabels(gsl_vector *y,int *label,double *y0,double *y1,double *y2,size_t n);
extern int LayerLabels(gsl_vector *mvec,int *label,size_t *n0,size_t *n1,size_t *n2);


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


double LaminarMean(gsl_vector *ya,gsl_vector *mvec0,
		   size_t numperm,long seed,gsl_vector *beta,gsl_vector *zval)
{
  size_t i,j;
  size_t perm;
  double u,v,w,u0=0,v0=0,w0=0;

  gsl_vector_set_zero(beta);
  gsl_vector_set_zero(zval);

  size_t n = 0;
  for (i=0; i<mvec0->size; i++) {
    if (mvec0->data[i] > 0.05 && mvec0->data[i] < 0.95) n++;
  }
  gsl_vector *y = gsl_vector_calloc(n);
  gsl_vector *mvec = gsl_vector_calloc(n);

  j=0;
  for (i=0; i<mvec0->size; i++) {
    if (mvec0->data[i] > 0.05 && mvec0->data[i] < 0.95) {
      y->data[j] = ya->data[i];
      mvec->data[j] = mvec0->data[i];
      j++;
    }
  }

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
  
  beta->data[0] = gsl_stats_mean(y0,1,n0);
  beta->data[1] = gsl_stats_mean(y1,1,n1);
  beta->data[2] = gsl_stats_mean(y2,1,n2);

  /* skip permutations */
  if (numperm < 1) {
    int computeZ = 1;
    zval->data[0] = TwosampleTest(y0,y1,n0,n1,computeZ);
    zval->data[1] = TwosampleTest(y0,y2,n0,n2,computeZ);
    zval->data[2] = TwosampleTest(y1,y2,n1,n2,computeZ);
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

  double *tmp_01 = (double *)VCalloc(n0+n1,sizeof(double));
  double *tmp_02 = (double *)VCalloc(n0+n2,sizeof(double));
  double *tmp_12 = (double *)VCalloc(n1+n2,sizeof(double));

  double *xtmp0 = (double *)VCalloc(n0,sizeof(double));
  double *xtmp1 = (double *)VCalloc(n1,sizeof(double));
  double *xtmp2 = (double *)VCalloc(n2,sizeof(double));
  

  /* permute depth labels */
  double pos[3],neg[3],xm0=0,xm1=0,xm2=0;
  for (i=0; i<3; i++) pos[i] = neg[i] = 0;
  for (perm=0; perm<numperm; perm++) {
    
    /* deep-middle */
    j=0;
    for (i=0; i<n0; i++) tmp_01[j++] = y0[i];
    for (i=0; i<n1; i++) tmp_01[j++] = y1[i];
    gsl_ran_shuffle(rx,tmp_01,n0+n1,sizeof(double));
    for (i=0; i<n0; i++) tmp0[i] = tmp_01[i];
    for (i=0; i<n1; i++) tmp1[i] = tmp_01[i+n0];
    xm0 = gsl_stats_mean(tmp0,1,n0);
    xm1 = gsl_stats_mean(tmp1,1,n1);
    u = xm0 - xm1;

    /* deep-superficial */
    j=0;
    for (i=0; i<n0; i++) tmp_02[j++] = y0[i];
    for (i=0; i<n2; i++) tmp_02[j++] = y2[i];
    gsl_ran_shuffle(rx,tmp_02,n0+n2,sizeof(double));
    for (i=0; i<n0; i++) tmp0[i] = tmp_02[i];
    for (i=0; i<n2; i++) tmp2[i] = tmp_02[i+n2];
    xm0 = gsl_stats_mean(tmp0,1,n0);
    xm2 = gsl_stats_mean(tmp2,1,n2);
    v = xm0 - xm2;
    
    /* middle-superficial */
    j=0;
    for (i=0; i<n1; i++) tmp_12[j++] = y1[i];
    for (i=0; i<n2; i++) tmp_12[j++] = y2[i];
    gsl_ran_shuffle(rx,tmp_12,n1+n2,sizeof(double));
    for (i=0; i<n1; i++) tmp1[i] = tmp_12[i];
    for (i=0; i<n2; i++) tmp2[i] = tmp_12[i+n1];
    xm1 = gsl_stats_mean(tmp1,1,n1);
    xm2 = gsl_stats_mean(tmp2,1,n2);
    w = xm1 - xm2;

    /* permuted test statistic */
    if (u > u0) pos[0]++;
    if (v > v0) pos[1]++;
    if (w > w0) pos[2]++;

    if (u < u0) neg[0]++;
    if (v < v0) neg[1]++;
    if (w < w0) neg[2]++;
  }
  gsl_rng_free(rx);
  VFree(tmp0);
  VFree(tmp1);
  VFree(tmp2);
  VFree(xtmp0);
  VFree(xtmp1);
  VFree(xtmp2);
  VFree(tmp_01);
  VFree(tmp_02);
  VFree(tmp_12);
  
  
  /* compute z-stats */
  double nx = (double)numperm;
  double z=0,p0=0,p1=0;
  for (i=0; i<3; i++) {
    p0 = (pos[i]+1.0)/(nx+1.0);
    p1 = (neg[i]+1.0)/(nx+1.0);
    z=0;
    if (p0 < p1) z = gsl_cdf_ugaussian_Qinv(p0);
    else z = -gsl_cdf_ugaussian_Qinv(p1);
    if (gsl_finite(z)==0) z = 0;
    zval->data[i] = z;
  }

 noperm: ;
  gsl_vector_free(y);
  gsl_vector_free(mvec);
  VFree(y0);
  VFree(y1);
  VFree(y2);
  VFree(label);
  return 1;
}

