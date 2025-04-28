
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

void YNorm(gsl_vector *x)
{
  double mean = gsl_stats_mean(x->data,1,x->size);
  gsl_vector_add_constant(x,-mean);
}

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


/* contrast variance, c^t cov c */
double contrastVariance(gsl_matrix *cov,gsl_vector *c)
{
  size_t i,j;

  double cv = 0.0;
  for (i=0; i<cov->size1; i++) {
    for (j=0; j<cov->size2; j++) {
      cv += c->data[i] * gsl_matrix_get(cov,i,j) * c->data[j];
    }
  }
  return cv;
}


/* z-stats for three contrasts: b0-b1, b0-b2, b1-b2 */
void ZStats(gsl_vector *beta,gsl_matrix *cov,double edf,gsl_vector *zval)
{
  double cv=0,t=0,tiny=0.0001;
  gsl_vector *contrast = gsl_vector_calloc(beta->size);
  gsl_vector_set_zero(zval);
  
  /* contrast 0 */
  gsl_vector_set_zero(contrast);
  contrast->data[0] = 1;
  contrast->data[1] = -1;
  cv = contrastVariance(cov,contrast);
  if (cv > tiny) {
    t = (beta->data[0] - beta->data[1])/sqrt(cv);
    zval->data[0] = Xt2z(t,edf);
  }
  
  /* contrast 1 */
  gsl_vector_set_zero(contrast);
  contrast->data[0] = 1;
  contrast->data[2] = -1;
  cv = contrastVariance(cov,contrast);
  if (cv > tiny) {
    t = (beta->data[0] - beta->data[2])/sqrt(cv);
    zval->data[1] = Xt2z(t,edf);
  }
  
  /* contrast 2 */
  gsl_vector_set_zero(contrast);
  contrast->data[1] = 1;
  contrast->data[2] = -1;
  cv = contrastVariance(cov,contrast);
  if (cv > tiny) {
    t = (beta->data[1] - beta->data[2])/sqrt(cv);
    zval->data[2] = Xt2z(t,edf);
  }
  gsl_vector_free(contrast);
}



/* gaussian */
double gaussian(double x,int i)
{
  double mean[3] = {0.05,0.5,0.95};
  double sigma=0.2;
  double z = exp(-(x-mean[i])*(x-mean[i])/(2.0*sigma*sigma));
  return z;
}



int LaminarGLM(gsl_vector *y,gsl_vector *mvec,
	       size_t numperm,long seed,gsl_vector *beta,gsl_vector *zval)
{
  int rtcode = -1;
  size_t i,j,perm;
  size_t nzval = zval->size;
  size_t p = beta->size;
  size_t n = y->size;
  double u=0,v=0;
  gsl_matrix *X = NULL;
  gsl_multifit_linear_workspace *work = NULL;

  gsl_vector_set_zero(beta);
  gsl_vector_set_zero(zval);

  /* centering, no intercept */
  YNorm(y);

  /* set up GLM */
  double chisq=0;
  gsl_vector *tval = gsl_vector_calloc(nzval);
  gsl_matrix *cov = gsl_matrix_calloc(p,p);
  work = gsl_multifit_linear_alloc(n,p);
  X = gsl_matrix_calloc(n,p);
  gsl_matrix_set_all(X,1.0);
  
  for (i=0; i<n; i++) {
    u = mvec->data[i];
    for (j=0; j<beta->size; j++) {
      v = gaussian(u,(int)j);
      gsl_matrix_set(X,i,j,v);
    }
  }

  /* GLM */
  gsl_multifit_linear(X,y,beta,cov,&chisq,work);

  
  /* reciprocal condition number of design matrix */
  double rcond = gsl_multifit_linear_rcond(work);
  if (rcond < 0.001) goto noperm;


  /* goodness of fit, R^2 */
  double tss = gsl_stats_tss(y->data,1,y->size);
  double r2 = 1.0 - chisq/tss;
  if (r2 < 0) goto noperm;

  
  /* non-permuted stats */
  double edf = -1;
  ZStats(beta,cov,edf,tval);
  
  if (numperm < 1) {  /* no permutation, degrees of freedom over-estimated */
    edf = (double)(n-p);
    ZStats(beta,cov,edf,zval);
    goto noperm;
  }
  

  /* ini random generator */
  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *rx = gsl_rng_alloc(T);
  gsl_rng_set(rx,(unsigned long int)seed);

  /* permutation testing */
  gsl_matrix *Xperm = gsl_matrix_calloc(n,p);
  gsl_matrix_set_all(Xperm,1.0);
  gsl_vector *xbeta = gsl_vector_calloc(beta->size);
  gsl_vector *tperm = gsl_vector_calloc(nzval);
  double *kx = (double *)VCalloc(nzval,sizeof(double));
  double *xtmp = (double *)VCalloc(n,sizeof(double));

  
  /* permutations */
  for (perm=0; perm<numperm; perm++) {

    /* shuffle depth values */
    memcpy(xtmp,mvec->data,n*sizeof(double));
    gsl_ran_shuffle(rx,xtmp,n,sizeof(double));
    
    for (i=0; i<n; i++) {
      u = xtmp[i];
      for (j=0; j<p; j++) {
	v = gaussian(u,(int)j);
	gsl_matrix_set(Xperm,i,j,v);
      }
    }
  
    gsl_multifit_linear(Xperm,y,xbeta,cov,&chisq,work);
    ZStats(xbeta,cov,edf,tperm);
       
    for (i=0; i<nzval; i++) {
      if (tperm->data[i] > tval->data[i]) kx[i]++;
    }
  }

  /* compute z-stats */
  double pmin = DBL_EPSILON;
  double pmax = 1.0-2.0*pmin;
  double nx = (double)numperm;
  double pval = 0;
  for (i=0; i<nzval; i++) {
    pval = kx[i]/nx;
    if (pval < pmin) pval = pmin;
    if (pval > pmax) pval = pmax;    
    zval->data[i] = gsl_cdf_ugaussian_Qinv(pval);
  }
  rtcode = 1;
  
  /* free memory */
  VFree(kx);
  VFree(xtmp);
  gsl_vector_free(xbeta);
  gsl_vector_free(tperm);
  gsl_matrix_free(Xperm);
  gsl_rng_free(rx);
  

 noperm: ;
  gsl_multifit_linear_free(work);
  gsl_vector_free(tval);
  gsl_matrix_free(cov);
  gsl_matrix_free(X);
  return rtcode;
}
