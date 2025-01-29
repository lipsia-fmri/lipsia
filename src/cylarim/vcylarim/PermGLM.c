
/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/Vlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
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
#include <math.h>

#include "../cylutils/cyl.h"

extern double gaussian(double x,int);
extern double EDF(gsl_multifit_linear_workspace *work);
extern void TStats(gsl_vector *beta,gsl_matrix *cov,gsl_vector *t);


double PermGLM(VImage zmap,VImage metric,
	       Cylinders *cyl,size_t cid,gsl_vector *beta,gsl_vector *zval,double *edf)
{
  int b,r,c;
  size_t i,j,k,perm,numperm=2000;
  size_t dim=beta->size;
  size_t n=cyl->addr[cid]->size;
  double u=0,v=0,w=0,p=0;
  double rtcode = -1;
  gsl_matrix *X = NULL;

  /* cylinder must be big enough for sufficient stats */
  if (n < dim*10) return -1;
  
  gsl_vector_set_zero(beta);
  gsl_vector *y = gsl_vector_calloc(n);
  gsl_vector *mvec = gsl_vector_calloc(n);

  /* fill y-vector with activation map values */
  for (i=0; i<n; i++) {
    k = cyl->addr[cid]->data[i];
    b = gsl_matrix_int_get(cyl->xmap,k,0);
    r = gsl_matrix_int_get(cyl->xmap,k,1);
    c = gsl_matrix_int_get(cyl->xmap,k,2);
    y->data[i] = VPixel(zmap,b,r,c,VFloat);
    mvec->data[i] = VPixel(metric,b,r,c,VFloat);
  }

  /* return if z-values in this cylinder have no variance */
  double tss = gsl_stats_tss(y->data,1,y->size);
  if (tss < 0.001) goto ende;

  
  /* random generator */
  gsl_rng_env_setup();
  unsigned long int seed = 555;
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *rx = gsl_rng_alloc(T);
  gsl_rng_set(rx,(unsigned long int)seed);


  double chisq=0;
  double *base = (double *)VCalloc(n,sizeof(double));
  gsl_vector *residuals = gsl_vector_calloc(n);
  gsl_vector *xbeta = gsl_vector_calloc(dim);
  gsl_vector *yx = gsl_vector_calloc(n);
  gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n,dim);
  gsl_matrix *cov = gsl_matrix_calloc(dim,dim);
  X = gsl_matrix_calloc(n,dim);
  gsl_matrix_set_all(X,1.0);
  
  for (i=0; i<n; i++) {
    u = mvec->data[i];
    for (j=0; j<3; j++) {
      v = gaussian(u,(int)j);
      gsl_matrix_set(X,i,j,v);
    }
  }
  gsl_multifit_linear(X,y,beta,cov,&chisq,work);

  
  /* reciprocal condition number of design matrix */
  double rcond = gsl_multifit_linear_rcond(work);
  if (rcond < 0.001) return -1;

    
  /* goodness of fit, R^2 */
  double r2 = 1.0 - chisq/tss;
  if (r2 < 0) goto ende;

      
  /* estimate noise variance, low noise results have higher impact on averages */
  (*edf) = EDF(work);
  double noise_variance = chisq/(*edf);
  w=0;
  if (fabs(noise_variance) > 0) w = 1.0/noise_variance;
  rtcode = w;


  /* residual permutation test */
  gsl_vector_int *xperm = gsl_vector_int_calloc(zval->size);
  gsl_vector *tval = gsl_vector_calloc(3);
  gsl_vector *tperm = gsl_vector_calloc(3);
  TStats(beta,cov,tval);

  
  gsl_multifit_linear_residuals(X,y,beta,residuals);
  for (perm=0; perm<numperm; perm++) {

    for (i=0; i<n; i++) base[i] = residuals->data[i];
    gsl_ran_shuffle(rx,base,n,sizeof(double));
    for (i=0; i<n; i++) yx->data[i] = y->data[i]+base[i];
    gsl_multifit_linear(X,yx,xbeta,cov,&chisq,work);
    gsl_multifit_linear_residuals(X,yx,xbeta,residuals);

    TStats(xbeta,cov,tperm);

    for (i=0; i<3; i++) {
      if (tperm->data[i] > tval->data[i]) xperm->data[i]++;
    }
  }

  for (i=0; i<3; i++) {
    p = ((double)xperm->data[i])/((double)numperm);
    if (p < 0.0001) p = 0.0001;
    if (p > 0.9999) p = 0.9999;    
    zval->data[i] = gsl_cdf_ugaussian_Qinv(p);
  }



  /* free memory */
  gsl_multifit_linear_free(work);
  gsl_matrix_free(X);
  gsl_matrix_free(cov);
  gsl_vector_free(yx);
  gsl_vector_free(xbeta);
  gsl_vector_int_free(xperm);
  gsl_vector_free(residuals);
  gsl_rng_free(rx);
  VFree(base);

 ende: ;
  gsl_vector_free(y);
  gsl_vector_free(mvec);

  return rtcode;
}



