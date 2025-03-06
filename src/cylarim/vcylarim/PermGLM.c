
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
extern double contrastVariance(gsl_matrix *cov,gsl_vector *c);
extern void ZStats(gsl_vector *beta,gsl_matrix *cov,double edf,gsl_vector *zval);


double PermGLM(VImage zmap,VImage metric,VImage rim,
	       Cylinders *cyl,size_t cid,gsl_vector *beta,gsl_vector *zval)
{
  int b,r,c;
  size_t i,j,k,perm,numperm=1000;
  size_t nzval = zval->size;
  size_t nn = cyl->addr[cid]->size;
  size_t p = beta->size;
  double u=0,v=0;
  double rtcode = -1;
  gsl_matrix *X = NULL;
  gsl_multifit_linear_workspace *work = NULL;

  gsl_vector_set_zero(beta);
  gsl_vector_set_zero(zval);

  /* cylinder must be big enough for sufficient stats */
  if (nn < p*5) return -1;

  size_t n=0;
  for (i=0; i<nn; i++) {
    k = cyl->addr[cid]->data[i];
    b = gsl_matrix_int_get(cyl->xmap,k,0);
    r = gsl_matrix_int_get(cyl->xmap,k,1);
    c = gsl_matrix_int_get(cyl->xmap,k,2);
    if (VPixel(rim,b,r,c,VUByte) == 3) n++;
  }
  
  gsl_vector *y = gsl_vector_calloc(n);
  gsl_vector *mvec = gsl_vector_calloc(n);

  /* fill y-vector with activation map values */
  j=0;
  for (i=0; i<nn; i++) {
    k = cyl->addr[cid]->data[i];
    b = gsl_matrix_int_get(cyl->xmap,k,0);
    r = gsl_matrix_int_get(cyl->xmap,k,1);
    c = gsl_matrix_int_get(cyl->xmap,k,2);
    if (VPixel(rim,b,r,c,VUByte) != 3) continue;
    y->data[j] = VPixel(zmap,b,r,c,VFloat);
    mvec->data[j] = VPixel(metric,b,r,c,VFloat);
    j++;
  }

  /* return if zmap-values in this cylinder have no variance */
  double tss = gsl_stats_tss(y->data,1,y->size);
  if (tss < 0.001) goto ende;

  double edf = (double)(n - beta->size);
  double chisq=0;
  gsl_matrix *cov = gsl_matrix_calloc(p,p);
  work = gsl_multifit_linear_alloc(n,p);
  X = gsl_matrix_calloc(n,p);
  gsl_matrix_set_all(X,1.0);
  
  for (i=0; i<n; i++) {
    u = mvec->data[i];
    for (j=0; j<p-1; j++) {
      v = gaussian(u,(int)j);
      gsl_matrix_set(X,i,j,v);
    }
  }
  gsl_multifit_linear(X,y,beta,cov,&chisq,work);

  
  /* reciprocal condition number of design matrix */
  double rcond = gsl_multifit_linear_rcond(work);
  if (rcond < 0.001) goto ende;


  /* goodness of fit, R^2 */
  double r2 = 1.0 - chisq/tss;
  if (r2 < 0) goto ende;
  

  /* ini random generator */
  gsl_rng_env_setup();
  unsigned long int seed = 555;
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *rx = gsl_rng_alloc(T);
  gsl_rng_set(rx,(unsigned long int)seed);


  /* permutation testing */
  gsl_matrix *Xperm = gsl_matrix_calloc(n,p);
  gsl_matrix_set_all(Xperm,1.0);
  size_t pn = p-1;  /* do not permute the intercept */
  double *xx = (double *) VCalloc(pn,sizeof(double));

  gsl_vector *tval = gsl_vector_calloc(nzval);
  gsl_vector *tperm = gsl_vector_calloc(nzval);
  double *kx = (double *)VCalloc(nzval,sizeof(double));


  /* non-permuted stats */
  ZStats(beta,cov,edf,tval);
  
  /* permutations */
  for (perm=0; perm<numperm; perm++) {

    /* shuffle each row of the design matrix, keeping spatial correlations intact */
    for (i=0; i<n; i++) {
      for (j=0; j<pn; j++) xx[j] = gsl_matrix_get(X,i,j);
      gsl_ran_shuffle(rx,xx,pn,sizeof(double));
      for (j=0; j<pn; j++) gsl_matrix_set(Xperm,i,j,xx[j]);
    }
    gsl_multifit_linear(Xperm,y,beta,cov,&chisq,work);
    ZStats(beta,cov,edf,tperm);
       
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
  rtcode=r2;

  
  /* free memory */
  VFree(kx);
  VFree(xx);
  gsl_matrix_free(Xperm);
  gsl_matrix_free(cov);
  gsl_rng_free(rx);


 ende: ;
  if (work != NULL) gsl_multifit_linear_free(work);
  if (X != NULL) gsl_matrix_free(X);
  gsl_vector_free(y);
  gsl_vector_free(mvec);

  return rtcode;
}



