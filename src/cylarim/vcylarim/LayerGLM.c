
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


/* From the standard C library: */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../cylutils/cyl.h"

extern double gaussian(double x,int);


/* no model, simple averaging */
void Simple(gsl_vector *y,gsl_vector *mvec,gsl_vector *beta)
{
  size_t i;
  double s0=0,s1=0,s2=0,nx0=0,nx1=0,nx2=0; 
  for (i=0; i<y->size; i++){
    if (mvec->data[i] < 0.333) { s0 += y->data[i]; nx0++; }
    if (mvec->data[i] >= 0.333 && mvec->data[i] < 0.666) { s1 += y->data[i]; nx1++; }
    if (mvec->data[i] >= 0.666) { s2 += y->data[i]; nx2++; }
  }
  gsl_vector_set_zero(beta);
  if (nx0 > 0.5) beta->data[0] = s0/nx0;
  if (nx1 > 0.5) beta->data[1] = s1/nx1;
  if (nx2 > 0.5) beta->data[2] = s2/nx2;
}



/* effective degrees of freedom as trace(I-H), where H is the "hat"-matrix */
double EDF(gsl_multifit_linear_workspace *work)
{
  size_t i;
  size_t p = work->p;
  size_t n = work->n;
  gsl_matrix *H = gsl_matrix_calloc(p,p);
  gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,work->Q,work->Q,0.0,H);
  double trace = 0;
  for (i=0; i<p; i++) trace += gsl_matrix_get(H,i,i);
  gsl_matrix_free(H);
  double edf = (double)(n) - trace;
  return edf;
}


double LayerGLM(VImage zmap,VImage metric,
		Cylinders *cyl,size_t cid,gsl_vector *beta,gsl_vector *bcov,double *edf,int model)
{
  int b,r,c;
  size_t i,j,k;
  size_t dim=beta->size;
  size_t n=cyl->addr[cid]->size;
  double u=0,v=0,w=0;
  double rtcode = -1;
  gsl_matrix *X = NULL;

  /* cylinder must be big enough for sufficient stats */
  if (n < dim*10) return -1;
  
  gsl_vector_set_zero(beta);
  gsl_vector_set_zero(bcov);
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


  switch (model) {

  case 0:  /* no model, simple averaging */
    Simple(y,mvec,beta);
    w = (double)n;
    break;
    

  case 1:     /* GLM */
    ;
    double chisq=0;
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

    /* covariances */
    bcov->data[0] = gsl_matrix_get(cov,0,0);
    bcov->data[1] = gsl_matrix_get(cov,0,1);
    bcov->data[2] = gsl_matrix_get(cov,0,2);
    bcov->data[3] = gsl_matrix_get(cov,1,1);
    bcov->data[4] = gsl_matrix_get(cov,1,2);
    bcov->data[5] = gsl_matrix_get(cov,2,2);

      
    /* beta variance, t-test */
    /*
    for (j=0; j<3; j++) betavar->data[j] = gsl_matrix_get(cov,j,j);
    t = (beta->data[0]-beta->data[1]) /
      sqrt(gsl_matrix_get(cov,0,0) + gsl_matrix_get(cov,1,1) - 2.0*gsl_matrix_get(cov,0,1));
    */
   
    gsl_multifit_linear_free(work);
    gsl_matrix_free(cov);
    break;

  default:
    VError(" unknown model");
  }
  rtcode = w;

  
  /* free memory */
 ende: ;
  gsl_vector_free(y);
  gsl_vector_free(mvec);
  if (X !=NULL) gsl_matrix_free(X);
  return rtcode;
}



