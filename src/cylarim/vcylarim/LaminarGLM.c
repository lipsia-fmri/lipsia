/*
** cylarim, residual permutation GLM 
**
** G.Lohmann, Jan 2025, MPI-KYB
*/
/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/Vlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>


/* From the standard C library: */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "../cylutils/cyl.h"

extern double gaussian(double x,int);
extern void XWriteOutput(VImage image,VAttrList geolist,char *filename);
extern void GetResolution(VImage src,gsl_vector *reso);
extern void Leverage(gsl_matrix *X,gsl_vector *leverage);


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



/* Laminar GLM */
double LaminarGLM(VImage zmap,VImage metric,VImage rim,
		  Cylinders *cyl,size_t cid,double upsampling,gsl_vector *beta,gsl_vector *zval)
{
  int b,r,c;
  size_t i,j,k;
  size_t p=beta->size;
  size_t nn=cyl->addr[cid]->size;
  size_t m=zval->size;
  double u=0,v=0,edf=0,r2=-1;
 
  gsl_vector_set_zero(beta);
  gsl_vector_set_zero(zval);
   
  /* cylinder must be big enough for sufficient stats */
  if (nn < p*10) return -1;

  size_t n=0;
  for (i=0; i<nn; i++) {
    k = cyl->addr[cid]->data[i];
    b = gsl_matrix_int_get(cyl->xmap,k,0);
    r = gsl_matrix_int_get(cyl->xmap,k,1);
    c = gsl_matrix_int_get(cyl->xmap,k,2);
    if (VPixel(rim,b,r,c,VUByte) == 3) n++;
  }

  /* fill y-vector with activation map values */
  gsl_vector *y = gsl_vector_calloc(n);
  gsl_vector *mvec = gsl_vector_calloc(n);
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

  /* return if z-values in this cylinder have no variance */
  double tss = gsl_stats_tss(y->data,1,y->size);
  if (tss < 0.001) goto ende;
  

  /* GLM */
  double chisq=0;
  gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n,p);
  gsl_matrix *cov = gsl_matrix_calloc(p,p);
  gsl_matrix *X = gsl_matrix_calloc(n,p);
  gsl_matrix_set_all(X,1.0);
   
  for (i=0; i<n; i++) {
    u = mvec->data[i];
    for (j=0; j<m; j++) {
      v = gaussian(u,(int)j);
      gsl_matrix_set(X,i,j,v);
    }
  }
  gsl_multifit_linear(X,y,beta,cov,&chisq,work);


  /* reciprocal condition number of design matrix */
  double rcond = gsl_multifit_linear_rcond(work);
  if (rcond < 0.001) goto ende;
    
  /* goodness of fit, R^2 */
  r2 = 1.0 - chisq/tss;
  if (r2 < 0) goto ende;

 
  /* degrees of freedom */
  double nx = (double)n;
  double px = (double)p;
  double neff = nx;

  /* adjust number of observations using upsampling factor */
  if (upsampling > DBL_EPSILON) neff = nx/upsampling;
  edf = neff - px;


  /*  z-values */
  ZStats(beta,cov,edf,zval);
  gsl_matrix_free(cov);
  gsl_multifit_linear_free(work);
  gsl_matrix_free(X);  

 ende: ;
  gsl_vector_free(y);
  gsl_vector_free(mvec);
  return r2;
}



