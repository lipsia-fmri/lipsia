
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

extern void SIMPLS1(gsl_matrix *X,gsl_vector *y,gsl_vector *beta,size_t A);
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


/* 
** approximation: convert t to z values 
*/
double Xt2z(double t,double df)
{
  double z=0,u;

  u = df*log(1.0+t*t/df)*(1.0-0.5/df);
  if (u <= 0) return 0;
  z = sqrt(u);
  if (t < 0) z = -z;
  return z;
}


void TStats(gsl_vector *beta,gsl_matrix *cov,gsl_vector *t)
{
  double bx0=beta->data[0];
  double bx1=beta->data[1];
  double bx2=beta->data[2];
  
  double c00=gsl_matrix_get(cov,0,0);
  double c01=gsl_matrix_get(cov,0,1);
  double c02=gsl_matrix_get(cov,0,2);
  double c11=gsl_matrix_get(cov,1,1);
  double c12=gsl_matrix_get(cov,1,2);
  double c22=gsl_matrix_get(cov,2,2);
  
  t->data[0] = (bx0 - bx1)/sqrt(c00 + c11 - 2.0*c01);
  t->data[1] = (bx0 - bx2)/sqrt(c00 + c22 - 2.0*c02);
  t->data[2] = (bx1 - bx2)/sqrt(c11 + c22 - 2.0*c12);
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


double LayerGLM(VImage zmap,VImage metric,Cylinders *cyl,size_t cid,int model,
		gsl_vector *beta,gsl_vector *zval,double *edf)
{
  int b,r,c;
  size_t i,j,k;
  size_t dim=beta->size;
  size_t n=cyl->addr[cid]->size;
  double u=0,v=0,w=0;
  double rtcode = -1;
  gsl_matrix *X = NULL;

  *edf=0;
  gsl_vector_set_zero(beta);
  gsl_vector_set_zero(zval);

  
  /* cylinder must be big enough for sufficient stats */
  if (n < dim*10) return -1;
  
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


    /* z-values */
    gsl_vector *tval = gsl_vector_calloc(3);
    TStats(beta,cov,tval);

    for (i=0; i<3; i++) {      
      zval->data[i] = Xt2z(tval->data[i],(*edf));
    }

    gsl_multifit_linear_free(work);
    gsl_matrix_free(cov);
    gsl_vector_free(tval);
    break;

    
  default:
    VError(" unknown model");
  }
  rtcode = w;

  
  /* free memory */
 ende: ;
  gsl_vector_free(y);
  gsl_vector_free(mvec);
  if (X!=NULL) gsl_matrix_free(X);
  return rtcode;
}



