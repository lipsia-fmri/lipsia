/*
** PLS partial least squares, Algorithm SIMPLS
** X is multivariate (p>0), Y is univariate (m=1)
**
** Ref:
**  Sijmen de Jong, Chemometrics and Intelligent Lab Systems, 18 (1993), pp. 251-263
**
** G.Lohmann, Oct 2020
*/
#include "viaio/Vlib.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>




void VectorNorm(gsl_vector *x)
{
  double u = gsl_blas_dnrm2(x);
  gsl_vector_scale(x,1.0/u);
}


void OrthoVec(gsl_vector *vv,gsl_vector *sv)
{
  size_t i;
  double s,u,v;
  
  s=0;
  for (i=0; i<sv->size; i++) s += vv->data[i]*sv->data[i];      

  for (i=0; i<sv->size; i++) {
    u = sv->data[i];
    v = vv->data[i]*s;
    sv->data[i] = u-v;
  }
}


void FillMat(gsl_matrix *U,gsl_vector *x,size_t j)
{
  size_t i;
  for (i=0; i<U->size1; i++) {
    gsl_matrix_set(U,i,j,x->data[i]);
  }
}


double SIMPLS1(gsl_matrix *X,gsl_vector *y,gsl_vector *beta,size_t A)
{
  size_t i,a;
  double meant=0,normt=0,q=0;

  if (A < 1) A = X->size2;
  if (A >= X->size2) VError("SIMPLS1, A");

  size_t n = X->size1;
  size_t p = X->size2;
  
  gsl_vector *sv = gsl_vector_calloc(p);
  gsl_blas_dgemv(CblasTrans,1,X,y,0,sv);

  gsl_vector *rv = gsl_vector_calloc(p);
  gsl_vector *tv = gsl_vector_calloc(n);
  gsl_vector *pv = gsl_vector_calloc(p);
  gsl_vector *uv = gsl_vector_calloc(n);
  gsl_vector *vv = gsl_vector_calloc(p);

  gsl_vector *tmp = gsl_vector_calloc(A);
  gsl_vector *tmpv = gsl_vector_calloc(vv->size);
  gsl_vector *tmpu = gsl_vector_calloc(uv->size);

  gsl_vector *qv = gsl_vector_calloc(A);
  gsl_matrix *R = gsl_matrix_calloc(pv->size,A);
  gsl_matrix *V = gsl_matrix_calloc(vv->size,A);
  gsl_matrix *U = gsl_matrix_calloc(uv->size,A);
  gsl_matrix *P = gsl_matrix_calloc(pv->size,A);
  gsl_matrix *T = gsl_matrix_calloc(tv->size,A);


  for (a=0; a<A; a++) {

    gsl_blas_ddot(sv,sv,&q);
    for (i=0; i<p; i++) rv->data[i] = q*sv->data[i];
 
    gsl_blas_dgemv(CblasNoTrans,1,X,rv,0,tv);
    meant = gsl_stats_mean(tv->data,1,tv->size);
    gsl_vector_add_constant(tv,-meant);

    normt = gsl_blas_dnrm2(tv);
    gsl_vector_scale(tv,1.0/normt);
    gsl_vector_scale(rv,1.0/normt);

    gsl_blas_dgemv(CblasTrans,1,X,tv,0,pv);

    gsl_blas_ddot(y,tv,&q);
    for (i=0; i<n; i++) uv->data[i] = q * y->data[i];
 
    gsl_vector_memcpy(vv,pv);
    if (a > 0) {
      gsl_blas_dgemv(CblasTrans,1,V,pv,0,tmp);
      gsl_blas_dgemv(CblasNoTrans,1,V,tmp,0,tmpv);
      gsl_vector_sub(vv,tmpv);
      
      gsl_blas_dgemv(CblasTrans,1,T,uv,0,tmp);
      gsl_blas_dgemv(CblasNoTrans,1,T,tmp,0,tmpu);
      gsl_vector_sub(uv,tmpu);
    }
    VectorNorm(vv);
    OrthoVec(vv,sv);

    qv->data[a] = q;
    FillMat(R,rv,a);
    FillMat(T,tv,a);
    FillMat(P,pv,a);
    FillMat(U,uv,a);
    FillMat(V,vv,a);
  }
  gsl_vector_free(tmpv);
  gsl_vector_free(tmpu);
  gsl_vector_free(tmp);


  /* compute beta */
  gsl_blas_dgemv(CblasNoTrans,1,R,qv,0.0,beta);

 
  /* variance explained in Y */
  /*
    gsl_blas_ddot(qv,qv,&yvar);
    double yvar /= (double)(n-1);
  */
  
  /* variance explained in X, factor loadings in X */
  /*
    double xvar=0;
    double nx = (double)n;
    gsl_matrix *PX = gsl_matrix_calloc(A,A);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,P,P,0.0,PX);
    for (i=0; i<A; i++) xloads->data[i] = gsl_matrix_get(PX,i,i)/nx;
    gsl_matrix_free(PX);
    gsl_blas_ddot(xloads,xloads,&xvar);
    xvar /= (double)(n-1);
  */

  /* free */
  gsl_vector_free(sv);
  gsl_vector_free(qv);
  gsl_vector_free(rv);
  gsl_vector_free(tv);
  gsl_vector_free(pv);
  gsl_vector_free(uv);
  gsl_vector_free(vv);

  gsl_matrix_free(R);
  gsl_matrix_free(V);
  gsl_matrix_free(U);
  gsl_matrix_free(P);
  gsl_matrix_free(T);
  return 0;
}

