#include <viaio/Vlib.h>
#include <viaio/VImage.h>

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_fit.h>


extern void MatColumnCentering(gsl_matrix *X);
extern void ZNorm(gsl_vector *);
extern void YNorm(gsl_vector *);
extern double PearsonCorr(gsl_vector *ytest,gsl_vector *result);
extern double CoD(gsl_vector *ytest,gsl_vector *result);
extern void MatNorm0(gsl_matrix *X);
extern void MatNorm1(gsl_matrix *,gsl_matrix *);
extern void KNorm(gsl_matrix *X);
extern void printmat(gsl_matrix *A);
extern void RankNorm(gsl_vector *x);
extern double Accuracy(gsl_vector *ytest,gsl_vector *result);



void CorrMat(gsl_matrix *X,gsl_matrix *C)
{
  size_t i,j;
  double *ptr0,*ptr1,s=0;

  for (i=0; i<X->size1; i++) {
    ptr0 = gsl_matrix_ptr(X,i,0);
    for (j=0; j<i; j++) {
      ptr1 = gsl_matrix_ptr(X,j,0);
      s = gsl_stats_correlation(ptr0,1,ptr1,1,X->size2);
      gsl_matrix_set(C,i,j,s);
      gsl_matrix_set(C,j,i,s);
    }
  }
  for (i=0; i<C->size1; i++) gsl_matrix_set(C,i,i,1.0);
}


/* univariate fit */
void PredictFit(gsl_vector *xtmp,gsl_vector *ytrain,gsl_vector *xtrain,gsl_vector *xtest)
{
  size_t i;
  double c0,c1,cov00,cov01,cov11,sumsq,y,y_err;

  gsl_fit_linear(xtmp->data,1,ytrain->data,1,xtrain->size,&c0,&c1,&cov00,&cov01,&cov11,&sumsq);
  for (i=0; i<xtrain->size; i++) {
    gsl_fit_linear_est(xtrain->data[i],c0,c1,cov00,cov01,cov11,&y,&y_err);
    xtrain->data[i] = y;
  }
  for (i=0; i<xtest->size; i++) {
    gsl_fit_linear_est(xtest->data[i],c0,c1,cov00,cov01,cov11,&y,&y_err);
    xtest->data[i] = y;
  }
}


void Fit(gsl_vector *xdata,gsl_vector *ydata,gsl_vector *zdata)
{
  size_t i;
  double c0,c1,cov00,cov01,cov11,sumsq,y,y_err;
  if (xdata->size != ydata->size) VError(" Fit, dim err");
  gsl_fit_linear(xdata->data,1,ydata->data,1,xdata->size,&c0,&c1,&cov00,&cov01,&cov11,&sumsq);
  for (i=0; i<xdata->size; i++) {
    gsl_fit_linear_est(xdata->data[i],c0,c1,cov00,cov01,cov11,&y,&y_err);
    zdata->data[i] = y;
  }
}


double XFit(gsl_matrix *X,gsl_vector *ytrain,gsl_vector *ytest,
	    gsl_vector *xtrain,gsl_vector *xtest,int *train,int *test)
{
  if (X == NULL) return 0;
  
  size_t i,j;
  size_t ntrain = xtrain->size;
  size_t ntest = xtest->size;
  size_t dimX = X->size2;
  double chisq=0;
  
  if (X->size1 == 1) {
    gsl_vector *xtmp = gsl_vector_calloc(ntrain);
    for (i=0; i<ntrain; i++) xtmp->data[i] = gsl_matrix_get(X,train[i],0);
    PredictFit(xtmp,ytrain,xtrain,xtest);
    gsl_vector_free(xtmp);
    return 0;
  }
  
  
  gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(ntrain,dimX);
  gsl_matrix *cov = gsl_matrix_alloc (X->size2,X->size2);
  gsl_vector *beta = gsl_vector_calloc(X->size2);
  gsl_matrix *Xtrain = gsl_matrix_calloc(ntrain,dimX);
  gsl_matrix *Xtest = gsl_matrix_calloc(ntest,dimX);

  for (j=0; j<dimX; j++) {
    for (i=0; i<ntrain; i++) {
      gsl_matrix_set(Xtrain,i,j,gsl_matrix_get(X,train[i],j));
    }
    for (i=0; i<ntest; i++) {
      gsl_matrix_set(Xtest,i,j,gsl_matrix_get(X,test[i],j));
    }
  }
  MatColumnCentering(Xtrain);
  MatColumnCentering(Xtest);

  gsl_multifit_linear(Xtrain,ytrain,beta,cov,&chisq,work);
  gsl_blas_dgemv(CblasNoTrans,1,Xtrain,beta,0,xtrain);
  gsl_blas_dgemv(CblasNoTrans,1,Xtest,beta,0,xtest);

  gsl_multifit_linear_free(work);
  gsl_matrix_free(Xtrain);
  gsl_matrix_free(Xtest);
  gsl_matrix_free(cov);
  gsl_vector_free(beta);
  return 0;
}


gsl_vector *Multifit(gsl_matrix *X,gsl_vector *regressor)
{
  if (X==NULL) return NULL;
  
  size_t i,j,dimX=X->size2;
  size_t nsubj=regressor->size;
  double z=0,r2=0,chisq=0;


  double mean = gsl_stats_mean(regressor->data,1,regressor->size);
  gsl_vector_add_constant(regressor,-mean);

  gsl_vector *ztmp = gsl_vector_calloc(nsubj);
  for (j=0; j<X->size2; j++) {
    for (i=0; i<X->size1; i++) ztmp->data[i] = gsl_matrix_get(X,i,j);
    z = PearsonCorr(ztmp,regressor);
    r2 = CoD(regressor,ztmp);
  }


  gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(X->size1,X->size2);
  gsl_matrix *cov = gsl_matrix_alloc (X->size2,X->size2);
  gsl_vector *beta = gsl_vector_calloc(dimX);

  gsl_multifit_linear(X,regressor,beta,cov,&chisq,work);
  gsl_blas_dgemv(CblasNoTrans,1,X,beta,0,ztmp);
  z = PearsonCorr(ztmp,regressor);
  r2 = CoD(regressor,ztmp);
  fprintf(stderr,"# Overall correlation of SI with target variable,  corr %.4f,  r2: %.4f\n\n",z,r2);

  gsl_vector_add_constant(regressor,mean);
  gsl_matrix_free(cov);
  gsl_vector_free(beta);
  gsl_multifit_linear_free(work);
  ZNorm(ztmp);
  return ztmp;
}


void XResiduals(gsl_matrix *X,int *test,gsl_vector *y,gsl_vector *residuals)
{
  if (X==NULL) {
    gsl_vector_set_zero(residuals);
    return;
  }
  size_t i,j;
  size_t ntest = y->size;
  double chisq=0;

  gsl_matrix *XX = gsl_matrix_calloc(ntest,X->size2);
  for (i=0; i<XX->size1; i++) {
    for (j=0; j<XX->size2; j++) {
      gsl_matrix_set(XX,i,j,gsl_matrix_get(X,test[i],j));
    }
  }
  MatColumnCentering(XX);
  gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(XX->size1,XX->size2);
  gsl_matrix *cov = gsl_matrix_alloc (XX->size2,XX->size2);
  gsl_vector *beta = gsl_vector_calloc(XX->size2);
  
  gsl_multifit_linear(XX,y,beta,cov,&chisq,work);
  gsl_multifit_linear_residuals(XX,y,beta,residuals);

  gsl_matrix_free(XX);
  gsl_matrix_free(cov);
  gsl_vector_free(beta);
  gsl_multifit_linear_free(work);
}



void xresiduals(gsl_vector *xtest,gsl_vector *ytest,gsl_vector *residuals)
{
  size_t i;
  double c0,c1,cov00,cov01,cov11,sumsq,y=0,y_err=0;
  
  gsl_fit_linear(xtest->data,1,ytest->data,1,xtest->size,&c0,&c1,&cov00,&cov01,&cov11,&sumsq);

  for (i=0; i<xtest->size; i++) {
    gsl_fit_linear_est(xtest->data[i],c0,c1,cov00,cov01,cov11,&y,&y_err);
    residuals->data[i] = ytest->data[i]-y;
  }
}
