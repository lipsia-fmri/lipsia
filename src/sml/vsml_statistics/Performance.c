#include <viaio/Vlib.h>

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
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>

extern void RankNorm(gsl_vector *);
  
double RXResiduals(gsl_vector *xtest,gsl_vector *fitted)
{
  double c0,c1,cov00,cov01,cov11,sumsq,tss=0;
  gsl_fit_linear(xtest->data,1,fitted->data,1,xtest->size,&c0,&c1,&cov00,&cov01,&cov11,&sumsq);
  tss = gsl_stats_tss(fitted->data,1,fitted->size);
  double z = 1.0 - sumsq/tss;
  return z;
}


/* coefficient of determination */
double CoD(gsl_vector *observed,gsl_vector *fitted)
{
  size_t i,n=observed->size;   
  double *y = observed->data;
  double *f = fitted->data;
  
  double mean = gsl_stats_mean(y,1,n);
  double ss_total=0;
  for (i=0; i<n; i++) ss_total += (y[i]-mean)*(y[i]-mean);

  double ss_res=0;
  for (i=0; i<n; i++) ss_res += (y[i]-f[i])*(y[i]-f[i]);

  double r2 = 0;
  if (ss_total > 0) r2 = 1.0 - ss_res/ss_total;
  return r2;
}


/* weighted mean absolute percentage error */
double WMAPE(gsl_vector *observed,gsl_vector *fitted)
{
  size_t i,n=observed->size;
  double s1=0,s2=0,z=0;
  for (i=0; i<n; i++) {
    s1 += fabs(observed->data[i] - fitted->data[i]);
    s2 += fabs(observed->data[i]);
  }
  if (s2 > 0) z = s1/s2;
  return z;
}



/* Root Mean square error */
double RMSE(gsl_vector *observed,gsl_vector *fitted)
{
  size_t i,n=observed->size;
  double *y = observed->data;
  double *f = fitted->data;

  double s=0,nx=(double)n;
  for (i=0; i<n; i++){
    s += (y[i]-f[i])*(y[i]-f[i]);
  }
  s = sqrt(s/nx);
  return s;
}



/* Mean square error */
double MSE(gsl_vector *observed,gsl_vector *fitted)
{
  size_t i,n=observed->size;
  double *y = observed->data;
  double *f = fitted->data;

  double s=0,nx=(double)n;
  for (i=0; i<n; i++){
    s += (y[i]-f[i])*(y[i]-f[i]);
  }
  return s/nx;
}


/* Mean absolute error */
double MAE(gsl_vector *observed,gsl_vector *fitted)
{
  size_t i,n=observed->size;
  double *y = observed->data;
  double *f = fitted->data;

  double s=0,nx=(double)n;
  for (i=0; i<n; i++){
    s += fabs(y[i]-f[i]);
  }
  return s/nx;
}


/* Mean absolute error */
double XMAE(double *observed,double *fitted,size_t n)
{
  size_t i;
  double *y = observed;
  double *f = fitted;

  double s=0,nx=(double)n;
  for (i=0; i<n; i++){
    s += fabs(y[i]-f[i]);
  }
  return s/nx;
}


double SpearmanCorr(gsl_vector *ytest,gsl_vector *fitted)
{
  if (ytest==NULL || fitted==NULL) return 0;
  size_t ntest = ytest->size;
  double *work = (double *)VCalloc(2*ntest,sizeof(double));
  double spear = gsl_stats_spearman(ytest->data,1,fitted->data,1,ntest,work);
  VFree(work);
  return spear;
}

double PearsonCorr(gsl_vector *ytest,gsl_vector *fitted)
{
  if (ytest==NULL || fitted==NULL) return 0;
  return gsl_stats_correlation(ytest->data,1,fitted->data,1,ytest->size);
}

double Accuracy(gsl_vector *ytest,gsl_vector *fitted)
{
  size_t i;
  double tx=0,nx=(double)ytest->size;
  for (i=0; i<ytest->size; i++) {
    if (ytest->data[i] * fitted->data[i] > 0) tx++;
  }
  double acc = tx/nx;
  return acc;
}


/* explained variance */
double R2(gsl_vector *observed,gsl_vector *fitted)
{
  size_t i,n = observed->size;   
  double *y = observed->data;
  double *f = fitted->data;
  
  double ymean = gsl_stats_mean(y,1,n);

  double ssf=0;
  for (i=0; i<n; i++) ssf += (f[i]-ymean)*(f[i]-ymean);

  double sst=0;
  for (i=0; i<n; i++) sst += (y[i]-ymean)*(y[i]-ymean);

  double r2 = 0;
  if (sst > 0) r2 = ssf/sst;
  return r2;
}


double R2fit(gsl_vector *x,gsl_vector *y)
{
  double c0,c1,cov00,cov01,cov11,sumsq;

  gsl_fit_linear(x->data,1,y->data,1,x->size,&c0,&c1,&cov00,&cov01,&cov11,&sumsq);
  double tss = gsl_stats_tss(y->data,1,y->size);
  double r2 = 1.0-sumsq/tss;
  return r2;
}


void XResidualFit(gsl_vector *xtest,gsl_vector *ytest,gsl_vector *result,
		  double *corr,double *r2,double *res_corr,double *res_r2)
{
  size_t i,ntest=xtest->size;
  double c0,c1,cov00,cov01,cov11,chisq,y=0,tss=0;

  gsl_vector *ztmp = gsl_vector_calloc(ntest);
  gsl_vector *residuals = gsl_vector_calloc(ntest);
  (*corr)=(*r2)=(*res_corr)=(*res_r2)=0;

  /* how much of the variance in estimated IQ ('result') is explained by education ('xtest')? */
  gsl_fit_linear(xtest->data,1,result->data,1,xtest->size,&c0,&c1,&cov00,&cov01,&cov11,&chisq);
  
  for (i=0; i<xtest->size; i++) {
    y = c0 + c1*xtest->data[i];
    ztmp->data[i] = y;
    residuals->data[i] = result->data[i]-y;
  }
  tss = gsl_stats_tss(result->data,1,result->size);
  *r2 = (1.0-chisq/tss);
  *corr = PearsonCorr(result,ztmp);
  gsl_vector_free(ztmp);
  if (ytest == NULL) return;

  
  /* how much of the unexplained variance ('residuals') is explained by the true IQ ('ytest')? */
  gsl_fit_linear(ytest->data,1,residuals->data,1,ntest,&c0,&c1,&cov00,&cov01,&cov11,&chisq);
  tss = gsl_stats_tss(residuals->data,1,ntest);
  
  
  /* how much of the true IQ ('ytest') is explained by the residuals ('residuals')? */
  /*
  gsl_fit_linear(residuals->data,1,ytest->data,1,ntest,&c0,&c1,&cov00,&cov01,&cov11,&chisq);
  tss = gsl_stats_tss(ytest->data,1,ntest);
  */
  *res_r2 = (1.0-chisq/tss);
  *res_corr = PearsonCorr(ytest,residuals);
  gsl_vector_free(residuals);
}


double Resid(gsl_vector *ytest,gsl_vector *result,gsl_vector *xtest)
{
  size_t i;
  gsl_vector *residuals = gsl_vector_calloc(ytest->size);
  
  for (i=0; i<ytest->size; i++) {
    residuals->data[i] = ytest->data[i]-result->data[i];
  }
  double z = PearsonCorr(xtest,residuals);
  
  return z;
}


