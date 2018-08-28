/*
** prepare single subject analysis for vlisa
**
** G.Lohmann, Dec 2016
*/
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_spline.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define ABS(x) ((x) > 0 ? (x) : -(x))
#define SQR(x) ((x)*(x))

#define LEN     10000   /* buffer length        */
#define NTRIALS 10000   /* max number of trials */


/* standard parameter values for gamma function,Glover 99 */
double a1 = 6;
double b1 = 0.9;
double a2 = 12;
double b2 = 0.9;
double cc = 0.35;


typedef struct TrialStruct {
  int   id;
  float onset;
  float duration;
  float height;
} Trial;

void printmat(gsl_matrix *R,char *str)
{
  int i,j;
  fprintf(stderr," %s: \n",str);
  for (i=0; i<R->size1; i++) {
    for (j=0; j<R->size2; j++) {
      fprintf(stderr," %9.6f",gsl_matrix_get(R,i,j));
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");
}

void printvec(gsl_vector *x,char *str)
{
  int i;
  fprintf(stderr," %s: \n",str);
  for (i=0; i<x->size; i++) {
    fprintf(stderr," %f\n",x->data[i]);
  }
  fprintf(stderr,"\n");
}


/* output txt file for plotting */
void PlotDesign(gsl_matrix *X,double tr,VString filename)
{
  int i,j;
  FILE *fp = fopen(filename,"w");
  if (fp == NULL) VError("err opening plot file");
  double t = 0;
  for (i=0; i<X->size1; i++) {
    fprintf(fp," %f",t);
    for (j=0; j<X->size2; j++) {
      fprintf(fp," %f",gsl_matrix_get(X,i,j));
    }
    t += tr;
    fprintf(fp,"\n");
  }
  fclose(fp);
}


/* create design X */
gsl_matrix *VCreateDesign(int ntimesteps,int nevents,int hemomodel,gsl_matrix *covariates)
{
  int dim = nevents;
  if (hemomodel == 1) dim = (nevents-1)*2+1;
  if (hemomodel == 2) dim = (nevents-1)*3+1;
  if (covariates != NULL) dim += covariates->size2;  /* nuisance covariates without task labels */

  gsl_matrix *X = gsl_matrix_calloc(ntimesteps,dim);
  return X;
}


/* Glover kernel, gamma function */
double xgamma(double xx, double t0) 
{
  double x, y, scale = 20;
  double y1, y2;
  double d1, d2;
  x = xx - t0;
  if(x < 0 || x > 50)
    return 0;
  d1 = a1 * b1;
  d2 = a2 * b2;
  y1 = pow(x / d1, a1) * exp(-(x - d1) / b1);
  y2 = pow(x / d2, a2) * exp(-(x - d2) / b2);
  y = y1 - cc * y2;
  y /= scale;
  return y;
}


/* Glover kernel, gamma function, parameters changed for block designs */
double bgamma(double xx, double t0) 
{
  double x, y, scale = 120;
  double y1, y2;
  double d1, d2;
  double aa1 = 6;
  double bb1 = 0.9;
  double aa2 = 12;
  double bb2 = 0.9;
  double cx  = 0.1;
  x = xx - t0;
  if(x < 0 || x > 50)
    return 0;
  d1 = aa1 * bb1;
  d2 = aa2 * bb2;
  y1 = pow(x / d1, aa1) * exp(-(x - d1) / bb1);
  y2 = pow(x / d2, aa2) * exp(-(x - d2) / bb2);
  y = y1 - cx * y2;
  y /= scale;
  return y;
}


/* first derivative */
double deriv1_gamma(double x, double t0) 
{
  double d1, d2, y1, y2, y, xx;
  double scale = 20.0;
  xx = x - t0;
  if(xx < 0 || xx > 50) return 0;
  d1 = a1 * b1;
  d2 = a2 * b2;
  y1 = pow(d1, -a1) * a1 * pow(xx, (a1 - 1.0)) * exp(-(xx - d1) / b1)
    - (pow((xx / d1), a1) * exp(-(xx - d1) / b1)) / b1;
  y2 = pow(d2, -a2) * a2 * pow(xx, (a2 - 1.0)) * exp(-(xx - d2) / b2)
    - (pow((xx / d2), a2) * exp(-(xx - d2) / b2)) / b2;
  y = y1 - cc * y2;
  y /= scale;
  return y;
}


/* second derivative */
double deriv2_gamma(double x, double t0) 
{
  double d1, d2, y1, y2, y3, y4, y, xx;
  double scale = 20.0;
  xx = x - t0;
  if (xx < 0 || xx > 50) return 0;
  d1 = a1 * b1;
  d2 = a2 * b2;
  y1 = pow(d1, -a1) * a1 * (a1 - 1) * pow(xx, a1 - 2) * exp(-(xx - d1) / b1)
    - pow(d1, -a1) * a1 * pow(xx, (a1 - 1)) * exp(-(xx - d1) / b1) / b1;
  y2 = pow(d1, -a1) * a1 * pow(xx, a1 - 1) * exp(-(xx - d1) / b1) / b1
    - pow((xx / d1), a1) * exp(-(xx - d1) / b1) / (b1 * b1);
  y1 = y1 - y2;
  y3 = pow(d2, -a2) * a2 * (a2 - 1) * pow(xx, a2 - 2) * exp(-(xx - d2) / b2)
    - pow(d2, -a2) * a2 * pow(xx, (a2 - 1)) * exp(-(xx - d2) / b2) / b2;
  y4 = pow(d2, -a2) * a2 * pow(xx, a2 - 1) * exp(-(xx - d2) / b2) / b2
    - pow((xx / d2), a2) * exp(-(xx - d2) / b2) / (b2 * b2);
  y2 = y3 - y4;
  y = y1 - cc * y2;
  y /= scale;
  return y;
}



/* Gaussian function */
double xgauss(double xx, double t0) 
{
  double sigma=1.0;
  double x, y, z, a = 2.506628273;
  x = (xx - t0);
  z = x / sigma;
  y = exp(-0.5*z*z) / (sigma * a);
  return y;
}



void XConvolve(double *src,double *dst,double *kernel,int nt,int kernelsize)
{
  int i,j,jj,k;
  double sum=0;

  for (i=0; i<nt; i++) {

    sum = 0;
    k=0;
    for (j=i; j<i+kernelsize; j++) {
      jj = i-k;
      if (jj >= 0 && jj < nt) {
	sum += src[jj] * kernel[k];
      }
      k++;
    }
    dst[i] = sum;
  }
}



/* hemodynamic modelling */
void VHemoModel(Trial *trial,int ntrials,int nevents,int ntimesteps,double tr,int hemomodel,
		gsl_matrix *X,gsl_matrix *covariates)
{
  int i,j,k;
  double t,t0,t1,h;


  for(i = 0; i < nevents; i++) {
    double xmin = VRepnMaxValue(VFloatRepn);
    for(j = 0; j < ntrials; j++) {
      if(trial[j].id != i) continue;
      if(trial[j].duration < xmin) xmin = trial[j].duration;
    }
  }


  /* get kernels */
  t1 = 30.0;    /* kernel duration = 30 secs */
  int kernelsize = t1 / tr;
  if (tr < 0.05) VWarning(" implausible TR (%f seconds)",tr);


  double *kernel1=NULL,*kernel2=NULL;
  double *bkernel  = (double *)VCalloc(kernelsize,sizeof(double));
  double *kernel0  = (double *)VCalloc(kernelsize,sizeof(double));
  if (hemomodel == 1 || hemomodel == 2) {
    kernel1 = (double *)VCalloc(kernelsize,sizeof(double));
  }
  if (hemomodel == 2) {
    kernel2 = (double *)VCalloc(kernelsize,sizeof(double));
  }


  i = 0;
  for (t = 0; t < t1; t += tr) {
    if (i >= kernelsize) break;

    bkernel[i] = xgauss(t, 5.0*tr);
    kernel0[i] = xgamma(t, 0.0);

    if(hemomodel == 1 || hemomodel == 2)
      kernel1[i] = deriv1_gamma(t, 0);
    if(hemomodel == 2)
      kernel2[i] = deriv2_gamma(t, 0);
    i++;
  }



  /* tmp storage */
  double *x = (double *) VCalloc(ntimesteps,sizeof(double));
  double *y = (double *) VCalloc(ntimesteps,sizeof(double));


  /* constant in column 0 */
  gsl_matrix_set_zero(X);
  for (j=0; j<ntimesteps; j++) gsl_matrix_set(X,j,0,1.0);


  /* for each trial,event, do... */
  int col=1;
  for(i = 1; i < nevents; i++) {
    for (k=0; k<ntimesteps; k++) x[k] = y[k] = 0;

    /* read design info */
    int trialcount = 0;
    for (j = 0; j < ntrials; j++) {
      if(trial[j].id != i) continue;
      trialcount++;
      t0 = trial[j].onset;
      t1 = trial[j].onset + trial[j].duration;
      h  = trial[j].height;
      int k0 = (int) (t0/tr + 0.5);
      int k1 = (int) (t1/tr + 0.5);

      for (k=k0; k<=k1; k++) {
	if (k < 0 || k >= X->size1) continue;
	x[k] = h;
      }
      if(trialcount < 1)
	VError(" no trials in event %d, please re-number event-ids,  ntrials= %d", i + 1,ntrials);
    }


    /* convolve */
    if (hemomodel == 3)  {       /* block design, gaussian kernel */
      XConvolve(x,y,bkernel,ntimesteps,kernelsize);
    }
    else {                       /* gamma function kernel */
      XConvolve(x,y,kernel0,ntimesteps,kernelsize);
    }
    for (j=0; j<X->size1; j++) {
      gsl_matrix_set(X,j,col,y[j]);
    }
    col++;

    if (hemomodel == 1 || hemomodel == 2) {   /* gamma function kernel, first derivative */
      XConvolve(x,y,kernel1,ntimesteps,kernelsize);
      for (j=0; j<X->size1; j++) gsl_matrix_set(X,j,col,y[j]);
      col++;
    }
    if (hemomodel == 2) {   /* gamma function kernel, second derivative */
      XConvolve(x,y,kernel2,ntimesteps,kernelsize);
      for (j=0; j<X->size1; j++) gsl_matrix_set(X,j,col,y[j]);
      col++;
    }
  }


  /* add further covariates that have no task labels, will not be permuted */
  if (covariates != NULL) {
    for (i=0; i<covariates->size2; i++) {
      for (j=0; j<covariates->size1; j++) {
	gsl_matrix_set(X,j,col,gsl_matrix_get(covariates,j,i));
      }
      col++;
    }
  }


  /* cleanup */
  VFree(x);
  VFree(y);
  VFree(kernel0);
  VFree(bkernel);
  if (kernel1 != NULL) VFree(kernel1);
  if (kernel2 != NULL) VFree(kernel2);
}
