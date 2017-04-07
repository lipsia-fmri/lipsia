/*
** compute histogram density function
**
** G.Lohmann, Nov 2012
*/

#include <viaio/Vlib.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_math.h>

#define SQR(x) ((x)*(x))
#define ABS(x) ((x) > 0 ? (x) : -(x))


/* kernel density estimation */
double GaussKernel(double x)
{
  return exp(-0.5*x*x)/sqrt(2.0*M_PI);
}

double VKernelDensity(double x,double h,double nh,gsl_histogram *histogram)
{
  size_t i;
  double sum,upper,lower;

  sum = 0;
  for (i=0; i<gsl_histogram_bins(histogram); i++) {
    double nx = gsl_histogram_get (histogram,i);
    gsl_histogram_get_range (histogram,i,&lower,&upper);
    double xi = lower + (upper-lower)*0.5;
    double u = (x-xi)/h;
    sum += nx * GaussKernel(u);
  }
  return sum/nh;
}

void VPrintHistogram(gsl_histogram *histogram,int numperm,VString filename)
{
  size_t i;
  double z=0,mx=0,lower=0,upper=0,fz=0;
  extern double GaussKernel(double x);
  extern double VKernelDensity(double,double,double,gsl_histogram *);

  gsl_histogram_get_range (histogram,0,&lower,&upper);
  double alpha = 1.06;  /* silverman's rule */
  double sig = gsl_histogram_sigma(histogram);
  double kx = gsl_histogram_sum(histogram);
  double h  = alpha * sig * pow(kx,-0.2);
  double nh = h*kx;

  FILE *fph = fopen(filename,"w");
  if (!fph) VError(" err opening hist file");
  for (i=0; i<gsl_histogram_bins(histogram); i++) {
    gsl_histogram_get_range (histogram,i,&lower,&upper);  
    z  = lower;
    fz = VKernelDensity(z,h,nh,histogram);
    mx = gsl_histogram_get (histogram,i);
    fprintf(fph,"%lf %lf %.2lf\n",lower,fz,mx);
  }
  fclose(fph);
}


/* update histogram */
void HistoUpdate(float *A,size_t nvox,gsl_histogram *hist)
{
  double tiny = 1.0e-8;
  size_t i,j;  
  double hmin = gsl_histogram_min(hist);
  double hmax = gsl_histogram_max(hist);

  for (i=0; i<nvox; i++) {
    for (j=0; j<i; j++) {
      const size_t k=j+i*(i+1)/2;
      double u = (double)A[k];
      if (ABS(u) < tiny) continue;     
      if (u < hmin) u = hmin;
      if (u > hmax-tiny) u = hmax-tiny;
      gsl_histogram_increment (hist,u);
    }
  }
}
