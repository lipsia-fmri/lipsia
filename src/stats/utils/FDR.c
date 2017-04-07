#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <viaio/option.h>
#include <via/via.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_math.h>

#define SQR(x) ((x)*(x))
#define ABS(x) ((x) > 0 ? (x) : -(x))


extern void VCheckImage(VImage src);


void HistParams(gsl_histogram *histogram,double *hx,double *nhx)
{
  double alpha = 1.06;  /* silverman's rule */
  double sig = gsl_histogram_sigma(histogram);
  double kx = gsl_histogram_sum(histogram);
  double h  = alpha * sig * pow(kx,-0.2);
  double nh = h*kx;
  *hx = h;
  *nhx = nh;
}

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


double CFDR(size_t j,gsl_histogram_pdf *cdf0,gsl_histogram_pdf *cdfz,double p0)
{
  double F0 = 1.0-cdf0->sum[j];
  double Fz = 1.0-cdfz->sum[j];
  double Fdr = 1.0;
  if (Fz  > 0.0) Fdr = p0*F0/Fz;
  if (Fdr > 1.0) Fdr = 1.0;
  if (Fdr < 0.0) Fdr = 0.0;
  return Fdr;
}

/* get prior p0 */
double GetP0(gsl_histogram *nullhist,gsl_histogram *realhist)
{
  size_t i;
  double p0=1.0,h0=0,hz=0,d=0,sum=0,lower,upper;

  sum = gsl_histogram_sum(nullhist);
  gsl_histogram_scale(nullhist,1.0/sum);
  sum = gsl_histogram_sum(realhist);
  gsl_histogram_scale(realhist,1.0/sum);
  
  sum = 0;
  for (i=0; i<gsl_histogram_bins(nullhist); i++) {
    gsl_histogram_get_range (nullhist,i,&lower,&upper);
    /* if (upper < 0) continue; */
    h0 = gsl_histogram_get(nullhist,i);
    hz = gsl_histogram_get(realhist,i);
    d = hz-h0;
    if (d > 0) sum += d;
  }
  p0 = 1.0 - sum;
  return p0;
}


/* tail-area FDR */
double GetCutoff(gsl_histogram *realhist,gsl_histogram_pdf *cdf0,gsl_histogram_pdf *cdfz,double p0,double alpha)
{
  size_t i=0;
  size_t nbins = gsl_histogram_bins (realhist);
  int flag1 = 0;
  double lower=0,upper=0;
  double cutoff=VRepnMaxValue(VFloatRepn);

  for (i=0; i<nbins; i++) {
    if (gsl_histogram_get(realhist,i) > 0) {
      gsl_histogram_get_range (realhist,i,&lower,&upper);
      double z = lower + 0.5*(upper-lower); 
      double Fdr = CFDR(i,cdf0,cdfz,p0);
      if (Fdr < alpha && flag1 == 0 && z > 0) {
	flag1 = 1;
	cutoff = z;
      }
    }
  }
  return cutoff;
}


void PrintHistogram(gsl_histogram *nullhist,gsl_histogram *realhist,gsl_histogram_pdf *cdf0,gsl_histogram_pdf *cdfz,VString filename)
{
  size_t i=0;
  size_t nbins = gsl_histogram_bins (nullhist);
  double hr,h0,nhr,nh0;
  double lower=0,upper=0;
  HistParams(realhist,&hr,&nhr);
  HistParams(nullhist,&h0,&nh0);

  FILE *fp = fopen(filename,"w");
  if (!fp) VError(" err opening file %s",filename);
  fprintf(fp,"#       z            f0            fz            F0             Fz           FDR\n");
  fprintf(fp,"#----------------------------------------------------------------------------------\n");

  for (i=0; i<nbins; i++) {
    if (gsl_histogram_get(realhist,i) > 0) {
      gsl_histogram_get_range (realhist,i,&lower,&upper);
      double z = lower + 0.5*(upper-lower); 
      double fz = VKernelDensity(z,hr,nhr,realhist);
      double f0 = VKernelDensity(z,h0,nh0,nullhist);

      double Fz = 1.0-cdfz->sum[i];
      double F0 = 1.0-cdf0->sum[i];
      double Fdr = 1.0;
      if (Fz  > 1.0e-12) Fdr = F0/Fz;
      if (Fdr > 1.0) Fdr = 1.0;
      if (Fdr < 0.0) Fdr = 0.0;
      
      fprintf(fp," %12.8lf  %12.8lf  %12.8lf  %12.8lf  %12.8lf  %12.8lf\n",z,f0,fz,F0,Fz,Fdr);
    }
  }
  fclose(fp);
}


void FDR(VImage src,VImage dest,double alpha,gsl_histogram *nullhist,gsl_histogram *realhist,VString filename)
{
  size_t i,j;
  size_t nbins = gsl_histogram_bins (nullhist);

  /* cumulative distribution functions (CDF) */
  gsl_histogram_pdf *cdfz = gsl_histogram_pdf_alloc(nbins);
  gsl_histogram_pdf *cdf0 = gsl_histogram_pdf_alloc(nbins);
  gsl_histogram_pdf_init (cdfz,realhist);
  gsl_histogram_pdf_init (cdf0,nullhist);

  if (strlen(filename) > 2) {
    PrintHistogram(nullhist,realhist,cdf0,cdfz,filename);
  }

  double p0 = 1.0;
  p0 = GetP0(nullhist,realhist);
  fprintf(stderr,"\n p0= %f\n",p0);
  double cutoff = GetCutoff(realhist,cdf0,cdfz,p0,alpha);
  


  /* apply FDR */
  size_t ic=0,n=0;
  double u,tiny=1.0e-8;
  VFillImage(dest,VAllBands,0);
  VFloat *pp0 = VImageData(src);
  VFloat *pp1 = VImageData(dest);
  for (i=0; i<VImageNPixels(src); i++) {
    u = *pp0;

    if (ABS(u) > tiny) {
      if (u < gsl_histogram_min(realhist)) u = gsl_histogram_min(realhist)+0.0001;
      if (u > gsl_histogram_max(realhist)) u = gsl_histogram_max(realhist)-0.0001;
      if (gsl_histogram_find (realhist,(double)u,&j) == GSL_SUCCESS) {
	double xFdr = CFDR(j,cdf0,cdfz,p0);

	*pp1 = 0;
	if (alpha < 0.5) {
	  if (u >= cutoff) {
	    *pp1=(1.0-xFdr);
	    ic++;
	  }
	}
	else {
	  *pp1 = (1.0-xFdr);
	}
	n++;
      }
    }
    pp0++;
    pp1++;
  }
  if (alpha < 0.99) fprintf(stderr," significant voxels:  %lu  of  %lu,  cutoff: %f\n",ic,n, cutoff);
  return;
}
