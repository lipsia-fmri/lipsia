/* 
** Various utilities for LISA algorithm
**
** G.Lohmann, Jan 2017
*/
#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_histogram.h>


#define SQR(x) ((x)*(x))
#define ABS(x) ((x) > 0 ? (x) : -(x))

extern float kth_smallest(float *a, size_t n, size_t k);
#define Median(a,n) kth_smallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))



/* update histogram */
void HistoUpdate(VImage src1,gsl_histogram *hist)
{
  float u,tiny = 1.0e-6;
  size_t i;
  float xmin = gsl_histogram_min (hist);
  float xmax = gsl_histogram_max (hist);

  float *pp1 = VImageData(src1);
  for (i=0; i<VImageNPixels(src1); i++) {
    u = *pp1++;
    if (ABS(u) < tiny) continue;
    if (u > xmax) u = xmax-tiny;
    if (u < xmin) u = xmin+tiny;
    gsl_histogram_increment (hist,u);
  }
}



/* In moderately skewed or asymmetrical distribution (Pearson) */
float VGetMode(VImage src)
{
  size_t i,n=0;
  float u,tiny=1.0e-8;
  VFloat *pp = VImageData(src);
  for (i=0; i<VImageNPixels(src); i++) {
    u = (*pp++);
    if (fabs(u) > tiny) n++;
  }
  float *data = (float *) VCalloc(n,sizeof(float));
  
  double sum=0,nx=0;
  n=0;
  pp = VImageData(src);
  for (i=0; i<VImageNPixels(src); i++) {
    u = (*pp++);
    if (fabs(u) > tiny) { 
      data[n] = u;
      sum += u;
      n++;
    }
  }
  nx = (double)n;
  float median = Median(data,n);
  float mean = (float)(sum/nx);
  float mode = 3.0*median - 2.0*mean;
  /* fprintf(stderr," mean: %f,  median: %f,  mode: %f\n",mean,median,mode); */
  VFree(data);
  return mode;
}



/* scale z-values */
void VZScale(VImage src,float mode,float stddev)
{
  size_t i=0;
  float u=0,tiny=1.0e-8;
  VFloat *pp=VImageData(src);
  for (i=0; i<VImageNPixels(src); i++) {
    u = (*pp);
    if (fabs(u) > tiny) (*pp) = (u-mode)/stddev;
    pp++;
  }  
}

/* get max z-value */
float VZMax(VImage src)
{
  size_t i=0;
  float u=0,umax=0;
  VFloat *pp=VImageData(src);
  for (i=0; i<VImageNPixels(src); i++) {
    u = (*pp++);
    if (u > umax) umax = u;
  }
  return umax;
}



/* get image variance */
double VImageVar(VImage src)
{
  size_t i=0;
  double u=0,s1=0,s2=0,nx=0,tiny=1.0e-8;
  VFloat *pp=VImageData(src);
  for (i=0; i<VImageNPixels(src); i++) {
    u = (double)(*pp++);
    if (ABS(u) < tiny) continue;
    s1 += u;
    s2 += u*u;
    nx++;
  }
  if (nx < 3.0) VError(" number of nonzero voxels: %g\n",nx);
  double mean = s1/nx;
  double var = (s2 - nx * mean * mean) / (nx - 1.0);
  return var;
}


/* count nonzero voxels */
void VImageCount(VImage src)
{
  float u=0;
  size_t i=0,npos=0;
  VFloat *pp=VImageData(src);
  for (i=0; i<VImageNPixels(src); i++) {
    u = (float)(*pp++);
    if (u > 0) npos++;
  }
}


/* check histogram range, image is normalized so [-10,10] should be okay */
void VGetHistRange(VImage src,double *hmin,double *hmax)
{
  size_t i=0;
  double u=0,x=0,y=0;
  double zmin = -10.0;
  double zmax = 10.0;

  VFloat *pp=VImageData(src);
  size_t npix=0,npos=0,nneg=0;
  for (i=0; i<VImageNPixels(src); i++) {
    u = (*pp++);
    if (u < zmin) nneg++;
    if (u > zmax) npos++;
    if (fabs(u) > 0) npix++;
  }
  if (npix < 1) VError(" image has no non-zero pixels");

  if (npos > 0 || nneg > 0) {
    x = (double)npos/(double)npix;
    y = (double)nneg/(double)npix;
    if (x > 0.01 || y > 0.01) {
      zmin = -20.0;
      zmax = 20.0;
    }
  }
 
  /* output */
  *hmin = zmin;
  *hmax = zmax;
}





/* remove isolated voxels */
void VIsolatedVoxels(VImage src,float threshold)
{
  int b,r,c,bb,rr,cc;
  float u=0;
  size_t i;

  VImage tmp = VCreateImage(VImageNBands(src),VImageNRows(src),VImageNColumns(src),VBitRepn);
  VFillImage(tmp,VAllBands,0);

  VFloat *pp=VImageData(src);
  VBit *pa=VImageData(tmp);
  for (i=0; i<VImageNPixels(src); i++) {
    u = (*pp++);
    if (u > threshold) *pa = 1;
    pa++;
  }

  for (b=1; b<VImageNBands(src)-1; b++) {
    for (r=1; r<VImageNRows(src)-1; r++) {
      for (c=1; c<VImageNColumns(src)-1; c++) {

	int n=0;
	for (bb=b-1; bb<=b+1; bb++) {
	  for (rr=r-1; rr<=r+1; rr++) {
	    for (cc=c-1; cc<=c+1; cc++) {
	      if (VPixel(tmp,bb,rr,cc,VBit) > 0) n++;
	      if (n >= 2) goto skip;
	    }
	  }
	}
	if (n < 2) VPixel(src,b,r,c,VFloat) = 0;
      skip: ;
      }
    }
  }
  VDestroyImage(tmp);
}
