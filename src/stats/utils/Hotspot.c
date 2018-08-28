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
  if (nx < 3.0) VError(" nx: %f\n",nx);
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
  fprintf(stderr," positive voxels:  %lu\n",npos);
}


/* needed for qsort */
int compare_function(const void *a,const void *b) 
{
  float *x = (float *) a;
  float *y = (float *) b;
  float d = ((*x) - (*y));
  if (d > 0) return 1;
  else if (d < 0) return -1;
  else return 0;
}


/* get histogram range */
void VGetHistRange(VImage src,double *hmin,double *hmax)
{
  size_t i=0;
  double zmin = VRepnMaxValue(VFloatRepn);
  double zmax = VRepnMinValue(VFloatRepn);
  double u=0;

  VFloat *pp=VImageData(src);
  for (i=0; i<VImageNPixels(src); i++) {
    u = (*pp++);
    if (u < zmin) zmin = u;
    if (u > zmax) zmax = u;
  }
  fprintf(stderr," zmin,zmax: %f %f\n",zmin,zmax);

  /* correct implausible ranges */
  if (zmin < -20.0) zmin = -20.0;
  if (zmax > 20.0) zmax = 20.0;
  
  if (zmin > -10.0) zmin = -10.0;
  if (zmax < 10.0) zmax = 10.0;


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
