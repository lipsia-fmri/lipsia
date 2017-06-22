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

#define SQR(x) ((x)*(x))
#define ABS(x) ((x) > 0 ? (x) : -(x))

extern float kth_smallest(float *a, size_t n, size_t k);
#define Median(a,n) kth_smallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))


/* scale z-values */
void VZScale(VImage src,float stddev)
{
  size_t i=0;
  float u=0,tiny=1.0e-8;
  VFloat *pp=VImageData(src);
  for (i=0; i<VImageNPixels(src); i++) {
    u = (*pp);
    if (fabs(u) > tiny) (*pp) = u/stddev;
    pp++;
  }  
}


/* scaling and centering of z-values */
void VZNormalize(VImage src,float stddev)
{
  size_t i=0;
  float u=0,s1=0,nx=0,tiny=1.0e-8;

  VFloat *pp=VImageData(src);
  for (i=0; i<VImageNPixels(src); i++) {
    u = (*pp);
    if (fabs(u) > tiny) {
      s1 += u;
      nx++;
    }
    pp++;
  }
  if (nx < 3.0) VError(" nx: %f\n",nx);
  float mean = s1/nx;

  pp=VImageData(src);
  for (i=0; i<VImageNPixels(src); i++) {
    u = (*pp);
    if (fabs(u) > tiny) (*pp) = (u-mean)/stddev;
    pp++;
  }
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
  /* count number of nonzero voxels */
  float u=0,v=0,tiny=1.0e-8;
  size_t i=0,n=0;
  VFloat *pp=VImageData(src);
  for (i=0; i<VImageNPixels(src); i++) {
    u = (*pp++);
    if (fabs(u) < tiny) continue;
    n++;
  }
  fprintf(stderr," number of nonzero voxels: %lu\n",n);


  /* get every second nonzero data point */
  size_t nn = n/2;
  float *data = (float *) VCalloc(nn,sizeof(float));
  n=0;
  pp=VImageData(src);
  for (i=0; i<VImageNPixels(src); i+=2) {
    u = (*pp++);
    v = (*pp++);
    if (fabs(u) < tiny || fabs(v) < tiny) continue;
    if (n >= nn) break;
    data[n] = u;
    n++;
  }
  if (n < 10) VError(" not enough nonzero voxels: %lu",n);


  /* get upper and lower quantiles */
  size_t k0 = (size_t)(0.001 * (float)n);
  if (k0 < 3) k0 = 3;
  size_t k1 = (size_t)(0.999 * (float)n);
  if (k1 < n-3) k1 = n-3;
  if (k0 >= k1) VError(" not enough nonzero voxels: %lu",n);
  qsort(data,n,sizeof(float),compare_function);
  float u0 = data[k0];
  float u1 = data[k1];
  float eps = (fabs(u1-u0))*0.05;
  VFree(data);

  /* output */
  *hmin = (double)(u0-eps);
  *hmax = (double)(u1+eps);
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
