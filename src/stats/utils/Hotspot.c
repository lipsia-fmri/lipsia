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


/* get image statistics */
void ImageStats(VImage src,double *ave,double *var,double *hmin,double *hmax)
{
  size_t i=0;
  double umin = 99999.9;
  double umax = -99999.9;
  double u=0,s1=0,s2=0,nx=0,mean=0,tiny=1.0e-8;

  VFloat *pp=VImageData(src);
  for (i=0; i<VImageNPixels(src); i++) {
    u = (double)(*pp++);
    if (u < umin) umin = u;
    if (u > umax) umax = u;
    if (ABS(u) < tiny) continue;
    s1 += u;
    s2 += u*u;
    nx++;
  }
  *hmin = umin;
  *hmax = umax;  
  mean = s1/nx;
  *ave = mean;
  *var = (s2 - nx * mean * mean) / (nx - 1.0);
}


/* get image variance */
double VImageVar(VImage src)
{
  size_t i=0;
  double u=0,s1=0,s2=0,nx=0,mean=0,tiny=1.0e-8;
  VFloat *pp=VImageData(src);
  for (i=0; i<VImageNPixels(src); i++) {
    u = (double)(*pp++);
    if (ABS(u) < tiny) continue;
    s1 += u;
    s2 += u*u;
    nx++;
  }
  if (nx < 3.0) VError(" nx: %f\n",nx);
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



/* median filter in 6-adj neighbourhood */
void VMedian6(VImage src)
{
  int b,r,c,bb,rr,cc,m,k,l;
  float u=0,x=0,z=0,tiny=1.0e-10;
  float data[8];

  int nslices = VImageNBands(src);
  int nrows = VImageNRows(src);
  int ncols = VImageNColumns(src);

  VImage dst = VCreateImageLike(src);  
  VCopyImagePixels(src,dst,VAllBands);

  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {

	x = VPixel(src,b,r,c,VFloat);
	if (fabs(x) < tiny) continue;  /* outside of brain */

	size_t n=0;
	data[n++] = x;

	for (m=-1; m<=1; m++) {
	  bb = b+m;
	  if (bb < 0 || bb >= nslices) continue;

	  for (k=-1; k<=1; k++) {
	    rr = r+k;
	    if (rr < 0 || rr >= nrows) continue;

	    for (l=-1; l<=1; l++) {
	      cc = c+l;
	      if (cc < 0 || cc >= ncols) continue;

	      int d = (m*m + k*k + l*l);
	      if (d >= 2) continue;
	      u = VPixel(src,bb,rr,cc,VFloat);
	      if (fabs(u) < tiny) continue;  /* outside of brain */

	      data[n] = u;
	      n++;
	    }
	  }
	}
	if (n < 3) continue;
	z = Median(data,n);
	VPixel(dst,b,r,c,VFloat) = z;
      }
    }
  }
  VCopyImagePixels(dst,src,VAllBands);
  VDestroyImage(dst);
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
