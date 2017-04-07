/* 
** find hotspots in zmaps
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

extern void VBilateralFilter(VImage src,VImage dest,VImage,int radius,double var1,double var2,int);


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


/* median filter in 6-adj neighbourhood */
void VMedian6(VImage src,VImage dst)
{
  int b,r,c,bb,rr,cc,m,k,l;
  float u=0,x=0,z=0,tiny=1.0e-10;
  float data[7];

  VFillImage(dst,VAllBands,0);
  for (b=1; b<VImageNBands(src)-1; b++) {
    for (r=1; r<VImageNRows(src)-1; r++) {
      for (c=1; c<VImageNColumns(src)-1; c++) {

	x = VPixel(src,b,r,c,VFloat);
	if (fabs(x) < tiny) continue;  /* outside of brain */

	size_t n=0;
	for (m=-1; m<=1; m++) {
	  bb = b+m;
	  for (k=-1; k<=1; k++) {
	    rr = r+k;
	    for (l=-1; l<=1; l++) {
	      cc = c+l;
	      int d = (m*m + k*k + l*l);
	      if (d >= 2) continue;
	      u = VPixel(src,bb,rr,cc,VFloat);
	      if (fabs(u) < tiny) continue;  /* outside of brain */
	      data[n] = u;
	      n++;
	    }
	  }
	}
	if (n < 4) continue;
	z = Median(data,n);
	if (x < z) z = x;  /* minimum of median and current value */
	VPixel(dst,b,r,c,VFloat) = z;
      }
    }
  }
}


/* remove isolated voxels */
void VIsolatedVoxels(VImage src,VImage dst,float threshold)
{
  int b,r,c,bb,rr,cc;  
  for (b=1; b<VImageNBands(src)-1; b++) {
    for (r=1; r<VImageNRows(src)-1; r++) {
      for (c=1; c<VImageNColumns(src)-1; c++) {
	if (VPixel(src,b,r,c,VFloat) < threshold) {
	  VPixel(dst,b,r,c,VFloat) = 0;
	  continue;
	}

	int n=0;
	for (bb=b-1; bb<=b+1; bb++) {
	  for (rr=r-1; rr<=r+1; rr++) {
	    for (cc=c-1; cc<=c+1; cc++) {
	      if (VPixel(src,bb,rr,cc,VFloat) >= threshold) n++;
	      if (n >= 2) goto skip;
	    }
	  }
	}
	if (n < 2) VPixel(dst,b,r,c,VFloat) = 0;
      skip: ;
      }
    }
  }
}

/* hotspot analysis */
void Hotspot(VImage src,VImage dst,VImage tmp,int type,int radius,double var1,double var2,int numiter)
{
  VFillImage(dst,VAllBands,0);
  VFillImage(tmp,VAllBands,0);

  switch (type) {
  case 0:
    VBilateralFilter(src,dst,tmp,radius,var1,var2,numiter);
    break;
  case 1:
    VBilateralFilter(src,tmp,dst,radius,var1,var2,numiter);
    VMedian6(tmp,dst);   /* cleanup */
    break;
  default:
    VError(" illegal type ");
  }
}

