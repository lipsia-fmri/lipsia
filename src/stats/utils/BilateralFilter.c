/*
** bilateral filter
**
** G.Lohmann, Jan 2017
*/
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

extern float kth_smallest(float *a, size_t n, size_t k);
#define Median(a,n) kth_smallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))
#define ABS(x) ((x) > 0 ? (x) : -(x))


/* Median filter in 18-adj neighbourhood */
float XMedian18(VImage src,int b,int r,int c)
{ 
  int bb,rr,cc,m,k,l;
  float u=0,z=0,tiny=1.0e-6;
  float data[30];

  size_t n=0;
  for (m=-1; m<=1; m++) {
    bb = b+m;
    for (k=-1; k<=1; k++) {
      rr = r+k;
      for (l=-1; l<=1; l++) {
	cc = c+l;
	if (ABS(m) > 0 && ABS(k) > 0 && ABS(l) > 0) continue; 
	u = VPixel(src,bb,rr,cc,VFloat);
	if (fabs(u) < tiny) continue;
	data[n] = u;
	n++;
      }
    }
  }
  z = 0;
  if (n > 9) z = Median(data,n);
  return z;
}


/* Median filter in 6-adj neighbourhood */
float XMedian6(VImage src,int b,int r,int c)
{ 
  int bb,rr,cc,m,k,l;
  float u=0,z=0,tiny=1.0e-6;
  float data[8];

  size_t n=0;
  u = VPixel(src,b,r,c,VFloat);
  data[n++] = u;

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
  z = 0;
  if (n > 3) z = Median(data,n);
  return z;
}



/* https://codingforspeed.com/using-faster-exponential-approximation/ */
double fastexp(double x) 
{
  x = 1.0 + x / 2048.0;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x;
  return x;
}


/* bilateral filter */
void VBilateralFilter(VImage src,VImage dest,int radius,double var1,double var2,int numiter)
{
  float nx=0,mx=0,z=0;
  double u=0,x=0,w=0,d=0,s1=0,s2=0,tiny=1.0e-10;
  int b,r,c,bb,rr,cc,m,k,l,iter;

  /* no filtering */
  if (numiter==0) {
    VCopyImagePixels(src,dest,VAllBands);
    return;
  }


  /* ini dest image */
  int nbands = VImageNBands(src);
  int nrows  = VImageNRows(src);
  int ncols  = VImageNColumns(src);
  VFillImage(dest,VAllBands,0);



  /* get max neighbourhood size */
  nx=0;
  int wn = radius-1;
  for (m=-radius; m<=radius; m++) {
    for (k=-radius; k<=radius; k++) {
      for (l=-radius; l<=radius; l++) {
	if ((ABS(m) > wn && ABS(k) > wn && ABS(l) > wn)) continue;
	nx++;
      }
    }
  }


  /* loop through voxels */
  for (iter=0; iter<numiter; iter++) {

    for (b=radius; b<nbands-radius; b++) {
      for (r=radius; r<nrows-radius; r++) {
	for (c=radius; c<ncols-radius; c++) {

	  x = (double)VPixel(src,b,r,c,VFloat);
	  if (fabs(x) < tiny) continue;

	  mx = 0;
	  s1 = s2 = 0;
	  for (m=-radius; m<=radius; m++) {
	    bb = b+m;
	    for (k=-radius; k<=radius; k++) {
	      rr = r+k;
	      for (l=-radius; l<=radius; l++) {
		cc = c+l;
		if ((ABS(m) > wn && ABS(k) > wn && ABS(l) > wn)) continue;

		u = (double)VPixel(src,bb,rr,cc,VFloat);
		if (fabs(u) < tiny) continue;

		d = (double)(m*m + k*k + l*l);
		z = (x-u)*(x-u)/var1 + d/var2;
		/* w = exp(-z); */
		w = fastexp(-z);
		s1 += u*w;
		s2 += w;
		mx++;
	      }
	    }
	  }
	  z=0;

	  /* bilateral filter if local neighbourhood mostly inside the brain */
	  if ((s2 > 0) && (mx/nx > 0.5)) {
	    z = (float)(s1/s2);
	  }

	  /* median filter if local neighbourhood partly outside of brain */
	  else {  
	    z = XMedian18(src,b,r,c);
	  }
	  VPixel(dest,b,r,c,VFloat) = z;
	}
      }
    }
    if (iter < numiter-1 && numiter > 1) src = VCopyImagePixels(dest,src,VAllBands);
  }
}
