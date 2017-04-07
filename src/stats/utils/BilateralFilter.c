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


/* https://codingforspeed.com/using-faster-exponential-approximation/ */
double fastexp(double x) 
{
  x = 1.0 + x / 2048.0;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x;
  return x;
}

void VBilateralFilter(VImage src,VImage dest,VImage tmp,int radius,double var1,double var2,int numiter)
{
  float nx,n;
  double u=0,x=0,w=0,d=0,s1=0,s2=0,z=0,tiny=1.0e-6;
  int b,r,c,bb,rr,cc,m,k,l;
  int iter;

  int nbands = VImageNBands(src);
  int nrows  = VImageNRows(src);
  int ncols  = VImageNColumns(src);
  nx = (float)((2*radius+1)*(2*radius+1)*(2*radius+1));

  VFillImage(dest,VAllBands,0);
  if (numiter > 1) tmp = VCopyImagePixels(src,tmp,VAllBands);
  else tmp = src;

  for (iter=0; iter<numiter; iter++) {

    for (b=radius; b<nbands-radius; b++) {
      for (r=radius; r<nrows-radius; r++) {
	for (c=radius; c<ncols-radius; c++) {

	  x = (double)VPixel(tmp,b,r,c,VFloat);
	  if (fabs(x) < tiny) continue;

	  n = 0;
	  s1 = s2 = 0;
	  for (m=-radius; m<=radius; m++) {
	    bb = b+m;
	    for (k=-radius; k<=radius; k++) {
	      rr = r+k;
	      for (l=-radius; l<=radius; l++) {
		cc = c+l;
		u = (double)VPixel(src,bb,rr,cc,VFloat);
		if (fabs(u) < tiny) continue;

		d = (double)(m*m + k*k + l*l);
		z = (x-u)*(x-u)/var1 + d/var2;
		/* w = exp(-z); */
		w = fastexp(-z);
		s1 += u*w;
		s2 += w;
		n++;
	      }
	    }
	  }
	  z=0;
	  /* if ((s2 > 0) && (n/nx > 0.5)) z = s1/s2;  !!!! */
	  /* if ((s2 > 0) && (n/nx > 0.75)) z = s1/s2;*/
	  if ((s2 > 0) && (n/nx > 0.6)) z = s1/s2;
	  VPixel(dest,b,r,c,VFloat) = (float)z;
	}
      }
    }
    if (iter < numiter-1 && numiter > 1) tmp = VCopyImagePixels(dest,tmp,VAllBands);
  }
}
