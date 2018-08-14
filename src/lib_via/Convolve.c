/*! \file
  Convolution

This file contains code for 2D and 3D convolutions, and for
1D convolutions. 1D convolutions can be used to implement
separable 2D/3D filters.

\par Author:
Gabriele Lohmann, MPI-CBS
*/


/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/mu.h>

/* From the standard C libaray: */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>



/*!
\fn VImage VConvolveCol (VImage src,VImage dest,VImage kernel)
\brief 1D convolution in column-direction (for separable filters)
\param src    input image  (any repn)
\param dest   output image (float repn)
\param kernel raster image containing convolution kernel (float repn)
*/
VImage VConvolveCol(VImage src,VImage dest,float *kernel,int dim)
{
  int b,r,c,nbands,nrows,ncols;
  int c0,c1,cc,k=0;
  float sum,x;

  if (dim%2 == 0) VError("VConvolveCol: kernel dim must be an odd number, dim= %d",dim);
  int d = dim/2;

  if (VPixelRepn(src) != VFloatRepn) VError(" input pixel repn must be float");

  nrows  = VImageNRows (src);
  ncols  = VImageNColumns (src);
  nbands = VImageNBands (src);

  dest = VSelectDestImage("VConvolveCol",dest,nbands,nrows,ncols,VFloatRepn);
  VFillImage(dest,VAllBands,0);

  for (b=0; b<nbands; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	sum = 0;
	c0 = c-d;
	c1 = c+d;
	k=0;
	for (cc=c0; cc<=c1; cc++) {
	  if (cc >= 0 && cc < ncols) {
	    x = VGetPixel(src,b,r,cc);
	    sum += x * kernel[k];
	  }
	  k++;
	}
	VPixel(dest,b,r,c,VFloat) = sum;
      }
    }
  }
  return dest;
}



/*!
\fn VImage VConvolveRow (VImage src,VImage dest,VImage kernel)
\brief 1D convolution in row-direction (for separable filters)

\param src  input image  (any repn)
\param dest output image  (float repn)
\param kernel convolution kernel  (float repn)
*/
VImage VConvolveRow (VImage src,VImage dest,float *kernel,int dim)
{
  int b,r,c,k,nbands,nrows,ncols;
  int r0,r1,rr;
  float sum,x;

  if (dim%2 == 0) VError("VConvolveRow: kernel dim must be an odd number, dim= %d",dim);
  int d = dim/2;

  nrows  = VImageNRows (src);
  ncols  = VImageNColumns (src);
  nbands = VImageNBands (src);

  dest = VSelectDestImage("VConvolveRow",dest,nbands,nrows,ncols,VFloatRepn);
  VFillImage(dest,VAllBands,0);

  for (b=0; b<nbands; b++) {
    for (r=d; r<nrows-d; r++) {
      for (c=0; c<ncols; c++) {

	sum = 0;
	r0 = r-d;
	r1 = r+d;
	k=0;
	for (rr=r0; rr<=r1; rr++) {	 
	  if (rr >= 0 && rr < nrows) {
	    x = VGetPixel(src,b,rr,c);
	    sum += x * kernel[k];
	    k++;
	  }
	}
	VPixel(dest,b,r,c,VFloat) = sum;
      }
    }
  }
  return dest;
}


/*!
\fn VImage VConvolveBand (VImage src,VImage dest,VImage kernel)
\brief 1D convolution in slice-direction (for separable filters)
\param src  input image (any repn)
\param dest output image  (float repn)
\param kernel convolution kernel  (float repn)
*/
VImage VConvolveBand (VImage src,VImage dest,float *kernel,int dim)
{
  int b,r,c,nbands,nrows,ncols;
  int b0,b1,bb,k;
  float sum,x;

  if (dim%2 == 0) VError("VConvolveBand: kernel dim must be an odd number, dim= %d",dim);
  int d = dim/2;


  nrows  = VImageNRows (src);
  ncols  = VImageNColumns (src);
  nbands = VImageNBands (src);

  dest = VSelectDestImage("VConvolveBand",dest,nbands,nrows,ncols,VFloatRepn);
  VFillImage(dest,VAllBands,0);

  for (b=0; b<nbands; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {

	sum = 0;
	b0 = b-d;
	b1 = b+d;
	k=0;
	for (bb=b0; bb<=b1; bb++) {
	  if (bb >= 0 && bb < nbands) {
	    x = VGetPixel(src,bb,r,c);
	    sum += x * kernel[k];
	  }
	  k++;
	}
	VPixel(dest,b,r,c,VFloat) = sum;
      }
    }
  }
  return dest;
}


