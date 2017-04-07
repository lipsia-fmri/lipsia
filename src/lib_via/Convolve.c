/*! \file
  Convolution

This file contains code for 2D and 3D convolutions, and for
1D convolutions. 1D convolutions can be used to implement
separable 2D/3D filters.

\par Author:
Gabriele Lohmann, MPI-CBS
*/


/* From the Vista library: */
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <via/via.h>

/* From the standard C libaray: */
#include <stdio.h>
#include <stdlib.h>




/*!
\fn VImage VConvolve3d (VImage src,VImage dest,VImage kernel)
\brief 3D convolution 
\param src    input image  (any repn)
\param dest   output image (float repn)
\param kernel raster image containing convolution kernel (float repn)
*/
/*
** 3d convolution
*/
VImage
VConvolve3d (VImage src,VImage dest,VImage kernel)
{
  int b,r,c,nbands,nrows,ncols;
  int b0,b1,r0,r1,c0,c1,bb,rr,cc;
  VFloat sum,*float_pp;
  int dimb,dimr,dimc,db,dr,dc;


  if (VPixelRepn(kernel) != VFloatRepn) VError(" kernel pixel repn must be float");

  dimc = VImageNColumns(kernel);
  dimr = VImageNRows(kernel);
  dimb = VImageNBands(kernel);


  if (dimc%2 == 0) VError("VConvolve3d: kernel dim must be an odd number (%d)",dimc);
  if (dimr%2 == 0) VError("VConvolve3d: kernel dim must be an odd number (%d)",dimr);
  if (dimb%2 == 0) VError("VConvolve3d: kernel dim must be an odd number (%d)",dimb);

  dc = dimc/2;
  dr = dimr/2;
  db = dimb/2;


  nrows  = VImageNRows (src);
  ncols  = VImageNColumns (src);
  nbands = VImageNBands (src);

  dest   = VSelectDestImage("VConvolve3d",dest,nbands,nrows,ncols,VFloatRepn);
  VFillImage(dest,VAllBands,0);
  VCopyImageAttrs (src, dest);

  for (b=db; b<nbands-db; b++) {
    for (r=dr; r<nrows-dr; r++) {
      for (c=dc; c<ncols-dc; c++) {

	float_pp = (VFloat *) VImageData(kernel);
	sum = 0;
	b0 = b-db;
	b1 = b+db;
	for (bb=b0; bb<=b1; bb++) {

	  r0 = r-dr;
	  r1 = r+dr;
	  for (rr=r0; rr<=r1; rr++) {

	    c0 = c-dc;
	    c1 = c+dc;
	    for (cc=c0; cc<=c1; cc++) {
	      sum += VReadPixel(src,bb,rr,cc) * (*float_pp++);
	    }
	  }
	}
	VPixel(dest,b,r,c,VFloat) = sum;
      }
    }
  }
  return dest;
}



/*!
\fn VImage VConvolve2d (VImage src,VImage dest,VImage kernel)
\brief 2D convolution 
\param src    input image (any repn)
\param dest   output image (float repn)
\param kernel raster image containing convolution kernel (float repn)
*/
/*
** 2d convolution
*/
VImage
VConvolve2d (VImage src,VImage dest,VImage kernel)
{
  int b,r,c,nbands,nrows,ncols;
  int r0,r1,c0,c1,rr,cc;
  float sum;
  int dimr,dimc,dr,dc;
  VFloat *float_pp;

  if (VPixelRepn(kernel) != VFloatRepn) VError(" kernel pixel repn must be float");
  if (VImageNBands(kernel) > 1) VError(" kernel must be 2D");

  dimc = VImageNColumns(kernel);
  dimr = VImageNRows(kernel);
  if (dimc%2 == 0) VError("VConvolve2d: kernel dim must be an odd number");
  if (dimr%2 == 0) VError("VConvolve2d: kernel dim must be an odd number");

  dc = dimc/2;
  dr = dimr/2;

  nrows  = VImageNRows (src);
  ncols  = VImageNColumns (src);
  nbands = VImageNBands (src);
  
  dest = VSelectDestImage("VConvolve2d",dest,nbands,nrows,ncols,VFloatRepn);
  dest = VCopyImage(src,dest,VAllBands);

  for (b=0; b<nbands; b++) {
    for (r=dr; r<nrows-dr; r++) {
      for (c=dc; c<ncols-dc; c++) {

	float_pp = (VFloat *) VImageData(kernel);
	sum = 0;
	r0 = r-dr;
	r1 = r+dr;
	for (rr=r0; rr<=r1; rr++) {

	  c0 = c-dc;
	  c1 = c+dc;
	  for (cc=c0; cc<=c1; cc++) {
	    sum += VReadPixel(src,b,rr,cc) * (*float_pp++);
	  }
	}
	VPixel(dest,b,r,c,VFloat) = sum;
      }
    }
  }
  return dest;
}




/*!
\fn VImage VConvolveCol (VImage src,VImage dest,VImage kernel)
\brief 1D convolution in column-direction (for separable filters)
\param src    input image  (any repn)
\param dest   output image (float repn)
\param kernel raster image containing convolution kernel (float repn)
*/
VImage
VConvolveCol (VImage src,VImage dest,VImage kernel)
{
  int b,r,c,nbands,nrows,ncols;
  int c0,c1,cc;
  float sum,x;
  VFloat *float_pp;
  int dim,d;


  if (VPixelRepn(src) != VFloatRepn) VError(" input pixel repn must be float");

  nrows  = VImageNRows (src);
  ncols  = VImageNColumns (src);
  nbands = VImageNBands (src);

  dest = VSelectDestImage("VConvolveCol",dest,nbands,nrows,ncols,VFloatRepn);
  dest = VCopyImage(src,dest,VAllBands);


  dim = VImageNColumns(kernel);
  if (dim%2 == 0) VError("VConvolveCol: kernel dim must be an odd number");
  d = dim/2;

  for (b=0; b<nbands; b++) {
    for (r=0; r<nrows; r++) {
      for (c=d; c<ncols-d; c++) {

	float_pp = (VFloat *) VImageData(kernel);
	sum = 0;
	c0 = c-d;
	c1 = c+d;
	for (cc=c0; cc<=c1; cc++) {
	  x = VReadPixel(src,b,r,cc);
	  sum += x * (*float_pp++);
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
VImage
VConvolveRow (VImage src,VImage dest,VImage kernel)
{
  int b,r,c,nbands,nrows,ncols;
  int r0,r1,rr;
  float sum,x;
  VFloat *float_pp;
  int d,dim;

  if (VPixelRepn(kernel) != VFloatRepn) VError(" kernel pixel repn must be float");

  dim = VImageNColumns(kernel);
  if (dim%2 == 0) VError("VConvolveRow: kernel dim must be an odd number");
  d = dim/2;


  nrows  = VImageNRows (src);
  ncols  = VImageNColumns (src);
  nbands = VImageNBands (src);


  dest = VSelectDestImage("VConvolveRow",dest,nbands,nrows,ncols,VFloatRepn);
  dest = VCopyImage(src,dest,VAllBands);

  for (b=0; b<nbands; b++) {
    for (r=d; r<nrows-d; r++) {
      for (c=0; c<ncols; c++) {

	float_pp = (VFloat *) VImageData(kernel);
	sum = 0;
	r0 = r-d;
	r1 = r+d;
	for (rr=r0; rr<=r1; rr++) {
	  x = VReadPixel(src,b,rr,c);
	  sum += x * (*float_pp++);
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
VImage
VConvolveBand (VImage src,VImage dest,VImage kernel)
{
  int b,r,c,nbands,nrows,ncols;
  int b0,b1,bb;
  float sum,x;
  VFloat *float_pp;
  int d,dim;

  if (VPixelRepn(kernel) != VFloatRepn) VError(" kernel pixel repn must be float");

  dim = VImageNColumns(kernel);
  if (dim%2 == 0) VError("VConvolveBand: kernel dim must be an odd number");
  d = dim/2;


  nrows  = VImageNRows (src);
  ncols  = VImageNColumns (src);
  nbands = VImageNBands (src);

  dest = VSelectDestImage("VConvolveBand",dest,nbands,nrows,ncols,VFloatRepn);
  dest = VCopyImage(src,dest,VAllBands);

  for (b=d; b<nbands-d; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {

	float_pp = (VFloat *) VImageData(kernel);
	sum = 0;
	b0 = b-d;
	b1 = b+d;
	for (bb=b0; bb<=b1; bb++) {
	  x = VReadPixel(src,bb,r,c);
	  sum += x * (*float_pp++);
	}
	VPixel(dest,b,r,c,VFloat) = sum;
      }
    }
  }
  return dest;
}


