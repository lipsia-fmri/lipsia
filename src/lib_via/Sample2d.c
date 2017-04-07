/*! \file
  Geometric transformations using bilinear interpolation.


This file contains functions for bilinear resampling
and 2D rotations. The transformation equation is:

   y = A(x-x0) + b

where x,x0,b,y are 1x2 vectors and A is a 2x2 matrix.
The vector x0 can be used to specify a position that
remains unchanged by the transformation.


\par Author:
Gabriele Lohmann, MPI-CBS
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

/* From the Vista library: */
#include <viaio/Vlib.h>
#include <viaio/file.h>
#include <viaio/mu.h>
#include <via/via.h>


#define ABS(x) ((x) > 0 ? (x) : -(x))


/*!
\fn VImage VBiLinearSample2d (VImage src,VImage dest,VImage transform,
     float r0,float c0,int dst_nrows,int dst_ncolumns)
\brief Resample a 2D image using bilinear interpolation.

\param src   input image (any repn)
\param dest  output image (any repn)
\param transform  3x2 transformation image (float or double repn).
The first column of <transform> contains the translation vector.
The remaining three columns contains the 2x2 linear transformation matrix.
\param r0            row address that remains fixed 
\param c0            column address that remains fixed 
\param dst_nrows     number of output rows
\param dst_ncolumns  number of output columns
*/
VImage 
VBiLinearSample2d(VImage src,VImage dest,VImage transform,
		   float r0,float c0,
		   int dst_nrows,int dst_ncolumns)
{
  int src_nrows, src_ncols, src_nbands; 
  VRepnKind repn;
  int b,r,c;
  float rp,cp,rx,cx;
  float a[2][2],x[2][2];
  int i,j,n;
  float shift[2];
  float val;
  int   rq,cq;
  float ra,ca,rb,cb;


  if (VPixelRepn(transform) != VFloatRepn && VPixelRepn(transform) != VDoubleRepn)
    VError("transform image must be float or double repn");
  
  /* Extract data from source image */
  src_nrows   = VImageNRows(src);
  src_ncols   = VImageNColumns(src);
  src_nbands  = VImageNBands(src);
  repn = VPixelRepn(src);


  /* get inverse of transformation matrix : */
  n = 2;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      x[i][j] = VGetPixel(transform,0,i,j+1);
    }
  }

  /* invert transformation matrix, Fischer, p. 152 */
  val = (x[0][0]*x[1][1] - x[0][1]*x[1][0]);
  if (val == 0) VError(" VBiLinearSample2d: matrix not invertible");
  a[0][0] =  x[1][1]/val;
  a[0][1] = -x[0][1]/val;
  a[1][0] = -x[1][0]/val;
  a[1][1] =  x[0][0]/val;
  
  shift[0] = VGetPixel(transform,0,0,0);
  shift[1] = VGetPixel(transform,0,1,0);


#define GetValues(type) \
{ \
  val  = rb * cb * VPixel(src,b,rq,cq,type); \
  val += rb * ca * VPixel(src,b,rq,cq+1,type); \
  val += ra * cb * VPixel(src,b,rq+1,cq,type); \
  val += ra * ca * VPixel(src,b,rq+1,cq+1,type); \
}


  /*
  ** create output image
  */
  dest = VSelectDestImage("VBiLinearSample2d",dest,src_nbands,dst_nrows,dst_ncolumns,repn);
  if (! dest) return NULL;
  VFillImage(dest,VAllBands,0);


  /* Determines the value of each pixel in the destination image: */
  for (b=0; b<src_nbands; b++) {
    for (r=0; r<dst_nrows; r++) {
      for (c=0; c<dst_ncolumns; c++) {

	rx = (float) r - shift[0];
	cx = (float) c - shift[1];

	rp = a[0][0] * rx + a[0][1] * cx;
	cp = a[1][0] * rx + a[1][1] * cx;
	
	rp += r0;
	cp += c0;

	rq = (int) rp;
	cq = (int) cp;

	/* check subcube */
	if ((cq < 0) || (cq >= src_ncols  - 1)) continue;
	if ((rq < 0) || (rq >= src_nrows  - 1)) continue;	

	/* compute fractions of subcube */
	ca = cp  - (float)cq;
	cb = 1.0 - ca;
	ra = rp  - (float)rq;
	rb = 1.0 - ra;


	val = 0;
	switch(repn) {
	case VShortRepn:
	  GetValues(VShort);
	  VPixel(dest,b,r,c,VShort) = VRint(val);
	  break;
	case VUByteRepn:
	  GetValues(VUByte);
	  val = VRint(val);
	  if (val <   0) val = 0;
	  if (val > 255) val = 255;
	  VPixel(dest,b,r,c,VUByte) = val;
	  break;
	case VFloatRepn:
	  GetValues(VFloat);
	  VPixel(dest,b,r,c,VUByte) = val;
	  break;
	case VSByteRepn:
	  GetValues(VSByte);
	  VPixel(dest,b,r,c,VUByte) = VRint(val);
	  break;
	default:
	  VError(" illegal pixel repn");
	}

      }
    }
  }

  VCopyImageAttrs (src, dest);
  return dest;
}

