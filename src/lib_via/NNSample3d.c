/*! \file
  Geometric transformations using nearest neighbour interpolation.


This file contains functions for resampling
and 3D rotations. The transformation equation is:

   y = A(x-x0) + b

where x,x0,b,y are 1x3 vectors and A is a 3x3 matrix.
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



/*!
\fn VImage VNNSample3d (VImage src,VImage dest,VImage transform,
     float b0,float r0,float c0,int dst_nbands,int dst_nrows,int dst_ncolumns)
\brief Resample a 3D image using nearest neighbour interpolation.

\param src   input image (any repn)
\param dest  output image (any repn)
\param transform  4x3 transformation image (float or double repn).
The first column of <transform> contains the translation vector.
The remaining three columns contains the 3x3 linear transformation matrix.
\param b0            slice address that remains fixed 
\param r0            row address that remains fixed 
\param c0            column address that remains fixed 
\param dst_nbands    number of output slices
\param dst_nrows     number of output rows
\param dst_ncolumns  number of output columns
*/
VImage 
VNNSample3d(VImage src,VImage dest,VImage transform,
		   float b0,float r0,float c0,
		   int dst_nbands,int dst_nrows,int dst_ncolumns)
{
  int src_nrows, src_ncolumns, src_nbands; 
  VRepnKind repn;
  int b,r,c;
  float bx,rx,cx;
  float a[3][3],ainv[3][3],detA;
  int i,j;
  int bb,rr,cc;
  float bp,rp,cp;
  float shift[3];
  float val;


  if (VPixelRepn(transform) != VFloatRepn && VPixelRepn(transform) != VDoubleRepn)
    VError("transform image must be float or double repn");
  
  /* Extract data from source image */
  src_nrows    = VImageNRows(src);
  src_ncolumns = VImageNColumns(src);
  src_nbands   = VImageNBands(src);
  repn = VPixelRepn(src);


  /*
  ** create output image
  */
  dest = VSelectDestImage("VNNSample3d",dest,dst_nbands,dst_nrows,dst_ncolumns,repn);
  if (! dest) return NULL;
  VFillImage(dest,VAllBands,0);


  /* get transformation matrix : */
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      a[i][j]  = VGetPixel(transform,0,i,j+1);
    }
  }

  /* get its inverse : */
  ainv[0][0] =  a[1][1]*a[2][2] - a[1][2]*a[2][1];
  ainv[1][0] = -a[1][0]*a[2][2] + a[1][2]*a[2][0];
  ainv[2][0] =  a[1][0]*a[2][1] - a[1][1]*a[2][0];

  ainv[0][1] = -a[0][1]*a[2][2] + a[0][2]*a[2][1];
  ainv[1][1] =  a[0][0]*a[2][2] - a[0][2]*a[2][0];
  ainv[2][1] = -a[0][0]*a[2][1] + a[0][1]*a[2][0];

  ainv[0][2] =  a[0][1]*a[1][2] - a[0][2]*a[1][1];
  ainv[1][2] = -a[0][0]*a[1][2] + a[0][2]*a[1][0];
  ainv[2][2] =  a[0][0]*a[1][1] - a[0][1]*a[1][0];

  /* determinant */
  detA = a[0][0]*ainv[0][0] + a[0][1]*ainv[1][0] + a[0][2]*ainv[2][0];
  if (detA == 0) VError(" VNNSample3d: transformation matrix is singular");

  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      ainv[i][j] /= detA;
    }
  }


  /* get shift vector */
  shift[0] = VGetPixel(transform,0,0,0);
  shift[1] = VGetPixel(transform,0,1,0);
  shift[2] = VGetPixel(transform,0,2,0);


  /* 
  ** resampling
  */
  for (b=0; b<dst_nbands; b++) {
    for (r=0; r<dst_nrows; r++) {
      for (c=0; c<dst_ncolumns; c++) {

	bx = (float) b - shift[0];
	rx = (float) r - shift[1];
	cx = (float) c - shift[2];

        bp = ainv[0][0] * bx + ainv[0][1] * rx + ainv[0][2] * cx;
        rp = ainv[1][0] * bx + ainv[1][1] * rx + ainv[1][2] * cx;
        cp = ainv[2][0] * bx + ainv[2][1] * rx + ainv[2][2] * cx;

	bp += b0;
	rp += r0;
	cp += c0;
	
	if (bp <= 1 || bp >= src_nbands-1) continue;
	if (rp <= 1 || rp >= src_nrows-1) continue;
	if (cp <= 1 || cp >= src_ncolumns-1) continue;

	bb = (int) (bp + 0.5);
        rr = (int) (rp + 0.5);
        cc = (int) (cp + 0.5);
        if (bb < 0 || bb >= src_nbands) continue;
        if (rr < 0 || rr >= src_nrows) continue;
        if (cc < 0 || cc >= src_ncolumns) continue;

	val = VGetPixel(src,bb,rr,cc);
        if (repn == VUByteRepn) val = VRint(val);
        VSetPixel(dest,b,r,c,(VDouble)val);
      }
    }
  }

  VCopyImageAttrs (src, dest);
  return dest;
}

