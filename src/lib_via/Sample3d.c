/*! \file
  Geometric transformations using trilinear interpolation.


This file contains functions for trilinear resampling
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
\fn VImage VTriLinearSample3d (VImage src,VImage dest,VImage transform,
     float b0,float r0,float c0,int dst_nbands,int dst_nrows,int dst_ncolumns)
\brief Resample a 3D image using trilinear interpolation.

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
VTriLinearSample3d(VImage src,VImage dest,VImage transform,
		   float b0,float r0,float c0,
		   int dst_nbands,int dst_nrows,int dst_ncolumns)
{
  int src_nrows, src_ncols, src_nbands; 
  VRepnKind repn;
  int   b,r,c;
  int   i,j;
  float bp,rp,cp,bx,rx,cx;
  float a[3][3],ainv[3][3],detA;
  float shift[3];
  float val;

  int   sx, sy, sz;   /* origin of subcube    */
  float px, py, pz;   /* fractions of subcube */
  float qx, qy, qz;   /* fractions of subcube */
  int lx, ly, lz;   /* lengths */
  int ox, oy, oz;   /* offsets */


  if (VPixelRepn(transform) != VFloatRepn && VPixelRepn(transform) != VDoubleRepn)
    VError("transform image must be float or double repn");
  
  /* Extract data from source image */
  src_nrows   = VImageNRows(src);
  src_ncols   = VImageNColumns(src);
  src_nbands  = VImageNBands(src);
  repn = VPixelRepn(src);


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
  if (detA == 0) VError(" VTriLinearSample3d: transformation matrix is singular");

  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      ainv[i][j] /= detA;
    }
  }


  /* get translation vector */
  shift[0] = VGetPixel(transform,0,0,0);
  shift[1] = VGetPixel(transform,0,1,0);
  shift[2] = VGetPixel(transform,0,2,0);



#define GetValues(type) \
{ \
  type *src_pp; \
  src_pp = (type *) VPixelPtr (src, sz, sy, sx); \
  val += (float) pz * py * px * *src_pp; src_pp += ox; \
  val += (float) pz * py * qx * *src_pp; src_pp += oy; \
  val += (float) pz * qy * px * *src_pp; src_pp += ox; \
  val += (float) pz * qy * qx * *src_pp; src_pp += oz; \
  val += (float) qz * py * px * *src_pp; src_pp += ox; \
  val += (float) qz * py * qx * *src_pp; src_pp += oy; \
  val += (float) qz * qy * px * *src_pp; src_pp += ox; \
  val += (float) qz * qy * qx * *src_pp; \
  VPixel(dest,b,r,c,type) = val; \
}


  /*
  ** create output image
  */
  dest = VSelectDestImage("VTriLinearSample3d",dest,dst_nbands,dst_nrows,dst_ncolumns,repn);
  if (! dest) return NULL;
  VFillImage(dest,VAllBands,0);

  /* Determines the value of each pixel in the destination image: */
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
	
	if (bp < 0 || bp > src_nbands) continue;
	if (rp < 0 || rp > src_nrows) continue;
	if (cp < 0 || cp > src_ncols) continue;

	/* compute origin of subcube */
	sx = (int) (cp);
	sy = (int) (rp);
	sz = (int) (bp);

	/* check subcube */
	if ((sx < -1) || (sx > src_ncols  - 1)) continue;
	if ((sy < -1) || (sy > src_nrows  - 1)) continue;
	if ((sz < -1) || (sz > src_nbands - 1)) continue;

	/* compute fractions of subcube */
	qx = cp - sx; px = 1 - qx;
	qy = rp - sy; py = 1 - qy;
	qz = bp - sz; pz = 1 - qz;

	/* compute lengths and offsets */
	lx = 1;
	ly = src_ncols;
	lz = src_nrows * src_ncols;
	if (sx == -1) {sx = 0; lx = 0;};
	if (sy == -1) {sy = 0; ly = 0;};
	if (sz == -1) {sz = 0; lz = 0;};
	if (sx == src_ncols  - 1) lx = 0;
	if (sy == src_nrows  - 1) ly = 0;
	if (sz == src_nbands - 1) lz = 0;
	ox = lx;
	oy = ox + ly - 2 * lx;
	oz = oy + lz - 2 * ly;
	val = 0;

	switch(repn) {
	case VShortRepn:
	  GetValues(VShort);
	  break;
	case VUByteRepn:
	  GetValues(VUByte);
	  break;
	case VFloatRepn:
	  GetValues(VFloat);
	  break;
	case VSByteRepn:
	  GetValues(VSByte);
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

