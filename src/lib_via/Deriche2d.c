/*! \file
 2D Deriche filter for edge detection
 
\par Reference:
R. Deriche. Fast algorithms for low-level vision. 
IEEE Transactions on Pattern Analysis and Machine Intelligence, 
1(12):78-88, January 1990.

\par Author:
Gabriele Lohmann, MPI-CBS
*/


/* From the Vista library: */
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

/* From the standard C libaray: */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/*!
\fn void VDeriche2d(VImage src,VFloat alpha,VImage *gradr,VImage *gradc);
\param src     input image (ubyte repn)
\param alpha   parameter controlling edge strength
\param *gradr  output gradient in row direction (float repn)
\param *gradc  output gradient in column direction (float repn)
*/
void
VDeriche2d (VImage src,VFloat alpha,VImage *gradr,VImage *gradc)
{
  int b,r,c;
  int nbands,nrows,ncols,npixels,len;
  double s,a,a0,a1,a2,a3,b1,b2,exp_alpha;
  float *left,*right;
  
  nbands  = VImageNBands (src);
  nrows   = VImageNRows (src);
  ncols   = VImageNColumns (src);
  npixels = nbands * nrows * ncols;

  len = (ncols > nrows) ? ncols : nrows;

  left  = (float *) VMalloc(sizeof(float) * len);
  right = (float *) VMalloc(sizeof(float) * len);

 
  exp_alpha = exp((double) ( - alpha));
  a  = 1.0 * exp_alpha;
  b1 = -2.0 * exp_alpha;
  b2 = exp((double) (-2.0 * alpha));

  s = ((1.0 - exp_alpha) * (1.0 - exp_alpha)) / (1.0 + 2.0 * alpha * exp_alpha - b2);

  a0 = s;
  a1 = s * (alpha - 1.0) * exp_alpha;
  a2 = a1 - s * b1;
  a3 = - s * b2;


  /* col-grad */

  *gradc = VCreateImage(nbands,nrows,ncols,VFloatRepn);

  for (b=0; b<nbands; b++) {
    for (r=0; r<nrows; r++) {

      /* left-to-right */
      left[0] = 0;
      left[1] = 0;
      for (c=2; c<ncols; c++) {
	left[c] = (VFloat) VGetPixel(src,b,r,c-1) 
	  - b1 * left[c-1]  - b2 * left[c-2];
      }

      /* right-to-left */
      right[ncols-1] = 0;
      right[ncols-2] = 0;
      for (c=ncols-3; c>=0; c--) {
	right[c] = (VFloat) VGetPixel(src,b,r,c+1) 
	  - b1 * right[c+1]  - b2 * right[c+2];
      }

      /* combine */

      for (c=0; c<ncols; c++) {
	VPixel(*gradc,b,r,c,VFloat) =  a * (left[c] - right[c]);
      }
    }
  }

  /* row-smooth */

  for (b=0; b<nbands; b++) {
    for (c=0; c<ncols; c++) {

      /* left-to-right */
      left[0] = 0;
      left[1] = 0;
      for (r=2; r<nrows; r++) {
	left[r] = 
	    a0 * (VFloat) VPixel(*gradc,b,r,c,VFloat)
	  + a1 * (VFloat) VPixel(*gradc,b,r-1,c,VFloat)
	  - b1 * left[r-1] - b2 * left[r-2];
      }

      /* right-to-left */
      right[nrows-1] = 0;
      right[nrows-2] = 0;

      for (r=nrows-3; r>=0; r--) {
	right[r] = 
	    a2 * (VFloat) VPixel(*gradc,b,r+1,c,VFloat)
	  + a3 * (VFloat) VPixel(*gradc,b,r+2,c,VFloat)
	  - b1 * right[r+1] - b2 * right[r+2];
      }
      /* combine */

      for (r=0; r<nrows; r++) {
	VPixel(*gradc,b,r,c,VFloat) += left[r] + right[r];
      }
    }
  }



  /* 
  ** row-grad 
  */
  *gradr = VCreateImage(nbands,nrows,ncols,VFloatRepn);

  /* row-deriv */
  for (b=0; b<nbands; b++) {
    for (c=0; c<ncols; c++) {

      /* left-to-right */
      left[0] = 0;
      left[1] = 0;
      for (r=2; r<nrows; r++) {
	left[r] = (VFloat) VGetPixel(src,b,r-1,c) 
	  - b1 * left[r-1]  - b2 * left[r-2];
      }

      /* right-to-left */
      right[ncols-1] = 0;
      right[ncols-2] = 0;
      for (r=nrows-3; r>=0; r--) {
	right[r] = (VFloat) VGetPixel(src,b,r+1,c) 
	  - b1 * right[r+1]  - b2 * right[r+2];
      }

      /* combine */

      for (r=0; r<nrows; r++) {
	VPixel(*gradr,b,r,c,VFloat) =  a * (left[r] - right[r]);
      }
    }
  }

  /* col-smooth */

  for (b=0; b<nbands; b++) {
    for (r=0; r<nrows; r++) {

      /* left-to-right */
      left[0] = 0;
      left[1] = 0;
      for (c=2; c<ncols; c++) {
	left[c] = 
	    a0 * (VFloat) VPixel(*gradr,b,r,c,VFloat)
	  + a1 * (VFloat) VPixel(*gradr,b,r,c-1,VFloat)
	  - b1 * left[c-1] - b2 * left[c-2];
      }

      /* right-to-left */
      right[ncols-1] = 0;
      right[ncols-2] = 0;

      for (c=ncols-3; c>=0; c--) {
	right[r] = 
	    a2 * (VFloat) VPixel(*gradr,b,r,c+1,VFloat)
	  + a3 * (VFloat) VPixel(*gradr,b,r,c+2,VFloat)
	  - b1 * right[c+1] - b2 * right[c+2];
      }
      /* combine */

      for (c=0; c<ncols; c++) {
	VPixel(*gradr,b,r,c,VFloat) += left[c] + right[c];
      }
    }
  }


  VFree(left);
  VFree(right);
}

 
