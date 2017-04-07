/*! \file
Constrained distance transform.

Perform Chamfer distance transform in 3D while avoiding obstacles, i.e.
get distance from each "inside" point to closest "outside" point
while avoiding "obstacle" points. The distance is measured as the length
of the shortest path to the closest outside point that does not touch an
obstacle point.

\par Reference:
Verwer, B.J., Verbeek, P.W.,  Dekker, S.T. (1989).
An efficient uniform cost algorithm applied to distance transforms.
IEEE Trans. on Pattern Aanalysis and Machine Intelligence,
Vol.11, No.4, pp. 425--429.

\par Author:
Gabriele Lohmann, MPI-CBS
*/

/* From the Vista library: */
#include <viaio/Vlib.h>
#include <via/via.h>
#include <stdio.h>
#include <math.h>

#define IMIN(a,b) ((a) < (b) ? (a) : (b))
#define IMAX(a,b) ((a) > (b) ? (a) : (b))


/*!
\fn VImage VCDT3d (VImage src,VImage dest, VLong inside,VLong outside,VLong obstacle,VRepnKind repn)
\param src      input image (VUByte)
\param dest     output image (VFloat)
\param inside   ubyte value representing inside voxels
\param outside  ubyte value representing outside voxels
\param repn     output pixel repn (VShortRepn or VFloatRepn)
*/
VImage 
VCDT3d (VImage src,VImage dest,
	VLong inside,VLong outside,VRepnKind repn)
{
  int nbands,nrows,ncols,b,r,c;
  int i,n_new,n_old,iter,npixels;
  VFloat d1,d2,d3;
  VFloat u,x,y;
  VUByte *src_pp;
  VFloat *dest_pp,inf;
  VShort *short_pp;
  VImage result=NULL;


  if (repn != VShortRepn && repn != VFloatRepn)
    VError(" VCDT3d: illegal output pixel repn: %d,  short= %d, float= %d",
	   repn,VShortRepn,VFloatRepn);

  /* Check the input image */
  if (VPixelRepn(src) != VUByteRepn) 
    VError("VCDT3d: input image must be of type VUByte");

  
  nbands = VImageNBands(src);
  nrows  = VImageNRows (src);
  ncols  = VImageNColumns (src);
  npixels = nbands * nrows * ncols;

  dest = VSelectDestImage("VCDT3d",dest,nbands,nrows,ncols,VFloatRepn);
  VFillImage(dest,VAllBands,0);
  inf = VPixelMaxValue(dest) - 2;

  dest_pp = (VFloat *) VPixelPtr(dest,0,0,0);
  src_pp = (VUByte *) VPixelPtr(src,0,0,0);
  npixels = nbands * nrows * ncols;
  for (i=0; i<npixels; i++) {
    x = *src_pp++;
    if (x == (VFloat) inside)
      y = inf;
    else if (x == (VFloat) outside)
      y = 0;
    else
      y = -100.0;
    *dest_pp++ = y;
  }

 
  /* optimal chamfer distances (Borgefors,1984, p. 334) */
  d1 = 1.0;
  d2 = 1.314;
  d3 = 1.628;

  /* optimal distances (see Beckers, Smeulders, 1992). */

  d1 = 0.88750;
  d2 = 1.34224;
  d3 = 1.59772;

  /* optimal distances (see Verwer, 1991). */
/*
  d1 = 0.894;
  d2 = 1.3409;
  d3 = 1.5879;
  */

  /* optimal distances (see Kiryati, 1993). */
  /*
  d1 = 0.9016;
  d2 = 1.289;
  d3 = 1.615;
  */

  iter = 0;
  n_old = 1;
  n_new = 0;
  while (n_old != n_new && iter < 6) {

    /* forward scan */

    for (b=1; b<nbands; b++) {
      for (r=1; r<nrows-1; r++) {
	for (c=1; c<ncols-1; c++) {

	  if (VPixel(src,b,r,c,VUByte) != inside) continue;

	  u = VPixel(dest,b,r,c,VFloat);

	  x = VPixel(dest,b-1,r-1,c-1,VFloat) + d3;
	  if (x < u && x >= 0) u = x;

	  x = VPixel(dest,b-1,r-1,c,VFloat) + d2;
	  if (x < u && x >= 0) u = x;

	  x = VPixel(dest,b-1,r-1,c+1,VFloat) + d3;
	  if (x < u && x >= 0) u = x;

	  x = VPixel(dest,b-1,r,c-1,VFloat) + d2;
	  if (x < u && x >= 0) u = x;

	  x = VPixel(dest,b-1,r,c,VFloat) + d1;
	  if (x < u && x >= 0) u = x;

	  x = VPixel(dest,b-1,r,c+1,VFloat) + d2;
	  if (x < u && x >= 0) u = x;

	  x = VPixel(dest,b-1,r+1,c-1,VFloat) + d3;
	  if (x < u && x >= 0) u = x;

	  x = VPixel(dest,b-1,r+1,c,VFloat) + d2;
	  if (x < u && x >= 0) u = x;

	  x = VPixel(dest,b-1,r+1,c+1,VFloat) + d3;
	  if (x < u && x >= 0) u = x;


	  x = VPixel(dest,b,r-1,c-1,VFloat) + d2;
	  if (x < u && x >= 0) u = x;

	  x = VPixel(dest,b,r-1,c,VFloat) + d1;
	  if (x < u && x >= 0) u = x;

	  x = VPixel(dest,b,r-1,c+1,VFloat) + d2;
	  if (x < u && x >= 0) u = x;

	  x = VPixel(dest,b,r,c-1,VFloat) + d1;
	  if (x < u && x >= 0) u = x;

	  if (u >= 0)
	    VPixel(dest,b,r,c,VFloat) = u;
	}
      }
    }

    /* backward scan */

    for (b=nbands-2; b>=0; b--) {
      for (r=nrows-2; r>=1; r--) {
	for (c=ncols-2; c>=1; c--) {

	  if (VPixel(src,b,r,c,VUByte) != inside) continue;

	  u = VPixel(dest,b,r,c,VFloat);
	
	  x = VPixel(dest,b,r,c+1,VFloat) + d1;
	  if (x < u && x >= 0) u = x;

	  x = VPixel(dest,b,r+1,c-1,VFloat) + d2;
	  if (x < u && x >= 0) u = x;

	  x = VPixel(dest,b,r+1,c,VFloat) + d1;
	  if (x < u && x >= 0) u = x;

	  x = VPixel(dest,b,r+1,c+1,VFloat) + d2;
	  if (x < u && x >= 0) u = x;
	

	  x = VPixel(dest,b+1,r-1,c-1,VFloat) + d3;
	  if (x < u && x >= 0) u = x;

	  x = VPixel(dest,b+1,r-1,c,VFloat) + d2;
	  if (x < u && x >= 0) u = x;
	
	  x = VPixel(dest,b+1,r-1,c+1,VFloat) + d3;
	  if (x < u && x >= 0) u = x;


	  x = VPixel(dest,b+1,r,c-1,VFloat) + d2;
	  if (x < u && x >= 0) u = x;

	  x = VPixel(dest,b+1,r,c,VFloat) + d1;
	  if (x < u && x >= 0) u = x;

	  x = VPixel(dest,b+1,r,c+1,VFloat) + d2;
	  if (x < u && x >= 0) u = x;


	  x = VPixel(dest,b+1,r+1,c-1,VFloat) + d3;
	  if (x < u && x >= 0) u = x;

	  x = VPixel(dest,b+1,r+1,c,VFloat) + d2;
	  if (x < u && x >= 0) u = x;

	  x = VPixel(dest,b+1,r+1,c+1,VFloat) + d3;
	  if (x < u && x >= 0) u = x;

	  if (u >= 0)
	    VPixel(dest,b,r,c,VFloat) = u;
	}
      }
    }

    dest_pp = (VFloat *) VPixelPtr(dest,0,0,0);
    n_old = n_new;
    n_new = 0;
    for (i=0; i<npixels; i++) {
      if (*dest_pp++ >= 9999) n_new++;
    }
    iter++;
  }

  dest_pp = (VFloat *) VPixelPtr(dest,0,0,0);
  for (i=0; i<npixels; i++) {
    if (*dest_pp >= 9999 || *dest_pp < 0) *dest_pp = 0;
    dest_pp++;
  }

  /* copy to a short image */
  if (repn == VShortRepn) {
    result = VCreateImage(nbands,nrows,ncols,VShortRepn);
    dest_pp = (VFloat *) VPixelPtr(dest,0,0,0);
    short_pp = (VShort *) VPixelPtr(result,0,0,0);
    for (i=0; i<npixels; i++) {
      *short_pp = VRint((float)(*dest_pp * 10.0f));
      dest_pp++;
      short_pp++;
    }
    VDestroyImage(dest);
    VCopyImageAttrs (src, result);
    return result;
  }
  else if (repn == VFloatRepn) {
    VCopyImageAttrs (src, dest);
    return dest;
  }
  return NULL;
}
