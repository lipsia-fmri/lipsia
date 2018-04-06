
/*! \file
 A non-linear 3d smoothing filter

The smoothing filter computes the most frequent value in a 6 (or 18)
neighbourhood and replaces the center pixel with this value.
The pixels within the neighbourhood are weighted.


\par Author:
Gabriele Lohmann, MPI-CBS
*/


/* From the Vista library: */
#include <viaio/Vlib.h>
#include <via/via.h>
#include <viaio/mu.h>


/* From the standard C library: */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 *  Smooth3d
 * 
 *  apply a voting filter in a small neighbourhood
 *
 */
#define Smooth3d(type) \
{ \
  type *dest_pp,*src_pp; \
  dest_pp = (type *) VPixelPtr(dest,0,0,0); \
  for (b = 0; b < nbands; b++) { \
    b0 = (b < 1) ? 0 : b-1; \
    b1 = (b > nbands - 2) ?  nbands - 1 : b+1; \
    for (r = 0; r < nrows; r++) { \
      r0 = (r < 1) ? 0 : r-1; \
      r1 = (r > nrows - 2) ?  nrows - 1 : r+1; \
      for (c = 0; c < ncols; c++) { \
	c0 = (c < 1) ? 0 : c-1; \
	c1 = (c > ncols - 2) ?  ncols - 1 : c+1; \
        isum = 0; \
        i1 = VPixel(src,b,r,c,type); \
	isum += (int) i1 * n1; \
	isum += (int) VPixel(src,b0,r,c,type) * n6; \
	isum += (int) VPixel(src,b,r0,c,type) * n6; \
	isum += (int) VPixel(src,b,r,c0,type) * n6; \
	isum += (int) VPixel(src,b1,r,c,type) * n6; \
	isum += (int) VPixel(src,b,r1,c,type) * n6; \
	isum += (int) VPixel(src,b,r,c1,type) * n6; \
        if (neighb == 1) { \
	  isum += VPixel(src,b0,r0,c,type) * n18; \
	  isum += (int) VPixel(src,b,r0,c0,type) * n18; \
	  isum += (int) VPixel(src,b0,r,c0,type) * n18; \
	  isum += (int) VPixel(src,b1,r1,c,type) * n18; \
	  isum += (int) VPixel(src,b,r1,c1,type) * n18; \
	  isum += (int) VPixel(src,b1,r,c1,type) * n18; \
	  isum += (int) VPixel(src,b1,r0,c,type) * n18; \
	  isum += (int) VPixel(src,b,r1,c0,type) * n18; \
	  isum += (int) VPixel(src,b1,r,c0,type) * n18; \
	  isum += (int) VPixel(src,b0,r1,c,type) * n18; \
	  isum += (int) VPixel(src,b,r0,c1,type) * n18; \
	  isum += (int) VPixel(src,b0,r,c1,type) * n18; \
        } \
        if (neighb == 2) { \
	  isum += (int) VPixel(src,b0,r0,c0,type) * n26; \
	  isum += (int) VPixel(src,b1,r0,c0,type) * n26; \
	  isum += (int) VPixel(src,b0,r1,c0,type) * n26; \
	  isum += (int) VPixel(src,b0,r0,c1,type) * n26; \
	  isum += (int) VPixel(src,b1,r1,c0,type) * n26; \
	  isum += (int) VPixel(src,b1,r0,c1,type) * n26; \
	  isum += (int) VPixel(src,b0,r1,c1,type) * n26; \
	  isum += (int) VPixel(src,b1,r1,c1,type) * n26; \
        } \
        sum = (double) isum / (double) norm; \
        i0 = (int) VRint((double)sum); \
	*dest_pp++ = (type) i0; \
        if (i0 != i1) n++; \
      } \
    } \
  } \
  if (numiter > 1) { \
    src_pp = (type *) VPixelPtr(src,0,0,0); \
    dest_pp = (type *) VPixelPtr(dest,0,0,0); \
    for (i=0; i<npixels; i++) *src_pp++ = *dest_pp++; \
  } \
}



/*!
\fn VImage VSmoothImage3d (VImage src, VImage dest, VLong neighb, VLong numiter)
\param src  input image (any repn)
\param dest output image (any repn)
\param neighb adjacency type (6,18, or 26)
\param numiter number of iterations (filtering may be applied repeatedly)
*/
VImage 
VSmoothImage3d (VImage src, VImage dest, VLong neighb, VLong numiter)
{
  long nbands,nrows,ncols,npixels;
  VRepnKind repn;
  long i,i0,i1,n,iter;
  long b0,b1,r0,r1,c0,c1,b,r,c;
  long n1,n6,n18,n26;
  int isum;
  double sum=0,norm=0;

  repn   = VPixelRepn (src);
  nbands = VImageNBands (src);
  nrows  = VImageNRows (src);
  ncols  = VImageNColumns (src);
  npixels = nbands * nrows * ncols;
  if (dest == NULL) 
    dest = VCreateImage (nbands,nrows,ncols,repn);
  if (! dest) return NULL;

  n1  = 8;
  n6  = 4;
  n18 = 2;
  n26 = 1;

  switch(neighb) {
  case 0:
    norm = n1 + 6 * n6;
    break;
  case 1:
    norm = n1 + 6 * n6 + 12 * n18;
    break;
  case 2:
    norm = n1 + 6 * n6 + 12 * n18 + 8 * n26;
    break;
  default: ;
  }

  iter = 0;
  n = 100;
  while (n > 1 && iter < numiter) {
    iter++;
    n = 0;
    switch (repn) {

    case VBitRepn:
      Smooth3d(VBit);
      break;
		
    case VUByteRepn:
      Smooth3d(VUByte);
      break;
		
    case VSByteRepn:
      Smooth3d(VSByte);
      break;
		
    case VShortRepn:
      Smooth3d(VShort);
      break;

    case VUShortRepn:
      Smooth3d(VUShort);
      break;

    case VIntegerRepn:
      Smooth3d(VInteger);
      break;

    case VUIntegerRepn:
      Smooth3d(VUInteger);
      break;
		
    case VLongRepn:
      Smooth3d(VLong);
      break;

    case VULongRepn:
      Smooth3d(VULong);
      break;

    case VFloatRepn:
      Smooth3d(VFloat);
      break;

    case VDoubleRepn:
      Smooth3d(VDouble);
      break;

    default:
      VError("Illegal representation type");
    }
  }


  VCopyImageAttrs (src, dest);
  VSetAttr (VImageAttrList(dest), "component_interp", NULL, VStringRepn, "image");
  return dest;
}
