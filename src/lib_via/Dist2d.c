/*! \file
2d chamfer distance transform.

For each foreground pixel, the length of the shortest
2D path to the nearest background pixel is computed.
The chamfer distance metric is an approximation to the
Euclidian distance.


\par References:
G. Borgefors (1984).
"Distance Transforms in arbitrary dimensions",
CVGIP 27, pp.321-345.

\par Author:
Gabriele Lohmann, MPI-CBS
*/


/* From the Vista library: */
#include <viaio/Vlib.h>
#include <viaio/file.h>
#include <viaio/mu.h>
#include <viaio/option.h>
#include <viaio/VImage.h>

/* From the standard C library: */
#include <stdio.h>
#include <math.h>



/*!
\fn VImage VChamferDist2d (VImage src,VImage dest,VBand band)
\param src   input image (bit repn)
\param dest  output image (float repn)
\param band  id of band to be processed. To select all bands, use 'VAllBands'.
*/
VImage VChamferDist2d(VImage src,VImage dest,VBand band)
{
  int b,r,c,nbands,ncols,nrows,npixels,i;
  VBit *src_pp;
  VFloat *dest_pp,inf;
  VFloat x,z,d1,d2;
  
  /* Check the destination image
     If it is NULL, then create one of the appropriate size and type. */

  if (VPixelRepn(src) != VBitRepn)
    VError("VChamferDist2d: input image must be of type bit.");

  nrows    = VImageNRows (src);
  ncols    = VImageNColumns (src);
  npixels  = nrows * ncols;

  nbands = (band == VAllBands) ? VImageNBands (src) : 1;
  dest = VSelectDestImage ("VDist2dImage", dest, nbands,nrows,ncols,
			   VFloatRepn);
  if (! dest) return NULL;
  

  
  d1 = 1; 
  d2 = 1;        /* checkerboard metric         */
  
  d1 = 1;
  d2 = 1.351;    /* siehe Borgefors, 1984       */

  d1 = 0.95509;  /* nach Borgefors, 1986, p.351 */ 
  d2 = 1.36930;

  inf = VPixelMaxValue(dest);    /* infinity  */
  
  for (b = 0; b < nbands; b++) {
    
    src_pp  = (VBit *)   VPixelPtr (src, b, 0, 0);
    dest_pp = (VFloat *) VPixelPtr (dest, b, 0, 0);
    
    for (i=0; i<npixels; i++) 
      *dest_pp++ = (*src_pp++ > 0) ? 0 : inf;
    
    /* Forward pass: */
    for (r=1; r<nrows-1; r++) {
      for (c=1; c<ncols-1; c++) {

	x = VPixel(dest,b,r,c,VFloat);
	if (x < 0.00001) continue;
	
	z = VPixel(dest,b,r-1,c-1,VFloat) + d2;
	if (z < x) x = z;
	z = VPixel(dest,b,r-1,c,VFloat) + d1;
	if (z < x) x = z;
	z = VPixel(dest,b,r-1,c+1,VFloat) + d2;
	if (z < x) x = z;
	z = VPixel(dest,b,r,c-1,VFloat) + d1;
	if (z < x) x = z;

	VPixel(dest,b,r,c,VFloat) = x;
      }
    }

    /* Backward pass: */
    for (r=nrows-2; r>1; r--) {
      for (c=ncols-2; c>1; c--) {

	x = VPixel(dest,b,r,c,VFloat);
	if (x < 0.00001) continue;

	z = VPixel(dest,b,r,c+1,VFloat) + d1;
	if (z < x) x = z;
	z = VPixel(dest,b,r+1,c-1,VFloat) + d2;
	if (z < x) x = z;
	z = VPixel(dest,b,r+1,c,VFloat) + d1;
	if (z < x) x = z;
	z = VPixel(dest,b,r+1,c+1,VFloat) + d2;
	if (z < x) x = z;

	VPixel(dest,b,r,c,VFloat) = x;
      }
    }
  }

  VCopyImageAttrs (src, dest);
  return dest;
}
