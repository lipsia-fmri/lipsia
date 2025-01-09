/*! \file
3d chamfer distance transform.

For each background voxel, the length of the shortest
3D path to the nearest foreground voxel is computed.
The chamfer distance metric is an approximation to the
Euclidian distance.

\par References:
G. Borgefors (1984).
"Distance Transforms in arbitrary dimensions",
CVGIP 27, pp.321-345.<br>
A.L.D. Beckers, A.W.M. Smeulders (1992),
"Optimization of Length measurements for isotropic distance transformations
in three dimensions",
CVGIP: Image understanding, Vol. 55, No.3, pp- 296-306.

\par Author:
Gabriele Lohmann, MPI-CBS
*/


/* From the Vista library: */
#include <viaio/Vlib.h>
#include <via/via.h>
#include <stdio.h>
#include <math.h>

#define IMIN(a,b) ((a) < (b) ? (a) : (b))

/*!
\fn VImage VChamferDist3d(VImage src,VImage dest,VRepnKind repn)
\param src   input image (bit repn)
\param dest  output image (short of float repn)
*/
VImage  VCDist3d(VImage src,VImage dest)
{
  int nbands,nrows,ncols,b,r,c;
  size_t i;
  VFloat d1,d2,d3;
  VFloat u,x;

  /* Check the input image */
  if (VPixelRepn(src) != VBitRepn) 
    VError("input image must be of type VBit");

  nbands = VImageNBands(src);
  nrows  = VImageNRows (src);
  ncols  = VImageNColumns (src);

  if (dest == NULL)
    dest = VCreateImage (nbands,nrows,ncols,VFloatRepn);
  if (! dest) return NULL;
  VFillImage(dest,VAllBands,0);


  /* ini dest image */
  size_t npixels = nbands * nrows * ncols;
  VFloat finf = VPixelMaxValue(dest) - 1;
  VBit *src_pp   = VImageData(src);
  VFloat *float_pp = VImageData(dest);
  for (i=0; i<npixels; i++) 
    float_pp[i] = (src_pp[i] > 0) ? 0 : finf;
    
  d1 = 0.9016; /* optimal distances (see Kiryati, 1993). */
  d2 = 1.289;
  d3 = 1.615;

      
  /* forward scan */
  for (b=0; b<nbands; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	u = VPixel(dest,b,r,c,VFloat);
	if (u == 0) continue;

	if (b < 1 || r < 1 || c < 1) x = finf;
	else x = VPixel(dest,b-1,r-1,c-1,VFloat) + d3;
	if (x < u) u = x;

	if (b < 1 || r < 1) x = finf;
	else x = VPixel(dest,b-1,r-1,c,VFloat) + d2;
	if (x < u) u = x;

	if (b < 1 || r < 1 || c >= ncols - 1) x = finf;
	else x = VPixel(dest,b-1,r-1,c+1,VFloat) + d3;
	if (x < u) u = x;

	if (b < 1 || c < 1) x = finf;
	else x = VPixel(dest,b-1,r,c-1,VFloat) + d2;
	if (x < u) u = x;

	if (b < 1) x = finf;
	else x = VPixel(dest,b-1,r,c,VFloat) + d1;
	if (x < u) u = x;

	if (b < 1 || c >= ncols - 1) x = finf;
	else x = VPixel(dest,b-1,r,c+1,VFloat) + d2;
	if (x < u) u = x;

	if (b < 1 || r >= nrows - 1 || c < 1) x = finf;
	else x = VPixel(dest,b-1,r+1,c-1,VFloat) + d3;
	if (x < u) u = x;

	if (b < 1 || r >= nrows - 1) x = finf;
	else x = VPixel(dest,b-1,r+1,c,VFloat) + d2;
	if (x < u) u = x;

	if (b < 1 || r >= nrows - 1 || c >= ncols - 1) x = finf;
	else x = VPixel(dest,b-1,r+1,c+1,VFloat) + d3;
	if (x < u) u = x;

	if (r < 1 || c < 1) x = finf;
	else x = VPixel(dest,b,r-1,c-1,VFloat) + d2;
	if (x < u) u = x;

	if (r < 1) x = finf;
	else x = VPixel(dest,b,r-1,c,VFloat) + d1;
	if (x < u) u = x;

	if (r < 1 || c >= ncols - 1) x = finf;
	else x = VPixel(dest,b,r-1,c+1,VFloat) + d2;
	if (x < u) u = x;

	if (c < 1) x = finf;
	else x = VPixel(dest,b,r,c-1,VFloat) + d1;
	if (x < u) u = x;

	VPixel(dest,b,r,c,VFloat) = u;
      }
    }
  }

  /* backward pass */
  for (b=nbands-1; b>=0; b--) {
    for (r=nrows-1; r>=0; r--) {
      for (c=ncols-1; c>=0; c--) {
	u = VPixel(dest,b,r,c,VFloat);
	if (u == 0) continue;
	
	if (c >= ncols - 1) x = finf;
	else x = VPixel(dest,b,r,c+1,VFloat) + d1;
	if (x < u) u = x;

	if (r >= nrows - 1 || c < 1) x = finf;
	else x = VPixel(dest,b,r+1,c-1,VFloat) + d2;
	if (x < u) u = x;

	if (r >= nrows - 1) x = finf;
	else x = VPixel(dest,b,r+1,c,VFloat) + d1;
	if (x < u) u = x;

	if (r >= nrows - 1 || c >= ncols - 1) x = finf;
	else x = VPixel(dest,b,r+1,c+1,VFloat) + d2;
	if (x < u) u = x;
	
	if (b >= nbands - 1 || r < 1 || c < 1) x = finf;
	else x = VPixel(dest,b+1,r-1,c-1,VFloat) + d3;
	if (x < u) u = x;

	if (b >= nbands - 1 || r < 1) x = finf;
	else x = VPixel(dest,b+1,r-1,c,VFloat) + d2;
	if (x < u) u = x;
	
	if (b >= nbands - 1 || r < 1 || c >= ncols - 1) x = finf;
	else x = VPixel(dest,b+1,r-1,c+1,VFloat) + d3;
	if (x < u) u = x;

	if (b >= nbands - 1 || c < 1) x = finf;
	else x = VPixel(dest,b+1,r,c-1,VFloat) + d2;
	if (x < u) u = x;

	if (b >= nbands - 1) x = finf;
	else x = VPixel(dest,b+1,r,c,VFloat) + d1;
	if (x < u) u = x;

	if (b >= nbands - 1 || c >= ncols - 1) x = finf;
	else x = VPixel(dest,b+1,r,c+1,VFloat) + d2;
	if (x < u) u = x;

	if (b >= nbands - 1 || r >= nrows - 1 || c < 1) x = finf;
	else x = VPixel(dest,b+1,r+1,c-1,VFloat) + d3;
	if (x < u) u = x;

	if (b >= nbands - 1 || r >= nrows - 1) x = finf;
	else x = VPixel(dest,b+1,r+1,c,VFloat) + d2;
	if (x < u) u = x;

	if (b >= nbands - 1 || r >= nrows - 1 || c >= ncols - 1) x = finf;
	else x = VPixel(dest,b+1,r+1,c+1,VFloat) + d3;
	if (x < u) u = x;

	VPixel(dest,b,r,c,VFloat) = u;
      }
    }
  }

  VCopyImageAttrs (src, dest);
  return dest;
}
