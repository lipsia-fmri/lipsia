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
\param repn  output pixel repn (VShortRepn or VFloatRepn). If 'short' is used, then
the distance values are multiplied by a factor of 10.
*/
VImage 
VChamferDist3d(VImage src,VImage dest,VRepnKind repn)
{
  int nbands,nrows,ncols,b,r,c;
  int i,npixels;
  VShort id1=3,id2=4,id3=5;
  VFloat d1,d2,d3;
  VShort iu,ix;
  VFloat u,x;
  VBit *src_pp;
  VShort *short_pp,sinf;
  VFloat *float_pp,finf;
  double a;

  /* Check the input image */
  if (VPixelRepn(src) != VBitRepn) 
    VError("input image must be of type VBit");
  
  if (repn != VShortRepn && repn != VFloatRepn)
    VError("output pixel repn must be short or float.");

  nbands = VImageNBands(src);
  nrows  = VImageNRows (src);
  ncols  = VImageNColumns (src);

  if (dest == NULL)
    dest = VCreateImage (nbands,nrows,ncols,repn);
  if (! dest) return NULL;


  /* ini dest image */
  npixels = nbands * nrows * ncols;

  if (repn == VShortRepn) {

    sinf = VPixelMaxValue(dest) - 10;
    src_pp  = (VBit *) VPixelPtr (src, 0, 0, 0);
    short_pp = (VShort *) VPixelPtr (dest, 0, 0, 0);
    for (i=0; i<npixels; i++) 
      *short_pp++ = (*src_pp++ > 0) ? 0 : sinf;

    id1 = 3; /* optimal distances (Borgefors, 1984). */
    id2 = 4;
    id3 = 5;
  }
  else {

    finf = VPixelMaxValue(dest) - 10;
    src_pp   = (VBit *) VPixelPtr (src, 0, 0, 0);
    float_pp = (VFloat *) VPixelPtr (dest, 0, 0, 0);
    for (i=0; i<npixels; i++) 
      *float_pp++ = (*src_pp++ > 0) ? 0 : finf;
    
    d1 = 0.9016; /* optimal distances (see Kiryati, 1993). */
    d2 = 1.289;
    d3 = 1.615;
  }


  if (repn == VShortRepn) {

    /* forward scan */
    for (b=0; b<nbands; b++) {
      for (r=0; r<nrows; r++) {
	for (c=0; c<ncols; c++) {
	  iu = VPixel(dest,b,r,c,VShort);
	  if (iu == 0) continue;

	  if (b < 1 || r < 1 || c < 1) ix = sinf;
	  else ix = VPixel(dest,b-1,r-1,c-1,VShort) + id3;
	  if (ix < iu) iu = ix;

	  if (b < 1 || r < 1) ix = sinf;
	  else ix = VPixel(dest,b-1,r-1,c,VShort) + id2;
	  if (ix < iu) iu = ix;

	  if (b < 1 || r < 1 || c >= ncols - 1) ix = sinf;
	  else ix = VPixel(dest,b-1,r-1,c+1,VShort) + id3;
	  if (ix < iu) iu = ix;

	  if (b < 1 || c < 1) ix = sinf;
	  else ix = VPixel(dest,b-1,r,c-1,VShort) + id2;
	  if (ix < iu) iu = ix;

	  if (b < 1) ix = sinf;
	  else ix = VPixel(dest,b-1,r,c,VShort) + id1;
	  if (ix < iu) iu = ix;

	  if (b < 1 || c >= ncols - 1) ix = sinf;
	  else ix = VPixel(dest,b-1,r,c+1,VShort) + id2;
	  if (ix < iu) iu = ix;

	  if (b < 1 || r >= nrows - 1 || c < 1) ix = sinf;
	  else ix = VPixel(dest,b-1,r+1,c-1,VShort) + id3;
	  if (ix < iu) iu = ix;

	  if (b < 1 || r >= nrows - 1) ix = sinf;
	  else ix = VPixel(dest,b-1,r+1,c,VShort) + id2;
	  if (ix < iu) iu = ix;

	  if (b < 1 || r >= nrows - 1 || c >= ncols - 1) ix = sinf;
	  else ix = VPixel(dest,b-1,r+1,c+1,VShort) + id3;
	  if (ix < iu) iu = ix;

	  if (r < 1 || c < 1) ix = sinf;
	  else ix = VPixel(dest,b,r-1,c-1,VShort) + id2;
	  if (ix < iu) iu = ix;

	  if (r < 1) ix = sinf;
	  else ix = VPixel(dest,b,r-1,c,VShort) + id1;
	  if (ix < iu) iu = ix;

	  if (r < 1 || c >= ncols - 1) ix = sinf;
	  else ix = VPixel(dest,b,r-1,c+1,VShort) + id2;
	  if (ix < iu) iu = ix;

	  if (c < 1) ix = sinf;
	  else ix = VPixel(dest,b,r,c-1,VShort) + id1;
	  if (ix < iu) iu = ix;

	  VPixel(dest,b,r,c,VShort) = iu;
	}
      }
    }

    /* backward pass */
    for (b=nbands-1; b>=0; b--) {
      for (r=nrows-1; r>=0; r--) {
	for (c=ncols-1; c>=0; c--) {
	  iu = VPixel(dest,b,r,c,VShort);
	  if (iu == 0) continue;
	
	  if (c >= ncols - 1) ix = sinf;
	  else ix = VPixel(dest,b,r,c+1,VShort) + id1;
	  if (ix < iu) iu = ix;

	  if (r >= nrows - 1 || c < 1) ix = sinf;
	  else ix = VPixel(dest,b,r+1,c-1,VShort) + id2;
	  if (ix < iu) iu = ix;

	  if (r >= nrows - 1) ix = sinf;
	  else ix = VPixel(dest,b,r+1,c,VShort) + id1;
	  if (ix < iu) iu = ix;

	  if (r >= nrows - 1 || c >= ncols - 1) ix = sinf;
	  else ix = VPixel(dest,b,r+1,c+1,VShort) + id2;
	  if (ix < iu) iu = ix;
	
	  if (b >= nbands - 1 || r < 1 || c < 1) ix = sinf;
	  else ix = VPixel(dest,b+1,r-1,c-1,VShort) + id3;
	  if (ix < iu) iu = ix;

	  if (b >= nbands - 1 || r < 1) ix = sinf;
	  else ix = VPixel(dest,b+1,r-1,c,VShort) + id2;
	  if (ix < iu) iu = ix;
	
	  if (b >= nbands - 1 || r < 1 || c >= ncols - 1) ix = sinf;
	  else ix = VPixel(dest,b+1,r-1,c+1,VShort) + id3;
	  if (ix < iu) iu = ix;

	  if (b >= nbands - 1 || c < 1) ix = sinf;
	  else ix = VPixel(dest,b+1,r,c-1,VShort) + id2;
	  if (ix < iu) iu = ix;

	  if (b >= nbands - 1) ix = sinf;
	  else ix = VPixel(dest,b+1,r,c,VShort) + id1;
	  if (ix < iu) iu = ix;

	  if (b >= nbands - 1 || c >= ncols - 1) ix = sinf;
	  else ix = VPixel(dest,b+1,r,c+1,VShort) + id2;
	  if (ix < iu) iu = ix;

	  if (b >= nbands - 1 || r >= nrows - 1 || c < 1) ix = sinf;
	  else ix = VPixel(dest,b+1,r+1,c-1,VShort) + id3;
	  if (ix < iu) iu = ix;

	  if (b >= nbands - 1 || r >= nrows - 1) ix = sinf;
	  else ix = VPixel(dest,b+1,r+1,c,VShort) + id2;
	  if (ix < iu) iu = ix;

	  if (b >= nbands - 1 || r >= nrows - 1 || c >= ncols - 1) ix = sinf;
	  else ix = VPixel(dest,b+1,r+1,c+1,VShort) + id3;
	  if (ix < iu) iu = ix;

	  VPixel(dest,b,r,c,VShort) = iu;
	}
      }
    }

    a = (double) 10.0 / (double) 3.0;
    short_pp = (VShort *) VPixelPtr(dest,0,0,0);
    for (i=0; i<npixels; i++) {
      *short_pp = (VShort) VRint((double) (a * (double) (*short_pp)));
      short_pp++; 
    }
  }
  

  else if (repn == VFloatRepn) {
      
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
  }

  /* Let the destination inherit any attributes of the source image: */
  VCopyImageAttrs (src, dest);
  return dest;
}
