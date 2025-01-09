/*! \file
3d Euclidean distance transform.

For each background voxel, the length of the shortest
3D path to the nearest foreground voxel is computed.


\par References:
Toyofumi Saito, Jun-Ichiro Toriwaki (1994).
"New algorithms for euclidean distance transformation of a n-dimensional
picture with applications",
Pattern Recognition, Vol.27, No.11, pp. 1551-1565.

\par Author:
Gabriele Lohmann, MPI-CBS
*/


/* From the Vista library: */
#include <viaio/Vlib.h>
#include <via/via.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define IMIN(a,b) ((a) < (b) ? (a) : (b))
#define IMAX(a,b) ((a) > (b) ? (a) : (b))

extern VImage VEDistShort3d(VImage,VImage);
extern VImage VEDistFloat3d(VImage,VImage);




/*
** Euclidean distance transform 
*/
VImage VEDist3d(VImage src,VImage dest)
{
  int b,r,c,bb,rr;
  int i,istart,iend;
  int nbands,nrows,ncols,npixels;
  int d,d1,d2,cc1,cc2;
  VFloat u,dmin,dmax,*destpix;
  VBit *srcpix;
  double g,*array;
  
  nbands = VImageNBands(src);
  nrows  = VImageNRows(src);
  ncols  = VImageNColumns(src);

  npixels = IMAX(nbands,nrows);
  array = (double *) VMalloc(sizeof(double) * npixels);

  npixels = VImageNPixels(src);

  dest = VSelectDestImage("VEDistFloat3d",dest,nbands,nrows,ncols,VFloatRepn);
  if (! dest) return NULL;
  VFillImage(dest,VAllBands,0);

  dmax = VPixelMaxValue(dest);

  /* first pass */
  for (b=0; b<nbands; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {

	if (VPixel(src,b,r,c,VBit) == 1) {
	  VPixel(dest,b,r,c,VFloat) = 0;
	  continue;
	}

	srcpix = (VBit *) VPixelPtr(src,b,r,c);
	cc1 = c;
	while (cc1 < ncols && *srcpix++ == 0) 
	  cc1++;
	d1 = (cc1 >= ncols ? ncols : (cc1 - c));

	srcpix = (VBit *) VPixelPtr(src,b,r,c);
	cc2 = c;
	while (cc2 >= 0  && *srcpix-- == 0) 
	  cc2--;
	d2 = (cc2 <= 0 ? ncols : (c - cc2));

	if (d1 <= d2) {
	  d  = d1;
	}
	else {
	  d  = d2;
	}
	VPixel(dest,b,r,c,VFloat) = (VFloat) (d * d);
      }
    }
  }

  /* second pass */
  for (b=0; b<nbands; b++) {
    for (c=0; c<ncols; c++) {

      for (r=0; r<nrows; r++)
	array[r] = (double) VPixel(dest,b,r,c,VFloat);

      for (r=0; r<nrows; r++) {
	if (VPixel(src,b,r,c,VBit) == 1) continue;
	
	dmin = dmax;
	g = sqrt(array[r]);
	istart = r - (int) g;
	if (istart < 0) istart = 0;
	iend = r + (int) g + 1;
	if (iend >= nrows)
	  iend = nrows;

	for (rr=istart; rr<iend; rr++) {
	  u = array[rr] + (r - rr) * (r - rr);
	  if (u < dmin) {
	    dmin = u;
	  }
	}
	VPixel(dest,b,r,c,VFloat) = dmin;
      }
    }
  }

  /* third pass */

  for (r=0; r<nrows; r++) {
    for (c=0; c<ncols; c++) {

      for (b=0; b<nbands; b++)
	array[b] = (double) VPixel(dest,b,r,c,VFloat);

      for (b=0; b<nbands; b++) {
	if (VPixel(src,b,r,c,VBit) == 1) continue;
	
	dmin = dmax;
	g = sqrt(array[b]);
	istart = b - (int) g - 1;
	if (istart < 0) istart = 0;
	iend = b + (int) g + 1;
	if (iend >= nbands)
	  iend = nbands;

	for (bb=istart; bb<iend; bb++) {
	  u = array[bb] + (b - bb) * (b - bb);
	  if (u < dmin) {
	    dmin = u;
	  }
	}
	VPixel(dest,b,r,c,VFloat) = dmin;
      }
    }
  }
  VFree(array);

  destpix = VImageData(dest);
  for (i=0; i<npixels; i++) {
    destpix[i] = sqrt((double) (destpix[i]));
  }

  VCopyImageAttrs (src, dest);
  return dest;
}
