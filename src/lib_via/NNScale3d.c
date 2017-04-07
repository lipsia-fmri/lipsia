/*! \file
  3D scaling using nearest neighbour interpolation.

Scale a 3D image using nearest neighbour interpolation.

\par Author:
Gabriele Lohmann, MPI-CBS
*/

/* From the Vista library: */
#include "viaio/Vlib.h"
#include "viaio/mu.h"
#include "viaio/os.h"
#include "viaio/VImage.h"

/* From the standard C library: */
#include <stdio.h>
#include <math.h>


/*!
\fn VImage VNNScale3d (VImage src, VImage dest,int dst_nbands,int dst_nrows,int dst_ncols,
                       float shift[3],float scale[3])
\brief 3D scaling using nearest neighbour interpolation.
\param src        input image (any repn)
\param dest       output image (any repn)
\param dst_nbands number of output slices
\param dst_nrows  number of output rows
\param dst_ncols  number of output columns
\param shift[3]   translation vector (band,row,column)
\param scale[3]   scaling vector (band,row,column)
*/
VImage 
VNNScale3d (VImage src,VImage dest,int dst_nbands,int dst_nrows,int dst_ncols,
	    float shift[3],float scale[3])
{
  int nrows,ncols,nbands; 
  VRepnKind repn;
  int b,r,c,bb,rr,cc;
  VDouble v;
  float xscale,yscale,zscale;
  float bp,rp,cp;


  /* Extract data from source image */
  nrows   = VImageNRows(src);
  ncols   = VImageNColumns(src);
  nbands  = VImageNBands(src);
  repn    = VPixelRepn(src);
    
  dest = VSelectDestImage("VNNScale3d",dest,dst_nbands,dst_nrows,dst_ncols,repn);
  if (! dest) return NULL;
  VFillImage(dest,VAllBands,0);

  zscale = 1.0 / scale[0];
  yscale = 1.0 / scale[1];
  xscale = 1.0 / scale[2];


  for (b=0; b<dst_nbands; b++) {
    bp = zscale * ((float) b - shift[0]);
    bb = (int) (bp + 0.5);
    if (bb < 0 || bb >= nbands) continue;

    for (r=0; r<dst_nrows; r++) {
      rp = yscale * ((float) r  - shift[1]);
      rr = (int) (rp + 0.5);
      if (rr < 0 || rr >= nrows) continue;

      for (c=0; c<dst_ncols; c++) {
	cp = xscale * ((float) c - shift[2]);
	cc = (int) (cp + 0.5);
	if (cc < 0 || cc >= ncols) continue;
	  
	v = VGetPixel(src,bb,rr,cc);
	if (repn == VUByteRepn) v = (int)(v + 0.4999);
	VSetPixel(dest,b,r,c,v);
      }
    }
  }

  VCopyImageAttrs (src, dest);
  return dest;
}
