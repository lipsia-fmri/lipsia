/*! \file
  3D border extraction. 


Extract 3D border voxels from a binary raster image.
A border voxel is a foreground voxel that is 6-connected
to a background voxel.
 
\par Author:
Gabriele Lohmann, MPI-CBS
*/

/* From the Vista library: */
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <stdio.h>
#include <math.h>


/*!
\fn VImage VBorderImage3d (VImage src,VImage dest)
\param src   input image (bit repn)
\param dest  output image (bit repn)
*/
VImage 
VBorderImage3d (VImage src,VImage dest)
{
  int nbands,nrows,ncols,b,r,c;

  if (VPixelRepn(src) != VBitRepn)
    VError("VBorderImage3d: input pixel repn must be bit");

  nbands  = VImageNBands (src);
  nrows   = VImageNRows (src);
  ncols   = VImageNColumns (src);
  

  dest = VSelectDestImage("VBorderImage3d",dest,nbands,nrows,ncols,VBitRepn);
  if (! dest) VError(" err creating dest image");
  VFillImage(dest,VAllBands,0);


  for (b=1; b<nbands-1; b++) {
    for (r=1; r<nrows-1; r++) {
      for (c=1; c<ncols-1; c++) {
		       
	if (VPixel(src,b,r,c,VBit) == 0) continue;

	if (VPixel(src,b-1,r,c,VBit) == 0) {
	  VPixel(dest,b,r,c,VBit) = 1;
	  continue;
	}
	if (VPixel(src,b,r-1,c,VBit) == 0) {
  	  VPixel(dest,b,r,c,VBit) = 1;
	  continue;
	}
	if (VPixel(src,b,r,c-1,VBit) == 0) {
	  VPixel(dest,b,r,c,VBit) = 1;
	  continue;
	}
	if (VPixel(src,b+1,r,c,VBit) == 0) {
	  VPixel(dest,b,r,c,VBit) = 1;
	  continue;
	}
	if (VPixel(src,b,r+1,c,VBit) == 0) {
	  VPixel(dest,b,r,c,VBit) = 1;
	  continue;
	}
	if (VPixel(src,b,r,c+1,VBit) == 0) {
	  VPixel(dest,b,r,c,VBit) = 1;
	  continue;
	}
      }
    }
  }

  VCopyImageAttrs (src, dest);
  return dest;
}
