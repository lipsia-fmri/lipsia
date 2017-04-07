/*! \file
  3D binary morphology using distance transforms

3D binary morphological operations are implemented using distance
transforms. This implementation is usually much faster than the 
Minkowski addition. However, only structuring elements of spherical
shape are permitted.

\par Reference:
G. Lohmann (1998). "Volumetric Image Analysis",
John Wiley & Sons, Chichester, England.

\par Author:
 Gabriele Lohmann, MPI-CBS
*/


/* From the Vista library: */
#include <viaio/Vlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <via/via.h>



/*!
  \fn VImage VDTErode(VImage src,VImage dest,VDouble radius)
  \brief 3D morphological erosion
  \param src   input image (bit repn)
  \param dest  output image (bit repn)
  \param radius radius of the spherical structural element
*/
VImage
VDTErode(VImage src,VImage dest,VDouble radius)
{
  VImage float_image;
  VBit *bin_pp,*src_pp;
  VFloat *float_pp;
  int i,nbands,nrows,ncols,npixels;

  if (VPixelRepn(src) != VBitRepn) 
    VError("Input image must be of type VBit");

  nbands  = VImageNBands(src);
  nrows   = VImageNRows(src);
  ncols   = VImageNColumns(src);
  npixels = nbands * nrows * ncols;

  dest = VSelectDestImage("VDTErode",dest,nbands,nrows,ncols,VBitRepn);
  if (! dest) return NULL;

  bin_pp    = (VBit *) VPixelPtr(dest,0,0,0);
  src_pp    = (VBit *) VPixelPtr(src,0,0,0);
  for (i=0; i<npixels; i++) {
    *bin_pp++ = (*src_pp++ > 0 ? 0 : 1);
  }

  float_image = VChamferDist3d(dest,NULL,VFloatRepn);
  if (! float_image) 
    VError(" VDTErode failed.\n");

  float_pp  = (VFloat *) VPixelPtr(float_image,0,0,0);
  bin_pp    = (VBit *) VPixelPtr(dest,0,0,0);
  for (i=0; i<npixels; i++)
    *bin_pp++ = ((*float_pp++ < radius) ? 0 : 1);

  VDestroyImage(float_image);

  VCopyImageAttrs (src, dest);
  return dest;
}




/*!
  \fn VImage VDTDilate(VImage src,VImage dest,VDouble radius)
  \brief 3D morphological dilation
  \param src   input image (bit repn)
  \param dest  output image (bit repn)
  \param radius radius of the spherical structural element
*/
VImage
VDTDilate(VImage src,VImage dest,VDouble radius)
{
  VImage float_image;
  VBit *bin_pp;
  VFloat *float_pp;
  int i,nbands,nrows,ncols,npixels;

  if (VPixelRepn(src) != VBitRepn) 
    VError("Input image must be of type VBit");

  nbands  = VImageNBands(src);
  nrows   = VImageNRows(src);
  ncols   = VImageNColumns(src);
  npixels = nbands * nrows * ncols;

  float_image = VChamferDist3d(src,NULL,VFloatRepn);
  if (! float_image) 
    VError("VDTDilate failed.\n");

  dest = VSelectDestImage("VDTDilate",dest,nbands,nrows,ncols,VBitRepn);
  if (! dest) return NULL;

  float_pp  = (VFloat *) VPixelPtr(float_image,0,0,0);
  bin_pp    = (VBit *) VPixelPtr(dest,0,0,0);
  for (i=0; i<npixels; i++)
    *bin_pp++ = ((*float_pp++ > radius) ? 0 : 1);

  VDestroyImage(float_image);

  VCopyImageAttrs (src,dest);
  return dest;
}


/*!
  \fn VImage VDTClose(VImage src,VImage dest,VDouble radius)
  \brief 3D morphological closing (dilation+erosion)
  \param src   input image (bit repn)
  \param dest  output image (bit repn)
  \param radius radius of the spherical structural element
*/
VImage
VDTClose(VImage src,VImage dest,VDouble radius)
{
  VImage float_image=NULL,tmp=NULL;
  VBit *bin_pp;
  VFloat *float_pp;
  int i,nbands,nrows,ncols,npixels,b,r,c;
  int border = 7;

  border = (int) (radius - 1);

  if (VPixelRepn(src) != VBitRepn) 
    VError("Input image must be of type VBit");

  nbands  = VImageNBands(src) + 2*border;
  nrows   = VImageNRows(src) + 2*border;
  ncols   = VImageNColumns(src) + 2*border;
  npixels = nbands * nrows * ncols;

  tmp = VCreateImage(nbands,nrows,ncols,VBitRepn);

  for (b=border; b<nbands-border; b++) {
    for (r=border; r<nrows-border; r++) {
      for (c=border; c<ncols-border; c++) {
	VPixel(tmp,b,r,c,VBit) = VPixel(src,b-border,r-border,c-border,VBit);
      }
    }
  }

  float_image = VChamferDist3d(tmp,NULL,VFloatRepn);
  if (! float_image) 
    VError("VDTClose failed.\n");

  float_pp  = (VFloat *) VPixelPtr(float_image,0,0,0);
  bin_pp    = (VBit *) VPixelPtr(tmp,0,0,0);
  for (i=0; i<npixels; i++)
    *bin_pp++ = ((*float_pp++ > radius) ? 1 : 0);

  float_image = VChamferDist3d(tmp,float_image,VFloatRepn);
  if (! float_image) 
    VError("VDTClose failed.\n");

  float_pp = (VFloat *) VPixelPtr(float_image,0,0,0);
  bin_pp   = (VBit *) VPixelPtr(tmp,0,0,0);
  for (i=0; i<npixels; i++)
    *bin_pp++ = ((*float_pp++ > radius) ? 1 : 0);

  VDestroyImage(float_image);


  dest = VSelectDestImage("VDTClose",dest,VImageNBands(src),VImageNRows(src),VImageNColumns(src),
			  VBitRepn);
			  
  if (! dest) return NULL;

  for (b=border; b<nbands-border; b++) {
    for (r=border; r<nrows-border; r++) {
      for (c=border; c<ncols-border; c++) {
	VPixel(dest,b-border,r-border,c-border,VBit) = VPixel(tmp,b,r,c,VBit);
      }
    }
  }
  VDestroyImage(tmp);

  VCopyImageAttrs (src, dest);
  return dest;
}



/*!
  \fn VImage VDTOpen(VImage src,VImage dest,VDouble radius)
  \brief 3D morphological closing (erosion+dilation)
  \param src   input image (bit repn)
  \param dest  output image (bit repn)
  \param radius radius of the spherical structural element
*/
VImage
VDTOpen(VImage src,VImage dest,VDouble radius)
{
  VImage float_image=NULL;
  VBit *bin_pp,*src_pp;
  VFloat *float_pp;
  int i,nbands,nrows,ncols,npixels;

  if (VPixelRepn(src) != VBitRepn) 
    VError("Input image must be of type VBit");

  nbands  = VImageNBands(src);
  nrows   = VImageNRows(src);
  ncols   = VImageNColumns(src);
  npixels = nbands * nrows * ncols;

  dest = VSelectDestImage("VDTOpen",dest,nbands,nrows,ncols,VBitRepn);
  if (! dest) return NULL;

  bin_pp    = (VBit *) VPixelPtr(dest,0,0,0);
  src_pp    = (VBit *) VPixelPtr(src,0,0,0);
  for (i=0; i<npixels; i++) {
    *bin_pp++ = (*src_pp++ > 0 ? 0 : 1);
  }

  float_image = VChamferDist3d(dest,NULL,VFloatRepn);
  if (! float_image) VError(" VDTOpen failed.\n");

  float_pp  = (VFloat *) VPixelPtr(float_image,0,0,0);
  bin_pp    = (VBit *) VPixelPtr(dest,0,0,0);
  for (i=0; i<npixels; i++)
    *bin_pp++ = ((*float_pp++ < radius) ? 0 : 1);

  float_image = VChamferDist3d(dest,float_image,VFloatRepn);
  if (! float_image) VError("VDTOpen failed.\n");

  float_pp = (VFloat *) VPixelPtr(float_image,0,0,0);
  bin_pp   = (VBit *) VPixelPtr(dest,0,0,0);
  for (i=0; i<npixels; i++)
    *bin_pp++ = ((*float_pp++ > radius) ? 0 : 1);

  VDestroyImage(float_image);

  VCopyImageAttrs (src, dest);
  return dest;
}

