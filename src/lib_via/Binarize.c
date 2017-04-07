/*! \file
  binarize an input image
    
\par Author:
Gabriele Lohmann, MPI-CBS
*/

/* From the Vista library: */
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <stdio.h>
#include <via/via.h>

/*!
\fn VImage VBinarizeImage (VImage src,VImage dest,VDouble xmin,VDouble xmax)
\param src  input image (any repn except bit)
\param dest  output image (bit repn)
\param xmin  lower threshold, voxels less than (or equal to) xmin are set to zero
\param xmax  upper threshold, voxels larger than (or equal to) xmax are set to zero
*/

#define binarize(type) \
{ \
  type *src_pp; \
  dest_pp = (VBit *) VImageData(dest); \
  src_pp  = (type *) VImageData(src); \
  for (i=0; i<npixels; i++) { \
    u = *src_pp++; \
    *dest_pp = (u >= xmin) && (u <= xmax); \
    dest_pp++; \
  } \
}


VImage VBinarizeImage (VImage src,VImage dest,VDouble xmin,VDouble xmax)
{
  VBit *dest_pp;
  int i,npixels;
  VFloat u;

  /*
  if (xmin < VPixelMinValue (src) || xmax > VPixelMaxValue (src))
    VWarning ("BinarizeImage: Thresholds (%g %g) are outside pixel range [%g,%g]",
    xmin,xmax,VPixelMaxValue (src), VPixelMaxValue (src));
  */
  dest = VSelectDestImage("VBinarizeImage",dest,
			  VImageNBands(src),VImageNRows(src),VImageNColumns(src),
			  VBitRepn);

  npixels = VImageNPixels(src);

  switch(VPixelRepn(src)) {

  case VUByteRepn:
    binarize(VUByte);
    break;
  case VSByteRepn:
    binarize(VSByte);
    break;
  case VShortRepn:
    binarize(VShort);
    break;
  case VUShortRepn:
    binarize(VShort);
    break;
  case VIntegerRepn:
    binarize(VShort);
    break;
  case VUIntegerRepn:
    binarize(VShort);
    break;
  case VLongRepn:
    binarize(VLong);
    break;
  case VULongRepn:
    binarize(VLong);
    break;
  case VFloatRepn:
    binarize(VFloat);
    break;
  case VDoubleRepn:
    binarize(VDouble);
    break;

  default:
    VError(" Binarize: illegal pixel repn");
  }

  VCopyImageAttrs (src, dest);
  return dest;
}
