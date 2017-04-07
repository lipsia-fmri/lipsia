/*
 *  $Id: Invert.c 3177 2008-04-01 14:47:24Z karstenm $
 *
 *  This file contains routines for inverting and negating images.
 */

/*
 *  Copyright 1993, 1994 University of British Columbia
 *
 *  Permission to use, copy, modify, distribute, and sell this software and its
 *  documentation for any purpose is hereby granted without fee, provided that
 *  the above copyright notice appears in all copies and that both that
 *  copyright notice and this permission notice appear in supporting
 *  documentation. UBC makes no representations about the suitability of this
 *  software for any purpose. It is provided "as is" without express or
 *  implied warranty.
 *
 *  Author: Arthur Pope, UBC Laboratory for Computational Intelligence
 */

/* From the Vista library: */
#include "viaio/Vlib.h"
#include "viaio/os.h"
#include "viaio/VImage.h"

/* File identification string: */
VRcsId ("$Id: Invert.c 3177 2008-04-01 14:47:24Z karstenm $");


/*  
 *  VInvertImage
 *
 *  Performs grey-scale inversion of an image:
 *
 *	dest-pixel = M - src-pixel
 *
 *  where M is 1 the maximum pixel value.
 *  Each src-pixel is tested to ensure that it is in the range [0, M].
 */

VImage VInvertImage (VImage src, VImage dest, VBand band)
{
  int npixels, i;
  VPointer src_pixels;
  VBoolean clipped = FALSE;

#define Invert(type, m)						\
	{								\
	    type *src_pp = src_pixels, *dest_pp = VImageData (dest);	\
	    for (i = 0; i < npixels; i++)				\
		*dest_pp++ = m - *src_pp++;				\
	}

#define InvertClipZero(type, m)					\
	{								\
	    type pixel, *src_pp = src_pixels, *dest_pp = VImageData (dest); \
	    for (i = 0; i < npixels; i++) {				\
		pixel = *src_pp++;					\
		if (pixel < 0) {					\
		    *dest_pp++ = m;					\
		    clipped = TRUE;					\
		} else *dest_pp++ = m - pixel;				\
	    }								\
	}

#define InvertClipBoth(type, m)					\
	{								\
	    type pixel, *src_pp = src_pixels, *dest_pp = VImageData (dest); \
	    for (i = 0; i < npixels; i++) {				\
		pixel = *src_pp++;					\
		if (pixel < 0) {					\
		    *dest_pp++ = m;					\
		    clipped = TRUE;					\
		} else if (pixel > m) {					\
		    *dest_pp++ = 0;					\
		    clipped = TRUE;					\
		} else *dest_pp++ = m - pixel;				\
	    }								\
	}

  /* Locate the source and destination pixels: */
  if (! VSelectBand ("VInvertImage", src, band, & npixels, & src_pixels))
    return NULL;
  dest = VSelectDestImage ("VInvertImage", dest,
			   band == VAllBands ? VImageNBands (src) : 1,
			   VImageNRows (src), VImageNColumns (src),
			   VPixelRepn (src));
  if (! dest)
    return NULL;

  /* Invert each pixel: */
  switch (VPixelRepn (src)) {

  case VBitRepn:
    Invert (VBit, 1);
    break;

  case VUByteRepn:
    Invert (VUByte, 255);
    break;

  case VSByteRepn:
    InvertClipZero (VSByte, 127);
    break;

  case VShortRepn:
    InvertClipZero (VShort, 32767);
    break;

  case VIntegerRepn:
    InvertClipZero (VInteger, 2147483647);
    break;

  case VUIntegerRepn:
    InvertClipZero (VInteger, 4294967295);
    break;

  case VLongRepn:
    InvertClipZero (VLong, 9223372036854775807);
    break;

  case VULongRepn:
    InvertClipZero (VULong, 1.844674e+19);
    break;

  case VFloatRepn:
    InvertClipBoth (VFloat, 1.0);
    break;

  case VDoubleRepn:
    InvertClipBoth (VDouble, 1.0);
    break;

  default:
    break;
  }

  if (clipped)
    VWarning ("VInvertImage: Source pixel(s) out of range");

  VCopyImageAttrs (src, dest);

  return dest;

#undef Invert
#undef InvertClipZero
#undef InvertClipBoth
}


/*  
 *  VNegateImage
 *
 *  Negates each pixel of an image. Valid only for images with signed
 *  pixel representation.
 */

VImage VNegateImage (VImage src, VImage dest, VBand band)
{
  int npixels, i;
  VPointer src_pixels;
  VBoolean clipped = FALSE;

#define Negate(type)						\
	{								\
	    type *src_pp = src_pixels, *dest_pp = VImageData (dest);	\
	    for (i = 0; i < npixels; i++)				\
		*dest_pp++ = - *src_pp++;				\
	}

#define NegateClip(type, m)						\
	{								\
	    type pixel, *src_pp = src_pixels, *dest_pp = VImageData (dest); \
	    for (i = 0; i < npixels; i++) {				\
		pixel = *src_pp++;					\
		if (pixel < -m) {					\
		    *dest_pp++ = m;					\
		    clipped = TRUE;					\
		} else *dest_pp++ = - pixel;				\
	    }								\
	}

  /* The image to be negated must have signed pixels: */
  if (VPixelRepn (src) == VBitRepn || VPixelRepn (src) == VUByteRepn) {
    VWarning ("VNegateImage: Source image must have signed pixels");
    return NULL;
  }

  /* Locate the source and destination pixels: */
  if (! VSelectBand ("VNegateImage", src, band, & npixels, & src_pixels))
    return NULL;
  dest = VSelectDestImage ("VNegateImage", dest,
			   band == VAllBands ? VImageNBands (src) : 1,
			   VImageNRows (src), VImageNColumns (src),
			   VPixelRepn (src));
  if (! dest)
    return NULL;

  /* Invert each pixel: */
  switch (VPixelRepn (src)) {

  case VSByteRepn:
    NegateClip (VSByte, 127);
    break;

  case VShortRepn:
    NegateClip (VShort, 32767);
    break;

  case VLongRepn:
    NegateClip (VLong, 2147483647);
    break;

  case VFloatRepn:
    Negate (VFloat);
    break;

  case VDoubleRepn:
    Negate (VDouble);
    break;

  default:
    break;
  }

  if (clipped)
    VWarning ("VNegateImage: Source pixel(s) out of range");

  VCopyImageAttrs (src, dest);

  return dest;

#undef Negate
#undef NegateClip
}
