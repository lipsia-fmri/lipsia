/*
 *  $Id: Transpose.c 3177 2008-04-01 14:47:24Z karstenm $
 *
 *  This file contains routines for transposing images
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
 *  Author: Daniel Ko, UBC Laboratory for Computational Intelligence
 */

/* From the Vista library: */
#include "viaio/Vlib.h"
#include "viaio/mu.h"
#include "viaio/os.h"
#include "viaio/VImage.h"

/* From the standard C libaray: */
#include <math.h>



/*
 * VTransposeImage
 *
 * Transpose an image such that row i, column j in the
 * source image becomes row j, column i in the destination
 * image.
 */

VImage VTransposeImage (VImage src, VImage dest, VBand band)
{
  int src_nrows, src_ncols, src_nbands;
  int dest_nrows, dest_ncols, dest_nbands;
  VRepnKind src_repn, dest_repn;
  int r, c, b, src_curband, dest_curband;
  VImage result;

  /* Read properties of "src": */
  src_nrows  = VImageNRows (src);
  src_ncols  = VImageNColumns (src);
  src_nbands = VImageNBands (src);
  src_repn   = VPixelRepn (src);

  /* Check to ensure that "band" exists: */
  if (band != VAllBands && (band < 0 || band >= src_nbands)) {
    VWarning ("VTransposeImage: Band %d does not exist", band);
    return NULL;
  }

  /* Determine properties of "dest": */
  dest_nrows  = src_ncols;
  dest_ncols  = src_nrows;
  dest_nbands = (band == VAllBands) ? src_nbands : 1;
  dest_repn   = src_repn;

  /* Create/select result image */
  result = VSelectDestImage ("VTransposeImage", dest,
			     dest_nbands, dest_nrows, dest_ncols, dest_repn);
    
  if (!result)   /* Cannot create destination image */
    return NULL;

  if (src == dest) {
    result = VCreateImage (dest_nbands, dest_nrows, dest_ncols, dest_repn);
    if (!result)
      return NULL;
  }
    
  /*
   * Transpose an image of a particular type: 
   */
#define Transpose(type)                               \
    {                                                     \
	type src_pixel, *dest_pixels;                     \
	dest_pixels = (type *) VPixelPtr (result,         \
					  dest_curband,   \
					  0,              \
					  0);             \
	for (r = 0; r < dest_nrows; r++) {                \
            for (c = 0; c < dest_ncols; c++) {            \
                src_pixel = VPixel (src,                  \
				    src_curband,          \
				    c,                    \
				    r,                    \
				    type);                \
	        *dest_pixels++ = src_pixel;               \
	    }                                             \
	}                                                 \
    }
    
  /* For each band in the dest image do: */
  for (b = 0; b < dest_nbands; b++) {

    src_curband = (band == VAllBands) ? b : band;
    dest_curband = (band == VAllBands) ? b : 0;
	
    /* Perform transposition according to pixel representation: */
    switch (src_repn) {
    case VBitRepn:          Transpose(VBit);    break;
    case VUByteRepn:	Transpose(VUByte);  break;
    case VSByteRepn:	Transpose(VSByte);  break;
    case VShortRepn:	Transpose(VShort);  break;
    case VUShortRepn:	Transpose(VUShort);  break;
    case VIntegerRepn:	Transpose(VInteger);  break;
    case VUIntegerRepn:	Transpose(VUInteger);  break;
    case VLongRepn:	Transpose(VLong);   break;
    case VULongRepn:	Transpose(VULong);   break;
    case VFloatRepn:	Transpose(VFloat);  break;
    case VDoubleRepn:	Transpose(VDouble); break;
    default: break;
    }
  }

  if (src == dest) {
    VCopyImagePixels (result, dest, VAllBands);
    VDestroyImage (result);
    return dest;
  } else {
    VCopyImageAttrs (src, result);
    return result;
  }

#undef Transpose
}
