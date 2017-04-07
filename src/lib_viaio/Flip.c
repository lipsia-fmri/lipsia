/*
 *  $Id: Flip.c 3177 2008-04-01 14:47:24Z karstenm $
 *
 *  This file contains routines for flipping images
 *
 * 10/09/1997: bug in horizonatal flipping of VShort images corrected (FK)
 *
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
 * VFlipImage
 *
 * Flip an image horizontally or vertically about its centerline.
 */

VImage VFlipImage (VImage src, VImage dest, VBand band,
		   VBooleanPromoted vertical)
{
  int dest_nbands, src_band, res_band, b, r, c, row_size;
  VImage result;

  /* Flip an image horizontally: */
#define Flip(type)                                   		\
    {									\
	for (r = 0; r < VImageNRows (src); r++) {			\
            type *src_pp = VPixelPtr (src, src_band, r, VImageNColumns(src)-1); \
            type *res_pp = VPixelPtr (result, res_band, r, 0);		\
            for (c = 0; c < VImageNColumns (src); c++)			\
            	*res_pp++ = *src_pp--; 					\
	}								\
    }
    
  /* Check to ensure that "band" exists: */
  if (band != VAllBands && (band < 0 || band >= VImageNBands (src))) {
    VWarning ("VFlipImage: Band %d referenced in image of %d bands",
	      band, VImageNBands (src));
    return NULL;
  }
    
  /* Create/select result image: */
  dest_nbands = (band == VAllBands) ? VImageNBands (src) : 1;
  result = VSelectDestImage ("VFlipImage", dest, dest_nbands,
			     VImageNRows (src), VImageNColumns (src),
			     VPixelRepn (src));
  if (! result)
    return NULL;	/* Cannot create destination image */

  if (src == dest) {
    result = VCreateImage (dest_nbands, VImageNRows (src),
			   VImageNColumns (src), VPixelRepn (src));
    if (! result)
      return NULL;
  }

  row_size = VImageNColumns (src) * VPixelSize (src);

  /* For each band in the dest image... */
  for (b = 0; b < dest_nbands; b++) {
	
    src_band = (band == VAllBands) ? b : band;
    res_band = (band == VAllBands) ? b : 0;
	
    if (vertical) {
      char *src_pp, *res_pp;

      /* Flip vertically by copying entire rows: */
      src_pp = VPixelPtr (src, src_band, VImageNRows (src) - 1, 0);
      res_pp = VPixelPtr (result, res_band, 0, 0);
      for (r = 0; r < VImageNRows (src); r++) {
	memcpy (res_pp, src_pp, row_size);
	src_pp -= row_size;
	res_pp += row_size;
      }

    } else {

      /* Flip horizontally by copying pixels: */
      switch (VPixelRepn (src)) {
      case VBitRepn:	Flip (VBit);	break;
      case VUByteRepn:	Flip (VUByte);	break;
      case VSByteRepn:	Flip (VSByte);	break;
      case VShortRepn:	Flip (VShort);	break;
      case VIntegerRepn: Flip (VInteger); break;
      case VUIntegerRepn: Flip (VUInteger); break;
      case VLongRepn:	Flip (VLong);	break;
      case VULongRepn:	Flip (VULong);	break;
      case VFloatRepn:	Flip (VFloat);	break;
      case VDoubleRepn:	Flip (VDouble);	break;
      default: break;
      }
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

#undef Flip
}

