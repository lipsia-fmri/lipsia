/*
 *  $Id: Fill.c 3177 2008-04-01 14:47:24Z karstenm $
 *
 *  This file contains a routine for filling an image with a pixel
 *  value constant.
 */

/*
 *  Copyright 1994 University of British Columbia
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
VRcsId ("$Id: Fill.c 3177 2008-04-01 14:47:24Z karstenm $");


/*
 *  VFillImage
 *
 *  Fills one or all image bands with a specified pixel value.
 *
 *  ASSUMPTION: Floating point 0.0 is binary all-bits-zero.
 *
 *  This is true of the IEEE standard floating point format, but it
 *  isn't assured by the ANSI C standard.
 */

VBoolean VFillImage (VImage image, VBand band, VDoublePromoted value)
{
  int32_t i, npixels;
  VPointer first_pixel;

#define Fill(type)					\
	{						\
	    type d = value, *pp = first_pixel;		\
	    for (i = npixels; i > 0; i--)		\
		*pp++ = d;				\
	}

  /* Locate the specified band(s) in the image: */
  if (! VSelectBand ("VFillImage", image, band, & npixels, & first_pixel))
    return FALSE;

  /* To avoid surprises when filling integer pixels, round the fill
     value to the nearest integer: */
  if (VIsIntegerRepn (VPixelRepn (image)))
    value += (value > 0.0) ? 0.5 : -0.5;

  /* The memset() routine is probably the fastest way of filling memory, but
     it can only be used here in certain situations. We only check for
     these, the most straightforward of the situations:
     (a) pixel values occupy a byte
     (b) the value to be filled is all 0's
     It is when the value to be filled is all 0's and the pixel
     representation is floating point that we take advantage of the
     assumption that 0.0 is represented as all-bits-zero. */
  if (VPixelSize (image) == 1 || value == 0.0)
    memset (first_pixel, (int) value, npixels * VPixelSize (image));

  /* Otherwise, fill by looping over all pixels: */
  else
    switch (VPixelRepn (image)) {
    case VBitRepn:		Fill (VBit);	break;
    case VUByteRepn:	Fill (VUByte);  break;
    case VSByteRepn:	Fill (VSByte);  break;
    case VShortRepn:	Fill (VShort);  break;
    case VFloatRepn:	Fill (VFloat);  break;
    case VDoubleRepn:	Fill (VDouble); break;
    case VUShortRepn:	Fill (VUShort);  break;
    case VIntegerRepn:  Fill (VInteger);   break;
    case VUIntegerRepn: Fill (VUInteger);   break;
    case VLongRepn:     Fill (VLong);   break;
    case VULongRepn:    Fill (VULong);   break;
    default: break;
    }

  return TRUE;

#undef Fill
}
