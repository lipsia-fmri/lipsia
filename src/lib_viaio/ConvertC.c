/*
 *  $Id: ConvertC.c 3177 2008-04-01 14:47:24Z karstenm $
 *
 *  This file contains a routine for converting an image from one pixel
 *  representation to another by simply copying (not mapping) pixel values.
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

/* From the standard C library: */
#include <math.h>

/* File identification string: */
VRcsId ("$Id: ConvertC.c 3177 2008-04-01 14:47:24Z karstenm $");

extern double rint(double);

/*
 *  Some macros for converting one type to another.
 */

#define Cast(src_type, result_type, op)					\
    {									\
	src_type pixel, *src_pp = (src_type *) src_first;		\
	result_type *res_pp = VImageData (result);			\
									\
	for (i = 0; i < npixels; i++) {					\
	    pixel = *src_pp++;						\
	    op;								\
	    *res_pp++ = pixel;						\
	}								\
    }

#define Clip(test, limit) if (pixel test limit) { pixel = limit; clipped = TRUE; }

#define Bracket(lower, upper) Clip (<, lower); Clip (>, upper);

#define Nothing

/*
 *  VConvertImageCopy
 *
 *  Converts an image from one pixel representation to another by copying
 *  (and perhaps rounding) pixel values, not mapping them from one range
 *  to another.
 *
 *  Note: This uses rint() to round a floating point value to the nearest
 *	  integer (under the default rounding mode). It isn't part of
 *	  ANSI C, but it's widely available.
 */

VImage VConvertImageCopy (VImage src, VImage dest, VBand band,
			  VRepnKind pixel_repn)
{
    VImage result;
    int npixels, i;
    VPointer src_first;
    VBoolean clipped = FALSE;

    /* If src already has the requested representation, simply copy: */
    if (pixel_repn == VPixelRepn (src))
	return (src == dest) ? src : VCopyImage (src, dest, band);
    
    /* Prepare to iterate over all source pixels: */
    if (! VSelectBand ("VConvertImageCopy", src, band, & npixels, & src_first))
	return NULL;

    /* Check dest if it exists; create it if it doesn't: */
    result = VSelectDestImage ("VConvertImageCopy", dest,
			       band == VAllBands ? VImageNBands (src) : 1,
			       VImageNRows (src), VImageNColumns (src),
			       pixel_repn);
    if (! result)
	return NULL;

    /* Copy pixels, casting or rounding them: */
    switch (VPixelRepn (src)) {
	
    case VBitRepn:
	switch (VPixelRepn (result)) {
	    
	case VUByteRepn: Cast (VBit, VUByte, Nothing); break;
	    
	case VSByteRepn: Cast (VBit, VSByte, Nothing); break;
	    
	case VShortRepn: Cast (VBit, VShort, Nothing); break;
	    
	case VLongRepn: Cast (VBit, VLong, Nothing); break;
	    
	case VFloatRepn: Cast (VBit, VFloat, Nothing); break;
	    
	case VDoubleRepn: Cast (VBit, VDouble, Nothing); break;

	default: break;
	}
	break;
	
    case VUByteRepn:
	switch (VPixelRepn (result)) {
	    
	case VBitRepn: Cast (VUByte, VBit, Clip (>, 1)); break;
	    
	case VSByteRepn: Cast (VUByte, VSByte, Clip (>, 127)); break;
	    
	case VShortRepn: Cast (VUByte, VShort, Nothing); break;
	    
	case VLongRepn: Cast (VUByte, VLong, Nothing); break;
	    
	case VFloatRepn: Cast (VUByte, VFloat, Nothing); break;
	    
	case VDoubleRepn: Cast (VUByte, VDouble, Nothing); break;

	default: break;
	}
	break;
	
    case VSByteRepn:
	switch (VPixelRepn (result)) {
	    
	case VBitRepn: Cast (VSByte, VBit, Bracket (0, 1)); break;
	    
	case VUByteRepn: Cast (VSByte, VUByte, Clip (<, 0)); break;
	    
	case VShortRepn: Cast (VSByte, VShort, Nothing); break;
	    
	case VLongRepn: Cast (VSByte, VLong, Nothing); break;
	    
	case VFloatRepn: Cast (VSByte, VFloat, Nothing); break;
	    
	case VDoubleRepn: Cast (VSByte, VDouble, Nothing); break;

	default: break;
	}
	break;
	
    case VShortRepn:
	switch (VPixelRepn (result)) {
	    
	case VBitRepn: Cast (VShort, VBit, Bracket (0, 1)); break;
	    
	case VUByteRepn: Cast (VShort, VUByte, Bracket (0, 255)); break;
	    
	case VSByteRepn: Cast (VShort, VSByte, Bracket (-128, 127)); break;
	    
	case VLongRepn: Cast (VShort, VLong, Nothing); break;
	    
	case VFloatRepn: Cast (VShort, VFloat, Nothing); break;
	    
	case VDoubleRepn: Cast (VShort, VDouble, Nothing); break;

	default: break;
	}
	break;
	
    case VLongRepn:
	switch (VPixelRepn (result)) {
	    
	case VBitRepn: Cast (VLong, VBit, Bracket (0, 1)); break;
	    
	case VUByteRepn: Cast (VLong, VUByte, Bracket (0, 255)); break;
	    
	case VSByteRepn: Cast (VLong, VSByte, Bracket (-128, 127)); break;
	    
	case VShortRepn: Cast (VLong, VShort, Bracket (-32768, 32767)); break;
	    
	case VFloatRepn: Cast (VLong, VFloat, Nothing); break;
	    
	case VDoubleRepn: Cast (VLong, VDouble, Nothing); break;

	default: break;
	}
	break;
	
    case VFloatRepn:
	switch (VPixelRepn (result)) {
	    
	case VBitRepn:
	    Cast (VFloat, VBit, Bracket (0, 1); pixel = rint (pixel));
	    break;
	    
	case VUByteRepn:
	    Cast (VFloat, VUByte,
		  Bracket (0, 255);  pixel = rint (pixel));
	    break;
	    
	case VSByteRepn:
	    Cast (VFloat, VSByte,
		  Bracket (-128, 127);  pixel = rint (pixel));
	    break;
	    
	case VShortRepn:
	    Cast (VFloat, VShort,
		  Bracket (-32768, 32767); pixel = rint (pixel));
	    break;
	    
	case VLongRepn:
	    Cast (VFloat, VLong,
		  Bracket (VRepnMinValue (VLongRepn),
			   VRepnMaxValue (VLongRepn)); );
	    break;
	    
	case VDoubleRepn: Cast (VFloat, VDouble, Nothing); break;

	default: break;
	}
	break;
	
    case VDoubleRepn:
	switch (VPixelRepn (result)) {
	    
	case VBitRepn:
	    Cast (VDouble, VBit, Bracket (0, 1);  pixel = rint (pixel));
	    break;
	    
	case VUByteRepn:
	    Cast (VDouble, VUByte, Bracket (0, 255);  pixel = rint (pixel));
	    break;
	    
	case VSByteRepn:
	    Cast (VDouble, VSByte, Bracket (-128, 127);  pixel = rint (pixel));
	    break;
	    
	case VShortRepn:
	    Cast (VDouble, VShort, Bracket (-32768, 32767);  pixel = rint (pixel));
	    break;
	    
	case VLongRepn:
	    Cast (VDouble, VLong,
		  Bracket (VRepnMinValue (VLongRepn),
			   VRepnMaxValue (VLongRepn)));
	    break;
	    
	case VFloatRepn: Cast (VDouble, VFloat, Nothing); break;

	default: break;
	}

    default: break;
    }

    if (clipped)
	VWarning ("VConvertImageCopy: "
		  "Source pixel value exceeds destination range");

    VCopyImageAttrs (src, result);

    return result;
}
