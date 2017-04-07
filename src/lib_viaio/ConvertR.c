/*
 *  $Id: ConvertR.c 3177 2008-04-01 14:47:24Z karstenm $
 *
 *  This file contains a routines for converting an image from one pixel
 *  representation to another by mapping the full range of source pixel
 *  values to the full range of destination pixel values.
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
#include "viaio/mu.h"

/* From the standard C library: */
#include <math.h>

/* File identification string: */
VRcsId ("$Id: ConvertR.c 3177 2008-04-01 14:47:24Z karstenm $");

extern double rint(double);

/*
 *  Some macros for converting one type to another.
 */

#define FillTable(table, op)						\
    for (i = 0; i < VNumber (table); i++)				\
        table[i] = op;

#define Convert(src_type, result_type, op)				\
    {									\
	src_type *src_pp = (src_type *) src_first;			\
	result_type *res_pp = VImageData (result);			\
									\
	for (i = 0; i < npixels; i++) {					\
	    op;								\
	}								\
    }

#define Clip								\
    if (t < 0) {							\
	t = 0;								\
	clipped = TRUE;							\
    }

#define Bracket(a,b)							\
    if (t < a) {							\
	t = a;								\
	clipped = TRUE;							\
    } else if (t > b) {							\
	t = b;								\
	clipped = TRUE;							\
    }

/*
 *  VConvertImageRange
 *
 *  Converts an image from one pixel representation to another using a
 *  mapping from all possible source pixel values to all destination pixel
 *  values.
 *
 *  Note: In ANSI C, left shifting a signed, negative value may or may
 *	  produce a negative result. Where we now the value's positive,
 *	  we use >>; where it could be negative, we use /.
 */

VImage VConvertImageRange (VImage src, VImage dest, VBand band,
			   VRepnKind pixel_repn)
{
    VImage result;
    int npixels, i;
    double v;
    VPointer src_first;
    VBoolean clipped = FALSE;

    /* If src already has the requested representation, simply copy: */
    if (pixel_repn == VPixelRepn (src))
	return (src == dest) ? src : VCopyImage (src, dest, band);
    
    /* Prepare to iterate over all source pixels: */
    if (! VSelectBand ("VConvertImageRange", src, band,
		       & npixels, & src_first))
	return NULL;

    /* Check dest if it exists; create it if it doesn't: */
    result = VSelectDestImage ("VConvertImageRange", dest,
			       band == VAllBands ? VImageNBands (src) : 1,
			       VImageNRows (src), VImageNColumns (src),
			       pixel_repn);
    if (! result)
	return NULL;

    /* Map pixels from one representation to the other: */
    switch (VPixelRepn (src)) {
	
    case VBitRepn:
	switch (VPixelRepn (result)) {
	    
	case VUByteRepn: 
	    Convert (VBit, VUByte, *res_pp++ = *src_pp++ ? 255 : 0); 
	    break;
	    
	case VSByteRepn:
	    Convert (VBit, VSByte, *res_pp++ = *src_pp++ ? 127 : 0);
	    break;
	    
	case VShortRepn:
	    Convert (VBit, VShort, *res_pp++ = *src_pp++ ? 32767 : 0);
	    break;
	    
	case VLongRepn:
	    Convert (VBit, VLong, *res_pp++ = *src_pp++ ? 0x7FFFFFFFl : 0);
	    break;
	    
	case VFloatRepn:
	    Convert (VBit, VFloat,
		     *res_pp++ = (*src_pp++ ?
				  VFloatConst (1.0) : VFloatConst (0.0)));
	    break;
	        
	case VDoubleRepn:
	    Convert (VBit, VDouble,
		     *res_pp++ = (*src_pp++ ?
				  VDoubleConst (1.0) : VDoubleConst (0.0)));
	    break;

	default:
	    break;
	}
	break;
	
    case VUByteRepn:
	switch (VPixelRepn (result)) {
	    
	case VBitRepn:
	    Convert (VUByte, VBit, *res_pp++ = *src_pp++ >> 7);
	    break;
	    
	case VSByteRepn:
	    Convert (VUByte, VSByte, *res_pp++ = *src_pp++ >> 1);
	    break;
	    	    
	case VShortRepn:
	    Convert (VUByte, VShort, *res_pp++ = *src_pp++ << 7);
	    break;
	    
	case VLongRepn:
	    Convert (VUByte, VLong, *res_pp++ = *src_pp++ << 23);
	    break;
	    
	case VFloatRepn:
	    {
		VFloat table[256];
		FillTable (table, (i / VFloatConst (255.0)));
		Convert (VUByte, VFloat, *res_pp++ = table[*src_pp++]);
	    }
	    break;
	    
	case VDoubleRepn:
	    {
		VDouble table[256];
		FillTable (table, (i / VDoubleConst (255.0)));
		Convert (VUByte, VDouble, *res_pp++ = table[*src_pp++]);
	    }
	    break;

	default:
	    break;
	}
	break;
	
    case VSByteRepn:
	switch (VPixelRepn (result)) {
	    
	case VBitRepn:
	    {
		VSByte t;
		Convert (VSByte, VBit,
			 t = *src_pp++; Clip; *res_pp++ = t >> 6);
	    }
	    break;
	    
	case VUByteRepn:
	    {
		VSByte t;
		Convert (VSByte, VUByte,
			 t = *src_pp++; Clip; *res_pp++ = t << 1);
	    }
	    break;

	case VShortRepn:
	    Convert (VSByte, VShort, *res_pp++ = ((VShort) *src_pp++ << 8));
	    break;
	    
	case VLongRepn:
	    Convert (VSByte, VLong, *res_pp++ = ((VLong) *src_pp++ << 24));
	    break;
	    
	case VFloatRepn:
	    {
		VFloat table[256];
		FillTable (table,
			   (i - VFloatConst (128.0)) / VFloatConst (128.0));
		Convert (VSByte, VFloat, *res_pp++ = table[*src_pp++ + 128]);
	    }
	    break;
	    
	case VDoubleRepn:
	    {
		VDouble table[256];
		FillTable (table,
			   (i - VDoubleConst (128)) / VDoubleConst (128.0));
		Convert (VSByte, VDouble, *res_pp++ = table[*src_pp++ + 128]);
	    }
	    break;

	default:
	    break;
	}
	break;
	
    case VShortRepn:
	v = -VPixelMinValue (src);
	switch (VPixelRepn (result)) {
	    
	case VBitRepn:
	    {
		VShort t;
		Convert (VShort, VBit,
			 t = *src_pp++; Clip; *res_pp++ = t >> 14);
	    }
	    break;
	    
	case VUByteRepn:
	    {
		VShort t;
		Convert (VShort, VUByte,
			 t = *src_pp++; Clip; *res_pp++ = t >> 7);
	    }
	    break;

	case VSByteRepn:
	    Convert (VShort, VSByte, *res_pp++ = (*src_pp++ / 256));
	    break;
	    
	case VLongRepn:
	    Convert (VShort, VLong, *res_pp++ = ((VLong) *src_pp++ << 8));
	    break;
	    
	case VFloatRepn:
	    Convert (VShort, VFloat, *res_pp++ = *src_pp++ / v);
	    break;
	    
	case VDoubleRepn:
	    Convert (VShort, VDouble, *res_pp++ = *src_pp++ / v);
	    break;

	default:
	    break;
	}
	break;
	
    case VLongRepn:
	v = -VPixelMinValue (src);
	switch (VPixelRepn (result)) {
	    
	case VBitRepn:
	    {
		VLong t;
		Convert (VLong, VBit,
			 t = *src_pp++; Clip; *res_pp++ = t >> 30);
	    }
	    break;
	    
	case VUByteRepn:
	    {
		VLong t;
		Convert (VLong, VUByte,
			 t = *src_pp++; Clip; *res_pp++ = t >> 23);
	    }
	    break;

	case VSByteRepn:
	    Convert (VLong, VSByte, *res_pp++ = (*src_pp++ / (1 << 23)));
	    break;
	    
	case VShortRepn:
	    Convert (VLong, VShort,
		     *res_pp++ = ((VLong) *src_pp++ / (1 << 16)));
	    break;
	    
	case VFloatRepn:
	    Convert (VLong, VFloat, *res_pp++ = *src_pp++ / v);
	    break;
	    
	case VDoubleRepn:
	    Convert (VLong, VDouble, *res_pp++ = *src_pp++ / v);
	    break;

	default:
	    break;
	}
	break;
	
    case VFloatRepn:
        {
	    VFloat t, m = VPixelMaxValue (result);
	    switch (VPixelRepn (result)) {
	    
	    case VBitRepn:
		Convert (VFloat, VBit,
			 t = *src_pp++; Bracket (0, 1);
			 *res_pp++ = t >= 0.5);
		break;
	    
	    case VUByteRepn:
		Convert (VFloat, VUByte,
			 t = *src_pp++; Bracket (0, 1);
			 *res_pp++ = rint (t * m));
		break;
	    
	    case VSByteRepn:
		Convert (VFloat, VSByte,
			 t = *src_pp++; Bracket (-1, 1);
			 *res_pp++ = rint (t * m));
		break;
	    
	    case VShortRepn:
		Convert (VFloat, VShort,
			 t = *src_pp++; Bracket (-1, 1);
			 *res_pp++ = rint (t * m));
		break;
	    
	    case VLongRepn:
		Convert (VFloat, VLong,
			 t = *src_pp++; Bracket (-1, 1);
			 *res_pp++ = rint (t * m));
		break;
	    
	    case VDoubleRepn:
		Convert (VFloat, VDouble, *res_pp++ = *src_pp++);
		break;

	    default:
		break;
	    }
	}
	break;
	
    case VDoubleRepn:
        {
	    VDouble t, m = VPixelMaxValue (result);
	    switch (VPixelRepn (result)) {
	    
	    case VBitRepn:
		Convert (VDouble, VBit,
			 t = *src_pp++; Bracket (0, 1);
			 *res_pp++ = t >= 0.5);
		break;
	    
	    case VUByteRepn:
		Convert (VDouble, VUByte,
			 t = *src_pp++; Bracket (0, 1);
			 *res_pp++ = rint (t * m));
		break;
	    
	    case VSByteRepn:
		Convert (VDouble, VSByte,
			 t = *src_pp++; Bracket (-1, 1);
			 *res_pp++ = rint (t * m));
		break;
	    
	    case VShortRepn:
		Convert (VDouble, VShort,
			 t = *src_pp++; Bracket (-1, 1);
			 *res_pp++ = rint (t * m));
		break;
	    
	    case VLongRepn:
		Convert (VDouble, VLong,
			 t = *src_pp++; Bracket (-1, 1);
			 *res_pp++ = rint (t * m));
		break;
	    
	    case VFloatRepn:
		Convert (VDouble, VFloat, *res_pp++ = *src_pp++);
		break;

	    default:
		break;
	    }

	}
	break;
	
    default:
	break;
    }

    if (clipped)
	VWarning ("VConvertImageRange: "
		  "Pixel value exceeds destination range");

    VCopyImageAttrs (src, result);

    return result;
}
