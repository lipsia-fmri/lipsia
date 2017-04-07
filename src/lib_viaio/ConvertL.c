/*
 *  $Id: ConvertL.c 3177 2008-04-01 14:47:24Z karstenm $
 *
 *  This file contains a routines for converting an image from one pixel
 *  representation to another by an abitrary lineary mapping
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
VRcsId ("$Id: ConvertL.c 3177 2008-04-01 14:47:24Z karstenm $");

extern double rint(double);

/*
 *  Some macros for converting one type to another.
 */

#define FillTable(table, op)						\
    for (i = 0; i < VNumber (table); i++) {				\
        op;								\
    }

#define Convert(src_type, result_type, op)				\
    {									\
	src_type *src_pp = (src_type *) src_first;			\
	result_type *res_pp = VImageData (result);			\
									\
	for (i = 0; i < npixels; i++) {					\
	    op;								\
	}								\
    }

#define Bracket(t, lower, upper)					\
    if (t < lower) {							\
	t = lower;							\
	clipped = TRUE;							\
    } else if (t > upper) {						\
	t = upper;							\
	clipped = TRUE;							\
    }


/*
 *  VConvertImageLinear
 *
 *  Converts an image from one pixel representation to another using a
 *  linear mapping between source and destination pixel values:
 *
 *	dest_pixel = src_pixel * a + b
 *
 *  for a pair a, b supplied as parameters.
 */

VImage VConvertImageLinear (VImage src, VImage dest, VBand band,
			    VRepnKind pixel_repn, double a, double b)
{
    VImage result;
    int npixels, i;
    VPointer src_first;
    VDouble a_min, a_max, d0, d1, t;
    VDouble d_min =
	VIsFloatPtRepn (pixel_repn) ? -1.0 : VRepnMinValue (pixel_repn);
    VDouble d_max = 
	VIsFloatPtRepn (pixel_repn) ? 1.0 : VRepnMaxValue (pixel_repn);
    VBoolean clipped = FALSE;

    /* Prepare to iterate over all source pixels: */
    if (! VSelectBand ("VConvertImageLinear",
		       src, band, & npixels, & src_first))
	return NULL;

    /* Check dest if it exists; create it if it doesn't: */
    result = VSelectDestImage ("VConvertImageLinear", dest,
			       band == VAllBands ? VImageNBands (src) : 1,
			       VImageNRows (src), VImageNColumns (src),
			       pixel_repn);
    if (! result)
	return NULL;

    /* Map pixels from one representation to the other: */
    switch (VPixelRepn (src)) {

    case VBitRepn:
	d0 = b;
	d1 = a + b;
	if (pixel_repn != VFloatRepn && pixel_repn != VDoubleRepn) {
	    Bracket (d0, d_min, d_max);
	    Bracket (d1, d_min, d_max);
	}
	switch (VPixelRepn (result)) {

	case VBitRepn:
	    {
		VBit dd0 = rint (d0), dd1 = rint (d1);
		Convert (VBit, VBit, *res_pp++ = *src_pp++ ? dd1 : dd0);
	    }
	    break;

	case VUByteRepn:
	    {
		VUByte dd0 = rint (d0), dd1= rint (d1);
		Convert (VBit, VUByte, *res_pp++ = *src_pp++ ? dd1 : dd0);
	    }
	    break;

	case VSByteRepn:
	    {
		VSByte dd0 = rint (d0), dd1= rint (d1);
		Convert (VBit, VUByte, *res_pp++ = *src_pp++ ? dd1 : dd0);
	    }
	    break;

	case VShortRepn:
	    {
		VShort dd0 = rint (d0), dd1= rint (d1);
		Convert (VBit, VUByte, *res_pp++ = *src_pp++ ? dd1 : dd0);
	    }
	    break;

	case VLongRepn:
	    {
		VLong dd0 = rint (d0), dd1= rint (d1);
		Convert (VBit, VUByte, *res_pp++ = *src_pp++ ? dd1 : dd0);
	    }
	    break;

	case VFloatRepn:
	    {
		Convert (VBit, VFloat, *res_pp++ = *src_pp++ ? d1 : d0);
	    }
	    break;

	case VDoubleRepn:
	    {
		Convert (VBit, VDouble, *res_pp++ = *src_pp++ ? d1 : d0);
	    }
	    break;

	default:
	    break;
	}
	break;

    case VUByteRepn:
	switch (VPixelRepn (result)) {

	case VBitRepn:
	    {
		VBit table[256];
		FillTable (table, t = i * a + b;
			   Bracket (t, d_min, d_max); table[i] = rint (t));
		Convert (VUByte, VBit, *res_pp++ = table[*src_pp++]);
	    }
	    break;

	case VUByteRepn:
	    {
		VUByte table[256];
		FillTable (table, t = i * a + b;
			   Bracket (t, d_min, d_max); table[i] = rint (t));
		Convert (VUByte, VUByte, *res_pp++ = table[*src_pp++]);
	    }
	    break;

	case VSByteRepn:
	    {
		VSByte table[256];
		FillTable (table, t = i * a + b;
			   Bracket (t, d_min, d_max); table[i] = rint (t));
		Convert (VUByte, VSByte, *res_pp++ = table[*src_pp++]);
	    }
	    break;

	case VShortRepn:
	    {
		VShort table[256];
		FillTable (table, t = i * a + b;
			   Bracket (t, d_min, d_max); table[i] = rint (t));
		Convert (VUByte, VShort, *res_pp++ = table[*src_pp++]);
	    }
	    break;

	case VLongRepn:
	    {
		VLong table[256];
		FillTable (table, t = i * a + b;
			   Bracket (t, d_min, d_max); table[i] = rint (t));
		Convert (VUByte, VLong, *res_pp++ = table[*src_pp++]);
	    }
	    break;

	case VFloatRepn:
	    {
		VFloat table[256];
		FillTable (table, table[i] = i * a + b);
		Convert (VUByte, VFloat, *res_pp++ = table[*src_pp++]);
	    }
	    break;

	case VDoubleRepn:
	    {
		VDouble table[256];
		FillTable (table, table[i] = i * a + b);
		Convert (VUByte, VFloat, *res_pp++ = table[*src_pp++]);
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
		VBit table[256];
		FillTable (table, t = (i - 128) * a + b;
			   Bracket (t, d_min, d_max); table[i] = rint (t));
		Convert (VSByte, VBit, *res_pp++ = table[*src_pp++ + 128]);
	    }
	    break;

	case VUByteRepn:
	    {
		VUByte table[256];
		FillTable (table, t = (i - 128) * a + b;
			   Bracket (t, d_min, d_max); table[i] = rint (t));
		Convert (VSByte, VUByte, *res_pp++ = table[*src_pp++ + 128]);
	    }
	    break;

	case VSByteRepn:
	    {
		VSByte table[256];
		FillTable (table, t = (i - 128) * a + b;
			   Bracket (t, d_min, d_max); table[i] = rint (t));
		Convert (VSByte, VSByte, *res_pp++ = table[*src_pp++ + 128]);
	    }
	    break;

	case VShortRepn:
	    {
		VShort table[256];
		FillTable (table, t = (i - 128) * a + b;
			   Bracket (t, d_min, d_max); table[i] = rint (t));
		Convert (VSByte, VShort, *res_pp++ = table[*src_pp++ + 128]);
	    }
	    break;

	case VLongRepn:
	    {
		VLong table[256];
		FillTable (table, t = (i - 128) * a + b;
			   Bracket (t, d_min, d_max); table[i] = rint (t));
		Convert (VSByte, VLong, *res_pp++ = table[*src_pp++ + 128]);
	    }
	    break;

	case VFloatRepn:
	    {
		VFloat table[256];
		FillTable (table, table[i] = (i - 128) * a + b);
		Convert (VSByte, VFloat, *res_pp++ = table[*src_pp++ + 128]);
	    }
	    break;

	case VDoubleRepn:
	    {
		VDouble table[256];
		FillTable (table, table[i] = (i - 128) * a + b);
		Convert (VSByte, VFloat, *res_pp++ = table[*src_pp++ + 128]);
	    }
	    break;

	default:
	    break;
	}
	break;

    case VShortRepn:
	switch (VPixelRepn (result)) {

	case VBitRepn:
	    Convert (VShort, VBit, t = *src_pp++ * a + b;
		     Bracket (t, d_min, d_max); *res_pp++ = rint (t));
	    break;

	case VUByteRepn:
	    Convert (VShort, VUByte, t = *src_pp++ * a + b;
		     Bracket (t, d_min, d_max); *res_pp++ = rint (t));
	    break;

	case VSByteRepn:
	    Convert (VShort, VSByte, t = *src_pp++ * a + b;
		     Bracket (t, d_min, d_max); *res_pp++ = rint (t));
	    break;

	case VShortRepn:
	    Convert (VShort, VShort, t = *src_pp++ * a + b;
		     Bracket (t, d_min, d_max); *res_pp++ = rint (t));
	    break;

	case VLongRepn:
	    Convert (VShort, VLong, t = *src_pp++ * a + b;
		     Bracket (t, d_min, d_max); *res_pp++ = rint (t));
	    break;

	case VFloatRepn:
	    Convert (VShort, VFloat, *res_pp++ = *src_pp++ * a + b);
	    break;

	case VDoubleRepn:
	    Convert (VShort, VDouble, *res_pp++ = *src_pp++ * a + b);
	    break;

	default:
	    break;
	}
	break;

    case VLongRepn:
	switch (VPixelRepn (result)) {

	case VBitRepn:
	    Convert (VLong, VBit, t = *src_pp++ * a + b;
		     Bracket (t, d_min, d_max); *res_pp++ = rint (t));
	    break;

	case VUByteRepn:
	    Convert (VLong, VUByte, t = *src_pp++ * a + b;
		     Bracket (t, d_min, d_max); *res_pp++ = rint (t));
	    break;

	case VSByteRepn:
	    Convert (VLong, VSByte, t = *src_pp++ * a + b;
		     Bracket (t, d_min, d_max); *res_pp++ = rint (t));
	    break;

	case VShortRepn:
	    Convert (VLong, VShort, t = *src_pp++ * a + b;
		     Bracket (t, d_min, d_max); *res_pp++ = rint (t));
	    break;

	case VLongRepn:
	    Convert (VLong, VLong, t = *src_pp++ * a + b;
		     Bracket (t, d_min, d_max); *res_pp++ = rint (t));
	    break;

	case VFloatRepn:
	    Convert (VLong, VFloat, *res_pp++ = *src_pp++ * a + b);
	    break;

	case VDoubleRepn:
	    Convert (VLong, VDouble, *res_pp++ = *src_pp++ * a + b);
	    break;

	default:
	    break;
	}
	break;

    case VFloatRepn:
	switch (VPixelRepn (result)) {

	case VBitRepn:
	    Convert (VFloat, VBit, t = *src_pp++ * a + b;
		     Bracket (t, d_min, d_max); *res_pp++ = rint (t));
	    break;

	case VUByteRepn:
	    Convert (VFloat, VUByte, t = *src_pp++ * a + b;
		     Bracket (t, d_min, d_max); *res_pp++ = rint (t));
	    break;

	case VSByteRepn:
	    Convert (VFloat, VSByte, t = *src_pp++ * a + b;
		     Bracket (t, d_min, d_max); *res_pp++ = rint (t));
	    break;

	case VShortRepn:
	    Convert (VFloat, VShort,  t = *src_pp++ * a + b;
		     Bracket (t, d_min, d_max); *res_pp++ = rint (t));
	    break;

	case VLongRepn:
	    Convert (VFloat, VLong, t = *src_pp++ * a + b;
		     Bracket (t, d_min, d_max); *res_pp++ = rint (t));
	    break;

	case VFloatRepn:
	    Convert (VFloat, VFloat, *res_pp++ = *src_pp++ * a + b);
	    break;

	case VDoubleRepn:
	    Convert (VFloat, VDouble, *res_pp++ = *src_pp++ * a + b);
	    break;

	default:
	    break;
	}
	break;

    case VDoubleRepn:
	switch (VPixelRepn (result)) {

	case VBitRepn:
	    Convert (VDouble, VBit, t = *src_pp++ * a + b;
		     Bracket (t, d_min, d_max); *res_pp++ = rint (t));
	    break;

	case VUByteRepn:
	    Convert (VDouble, VUByte, t = *src_pp++ * a + b;
		     Bracket (t, d_min, d_max); *res_pp++ = rint (t));
	    break;

	case VSByteRepn:
	    Convert (VDouble, VSByte, t = *src_pp++ * a + b;
		     Bracket (t, d_min, d_max); *res_pp++ = rint (t));
	    break;

	case VShortRepn:
	    Convert (VDouble, VShort, t = *src_pp++ * a + b;
		     Bracket (t, d_min, d_max); *res_pp++ = rint (t));
	    break;

	case VLongRepn:
	    Convert (VDouble, VLong, t = *src_pp++ * a + b;
		     Bracket (t, d_min, d_max); *res_pp++ = rint (t));
	    break;

	case VFloatRepn:
	    Convert (VDouble, VFloat, *res_pp++ = *src_pp++ * a + b);
	    break;

	case VDoubleRepn:
	    Convert (VDouble, VDouble, *res_pp++ = *src_pp++ * a + b);
	    break;

	default:
	    break;
	}

    default:
	break;
    }

    /* In the cases were a table was used for the conversion, we don't know
       for sure whether any pixel values were clipped. This determines
       with certainty: */
    if (clipped && (pixel_repn == VBitRepn ||
		    pixel_repn == VUByteRepn || pixel_repn == VSByteRepn)) {
	clipped = FALSE;
	VImageStats (src, band, & a_min, & a_max, NULL, NULL);
	t = a_min * a + b;
	if (t < d_min || t > d_max)
	    clipped = TRUE;
	t = a_max * a + b;
	if (t < d_min || t > d_max)
	    clipped = TRUE;
    }	   

    if (clipped)
	VWarning ("VConvertImageLinear: "
		  "Pixel value exceeds destination range");

    VCopyImageAttrs (src, result);

    return result;
}
