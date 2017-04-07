/*
 *  $Id: ConvertO.c 3177 2008-04-01 14:47:24Z karstenm $
 *
 *  This file contains a routines for converting an image from one pixel
 *  representation to another by mapping the actual range of source pixel
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

/* File identification string: */
VRcsId ("$Id: ConvertO.c 3177 2008-04-01 14:47:24Z karstenm $");

/* Later in this file: */
static VImage ConstantImage (VImage, VImage, VBand, VRepnKind, VDouble);


/*
 *  VConvertImageOpt
 *
 *  Converts an image from one pixel representation to another using a
 *  mapping from the actual range of source pixel values to full range of
 *  destination pixel values, thus preserving as much as possible of the
 *  information in the source image. If preserve_zero is TRUE, the mapping
 *  will be such that 0 maps to 0.
 */

VImage VConvertImageOpt (VImage src, VImage dest, VBand band,
			 VRepnKind pixel_repn, int method)
{
    VDouble a_min, a_max, d_min, d_max, v, a, b, e;

    /* Determine range of pixel values in the source and destination images: */
    VImageStats (src, band, & a_min, & a_max, NULL, NULL);
    d_min = VIsFloatPtRepn (pixel_repn) ? -1.0 : VRepnMinValue (pixel_repn);
    d_max = VIsFloatPtRepn (pixel_repn) ? 1.0 : VRepnMaxValue (pixel_repn);

    /* Set v = max ( | a_min |, | a_max | ): */
    v = a_max > 0.0 ? a_max : - a_max;
    if (a_min > v)
	v = a_min;
    else if ( -a_min > v)
	v = -a_min;

    /* Choose an appropriate linear mapping from [a_min, a_max] to
       [d_min, d_max]: */
    switch (method) {

    case 1:

	/* A mapping that preserves sign and maps 0 -> 0: */
	if (a_min == 0.0 && a_max == 0.0)
	    return ConstantImage (src, dest, band, pixel_repn, 0.0);
	if (a_min < 0 && d_min >= 0.0)
	    VWarning ("VConvertImageOpt: Negative pixel value converted "
		      "to unsigned destination");
	return VConvertImageLinear (src, dest, band, pixel_repn,
				    d_max / v, 0.0);

    case 2:

	/* A mapping that preserves the sign of the actual pixel values mapped
	   but doesn't necessarily map 0 -> 0: */
	if (a_min == a_max)
	    return ConstantImage (src, dest, band, pixel_repn,
				  a_max < 0.0 ? d_min :
				  a_min > 0.0 ? d_max : 0.0);
	if (a_min < 0 && d_min >= 0.0)
	    VWarning ("VConvertImageOpt: Negative pixel value converted "
		      "to unsigned destination");
	if (a_max < 0.0) {
	    e = 1.0 / (a_max - a_min + 1.0);
	    a = (-e - d_min) / (a_max - a_min);
	    b = d_min - a_min * a;
	} else if (a_max == 0.0) {
	    a = d_min / a_min;
	    b = 0.0;
	} else if (a_min == 0.0) {
	    a = d_max / a_max;
	    b = 0.0;
	} else if (a_min > 0.0) {
	    e = 1.0 / (a_max - a_min + 1.0);
	    a = (d_max - e) / (a_max - a_min);
	    b = d_max - a_max * a;
	} else {
	    a = d_max / v;
	    b = 0.0;
	}
	return VConvertImageLinear (src, dest, band, pixel_repn,
				    (double) a, (double) b);

    case 3:

	/* A mapping from the range of actual source pixel values to the
	   range of possible destination pixel values, not necessarily
	   preserving either sign or zero: */
	if (a_min == a_max)
	    return ConstantImage (src, dest, band, pixel_repn,
				  a_max < 0.0 ? d_min :
				  a_min > 0.0 ? d_max : 0.0);
	a = (d_max - d_min) / (a_max - a_min);
	b = d_max - a_max * a;
	return VConvertImageLinear (src, dest, band, pixel_repn,
				    (double) a, (double) b);

    default:
	VError ("VConvertImageopt: Unknown conversion method %d", method);
    }
    return NULL;
}


/*
 *  ConstantImage -- return a constant image on behalf of VConvertImageOpt.
 */

static VImage ConstantImage (VImage src, VImage dest, VBand band,
			     VRepnKind pixel_repn, VDouble c)
{
    dest = VSelectDestImage ("VConvertImageOpt", dest,
			     band == VAllBands ? VImageNBands (src) : 1,
			     VImageNRows (src), VImageNColumns (src),
			     pixel_repn);
    if (dest)
	VFillImage (dest, VAllBands, c);
    return dest;
}
