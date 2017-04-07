/*
 *  This file contains routines for computing image statistics.
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


/*
 *  TallyStats
 *
 *  Macro for computing the min, max, mean, and variance of pixel values.
 */

#define TallyStats(type)						\
    {									\
	type *pp = first_pixel;						\
									\
	tmin = tmax = *pp;						\
	for (i = 0; i < npixels; i++) {					\
	    tmean += pixel = *pp++;					\
	    tvar  += (pixel * pixel);					\
	    if (pixel < tmin)						\
		tmin = pixel;						\
	    else if (pixel > tmax)					\
		tmax = pixel;						\
	}								\
}


/*
 *  VImageStats
 *
 *  Compute the min, max, mean, and variance of an image's pixel values.
 */

VBoolean VImageStats (VImage src, VBand band, VDouble *pmin, VDouble *pmax,
		      VDouble *pmean, VDouble *pvar)
{
    int npixels, i;
    VPointer first_pixel;
    VDouble pixel, tmin, tmax, tmean, tvar;

    /* Prepare to iterate over the specified band(s) of source pixels: */
    if (! VSelectBand ("VImageStats", src, band, & npixels, & first_pixel))
	return FALSE;

    /* Initialize accumulators: */
    tmin = tmax = tmean = tvar = 0.0;

    /* Tally pixels: */
    switch (VPixelRepn (src)) {

    case VBitRepn:
	TallyStats (VBit);
	break;

    case VUByteRepn:
	TallyStats (VUByte);
	break;

    case VSByteRepn:
	TallyStats (VSByte);
	break;

    case VShortRepn:
	TallyStats (VShort);
	break;

    case VUShortRepn:
	TallyStats (VUShort);
	break;

    case VIntegerRepn:
	TallyStats (VInteger);
	break;

    case VUIntegerRepn:
	TallyStats (VUInteger);
	break;

    case VLongRepn:
	TallyStats (VLong);
	break;

    case VULongRepn:
	TallyStats (VULong);
	break;

    case VFloatRepn:
	TallyStats (VFloat);
	break;

    case VDoubleRepn:
	TallyStats (VDouble);
	break;

    default:
	break;
    }

    /* Compute and return the final results: */
    tmean /= npixels;
    tvar = tvar / npixels - (tmean * tmean);
    if (pmin)
	*pmin = tmin;
    if (pmax)
	*pmax = tmax;
    if (pmean)
	*pmean = tmean;
    if (pvar)
	*pvar = tvar;
    return TRUE;
}
