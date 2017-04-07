/*
 * $Id: Crop.c 3177 2008-04-01 14:47:24Z karstenm $
 * 
 * This file contains routines for inverting and negating images.
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
 *  Authors: Ralph Horstman, Art Pope,
 *	     UBC Laboratory for Computational Intelligence
 */

/* From the Vista library: */
#include "viaio/Vlib.h"
#include "viaio/os.h"
#include "viaio/VImage.h"
#include "viaio/mu.h"

/* File identification string: */
VRcsId ("$Id: Crop.c 3177 2008-04-01 14:47:24Z karstenm $");


/*
 *  VCropImage
 *
 *  Copy a rectangular region from an image.
 */

VImage VCropImage (VImage src, VImage dest, VBand band,
		   int top, int left, int height, int width)
{
    int nbands, src_band, dest_band, src_row, dest_row;
    int first_src_row, last_src_row, first_src_column;
    int first_dest_row, first_dest_column;
    size_t copy_len;
    VImage result;
    
    /* Check arguments: */
    if (height < 0 || width < 0) {
	VWarning("VCropImage: Cropped dimensions are negative");
	return NULL;
    }
    if (band == VAllBands) {
	src_band = 0;
	nbands = VImageNBands (src);
    } else if (band >= 0 && band < VImageNBands (src)) {
	src_band = band;
	nbands = 1;
    } else {
	VWarning ("VCropImage: Band %d referenced in image of %d band(s)",
		  band, VImageNBands (src));
	return NULL;
    }

    /* Create or check the destination image: */
    result = VSelectDestImage ("VCropImage", dest, nbands, 
			       height, width, VPixelRepn(src));
    if (! result)
	return NULL;
    if (src == result)
	result = VCreateImage (nbands, height, width, VPixelRepn (src));
    
    /* Clear the destination image: */
    VFillImage (result, VAllBands, 0.0);

    /* Copy over the cropped portion of each band of the source image: */
    first_src_row = VMax (top, 0);
    last_src_row = VMin (top + height, VImageNRows (src)) - 1;
    first_src_column = VMax (left, 0);
    first_dest_row = VMax (-top, 0);
    first_dest_column = VMax (-left, 0);
    copy_len = (VMin (left + width, VImageNColumns (src)) - VMax (left, 0)) *
	VPixelSize (src);
    for (dest_band = 0; dest_band < nbands; dest_band++, src_band++)
	for (src_row = first_src_row, dest_row = first_dest_row;
	     src_row <= last_src_row; src_row++, dest_row++)
	    memcpy (VPixelPtr (result, dest_band, dest_row, first_dest_column),
		    VPixelPtr (src, src_band, src_row, first_src_column),
		    copy_len);

    if (src == dest) {
	VCopyImagePixels (result, dest, VAllBands);
	VDestroyImage (result);
	return dest;
    } else {
	VCopyImageAttrs (src, result);
	return result;
    }
}
