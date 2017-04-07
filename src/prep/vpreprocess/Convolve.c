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
VRcsId("$Id: Convolve.c 3181 2008-04-01 15:19:44Z karstenm $");

extern double rint(double);

/* Message reporting band index out of range: */
static VStringConst bandMessage =
  "VConvolveImage: Band %d referenced in image of %d band(s)";

/* Test whether representation type includes signed values: */
#define SignedRepn(repn) ((repn) != VBitRepn && (repn) != VUByteRepn)

/*
 *  ClampIndex
 *
 *  This macro implements behavior near the borders of the source image.
 *  Index is the band, row or column of the pixel being convolved.
 *  Limit is the number of bands, rows or columns in the source image.
 *  Label is a label to jump to to break off computation of the current
 *  destination pixel.
 */

#define ClampIndex(index, limit, label)				\
  {								\
    if (index < 0)						\
      switch (pad_method) {					\
      case VConvolvePadBorder:    index = 0; break;		\
      case VConvolvePadNone:      a = 0; goto label;		\
      case VConvolvePadWrap:      index += limit; break;	\
      default:            continue;				\
      }								\
    else if (index >= limit)					\
      switch (pad_method) {					\
      case VConvolvePadBorder:    index = limit - 1; break;	\
      case VConvolvePadNone:      a = 0; goto label;		\
      case VConvolvePadWrap:      index -= limit; break;	\
      default:            continue;				\
      }								\
  }


/*
 *  Convolve
 *
 *  This macro performs the actual convolution, iterating over pixels
 *  of the destination and mask images.
 */

#define Convolve(mask_type, src_type, dest_type, op, label)		\
  {									\
    mask_type a, *mask_pp;						\
    dest_type *dest_pp = VImageData (result);				\
									\
    for (dest_band = 0; dest_band < dest_nbands; dest_band++)           \
      for (dest_row = 0; dest_row < dest_nrows; dest_row++)		\
	for (dest_col = 0; dest_col < dest_ncols; dest_col++) {		\
	  a = 0;							\
	  for (i = 0; i < VImageNBands (mask); i++) {			\
	    src_band = dest_band + i + band_offset;			\
	    ClampIndex (src_band, VImageNBands (src), label);		\
	    for (j = 0; j < VImageNRows (mask); j++) {			\
	      src_row = dest_row + j + row_offset;			\
	      ClampIndex (src_row, VImageNRows (src), label);		\
	      mask_pp = mask_cvt->band_index[i][j];			\
	      for (k = 0; k < VImageNColumns (mask);			\
		   k++, mask_pp++) {					\
		src_col = dest_col + k + col_offset;			\
		ClampIndex (src_col, VImageNColumns (src), label);	\
		a += VPixel (src, src_band, src_row, src_col,		\
			     src_type) * *mask_pp;			\
	      }								\
	    }								\
	  }								\
	  op;								\
	label:    *dest_pp++ = a;					\
	}								\
  }

/* Clip a value to a single limit: */
#define Clip(value, test, limit)		\
  if (value test limit) {			\
    value = limit;                              \
    clipped = TRUE;                             \
  }

/* Clip a value to a pair of limits: */
#define Bracket(value, lower, upper)		\
  Clip (value, <, lower);			\
  Clip (value, >, upper);


/*
 *  VConvolveImage
 *
 *  Convolves one image with another.
 *  Returns a pointer to the resulting image if successful, zero otherwise.
 *  The band parameter may be VAllBands, in which case all bands of each
 *  band of the source image is convolved with the mask in turn, or it may
 *  be a particular band number, in which case only that single band is
 *  convolved. The pad_method parameter specifies what to do when convolving
 *  near the border of the source image: extend it with zeros, extend it with
 *  the value of the nearest border pixel, or wrap it around on itself.
 */

VImage VConvolveImage(VImage src, VImage dest, VBand band, VImage mask,
                      VConvolvePadMethod pad_method, int shift) {
  int dest_nbands, dest_nrows, dest_ncols;
  VRepnKind dest_repn, mask_repn;
  int band_offset, row_offset, col_offset;
  int src_band, src_row, src_col, dest_band, dest_row, dest_col, i, j, k;
  VImage result, mask_cvt;
  VDouble scale_factor = 1.0 / (1 << shift);
  VBoolean clipped = FALSE;
  /* Determine what dimensions the destination image should have:
     o if band specifies a specific band, the destination image has
     a single band; otherwise it can have several.
     o if the pad method is Trim, the destination image will have smaller
     dimensions than the source image (by an amount corresponding to the
     mask size); otherwise it will have the same dimensions. */
  if(band == VAllBands) {
    dest_nbands = VImageNBands(src);
    if(pad_method == VConvolvePadTrim)
      dest_nbands -= (VImageNBands(mask) - 1);
  } else if(band >= 0 && band < VImageNBands(src)) {
    if(pad_method == VConvolvePadTrim) {
      i = band - VImageNBands(mask) / 2;
      if(i < 0) {
	VWarning(bandMessage, i, VImageNBands(src));
	return NULL;
      }
      i = band + (VImageNBands(mask) - 1) / 2;
      if(i >= VImageNBands(src)) {
	VWarning(bandMessage, i, VImageNBands(src));
	return NULL;
      }
    }
    dest_nbands = 1;
  } else {
    VWarning(bandMessage, band, VImageNBands(src));
    return NULL;
  }

  dest_nrows = VImageNRows(src);
  dest_ncols = VImageNColumns(src);
  if(pad_method == VConvolvePadTrim) {
    dest_nrows -= (VImageNRows(mask) - 1);
    dest_ncols -= (VImageNColumns(mask) - 1);
    if(dest_nrows <= 0 || dest_ncols <= 0) {
      VWarning("VConvolveImage: "
	       "Image (%dx%d) is smaller than mask (%dx%d)",
	       VImageNRows(src), VImageNColumns(src),
	       VImageNRows(mask), VImageNColumns(mask));
      return NULL;
    }
  }

  /* Choose an appropriate destination representation:
     (Usually the destination image has the same representation as the
     source. If the source is Bit or UByte, however, and the mask is
     some signed representation, then the destination is SByte.) */
  dest_repn =
    (! SignedRepn(VPixelRepn(src)) && SignedRepn(VPixelRepn(mask))) ?
    VSByteRepn : VPixelRepn(src);


  /* Locate the destination. Since the operation cannot be carried out in
     place, we create a temporary work image to serve as the destination if
     dest == src or dest == mask: */
  result = VSelectDestImage("VConvolveImage", dest, dest_nbands,
			    dest_nrows, dest_ncols, dest_repn);
  if(! result)
    return NULL;
  if(src == result || mask == result) {
    result = VCreateImage(dest_nbands, dest_nrows, dest_ncols, dest_repn);
    if(! result)
      return NULL;
  }

  /* Convert the mask image to pixels of type VLong or VDouble, if necessary: */
  if(VIsFloatPtRepn(VPixelRepn(src)) ||
     VIsFloatPtRepn(VPixelRepn(mask)))
    mask_repn = VDoubleRepn;
  else
    mask_repn = VLongRepn;
  mask_cvt = (mask_repn == VPixelRepn(mask)) ? mask :
    VConvertImageCopy(mask, NULL, VAllBands, mask_repn);
  if(! mask_cvt) {
    if(result != dest)
      VDestroyImage(result);
    return NULL;
  }

  /* Determine the mapping from destination coordinates + mask coordinates to
     source coordinates: */
  if(band == VAllBands && pad_method == VConvolvePadTrim)
    band_offset = 0;
  else
    band_offset = (band == VAllBands ? 0 : band) - VImageNBands(mask) / 2;
  if(pad_method == VConvolvePadTrim)
    row_offset = col_offset = 0;
  else {
    row_offset = - (VImageNRows(mask) / 2);
    col_offset = - (VImageNColumns(mask) / 2);
  }


  /* Perform the convolution over all destination bands, rows, columns: */
  switch((int) VPixelRepn(mask_cvt)) {
  case VLongRepn:
    switch((int) VPixelRepn(src)) {
    case VBitRepn:
      if(dest_repn == VSByteRepn) {
	Convolve(VLong, VBit, VSByte,
		 a >>= shift; Bracket(a, -128, 127), L1);
      } else {
	Convolve(VLong, VBit, VBit, a >>= shift; Clip(a, > , 1), L2);
      }
      break;

    case VUByteRepn:
      if(dest_repn == VSByteRepn) {
	Convolve(VLong, VUByte, VSByte,
		 a >>= shift; Bracket(a, -128, 127), L3);
      } else {
	Convolve(VLong, VUByte, VUByte,
		 a >>= shift; Clip(a, > , 255), L4);
      }
      break;

    case VSByteRepn:
      Convolve(VLong, VSByte, VSByte,
	       a >>= shift; Bracket(a, -128, 127), L5);
      break;
    case VShortRepn:
      Convolve(VLong, VShort, VShort,
	       a >>= shift; Bracket(a, -32768, 32767), L6);
      break;
    case VLongRepn:
      Convolve(VLong, VLong, VLong, a >>= shift, L7);
    }
    break;


  case VDoubleRepn:
    switch((int) VPixelRepn(src)) {
    case VBitRepn:
      Convolve(VDouble, VBit, VSByte,
	       a = rint(a * scale_factor);
	       Bracket(a, -128.0, 127.0), L8);
      break;
    case VUByteRepn:
      Convolve(VDouble, VUByte, VSByte,
	       a = rint(a * scale_factor);
	       Bracket(a, -128.0, 127.0), L9);
      break;
    case VSByteRepn:
      Convolve(VDouble, VSByte, VSByte,
	       a = rint(a * scale_factor);
	       Bracket(a, -128.0, 127.0), L10);
      break;
    case VShortRepn:
      Convolve(VDouble, VShort, VShort,
	       a = rint(a * scale_factor);
	       Bracket(a, -32768.0, 32767.0), L11);
      break;
    case VLongRepn:
      Convolve(VDouble, VLong, VLong,
	       a = rint(a * scale_factor);
	       Bracket(a, -2147483648.0, 2147483647.0), L12);
      break;
    case VFloatRepn:
      Convolve(VDouble, VFloat, VFloat, a *= scale_factor, L13);
      break;
    case VDoubleRepn:
      Convolve(VDouble, VDouble, VDouble, a *= scale_factor, L14);
    }
  }
  if(clipped)
    VWarning("VConvolveImage: Destination pixel value(s) clipped");
  /* Successful completion: */
  if(mask_cvt != mask)
    VDestroyImage(mask_cvt);
  if(dest == src || dest == mask) {
    VCopyImagePixels(result, dest, VAllBands);
    VDestroyImage(result);
    return dest;
  } else {
    VCopyImageAttrs(src, result);
    return result;
  }
}
