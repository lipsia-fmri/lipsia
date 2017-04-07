/*
 *  $Id: GaussConv.c 3181 2008-04-01 15:19:44Z karstenm $
 *
 *  This file contains routines for performing Gaussian convolution.
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
 *  Author: Daniel Ko, Art Pope, UBC Laboratory for Computational Intelligence
 */

/* From the Vista library: */
#include <math.h>
#include "viaio/Vlib.h"
#include "viaio/mu.h"
#include "viaio/os.h"
#include "viaio/VImage.h"

/* File identification string: */
VRcsId("$Id: GaussConv.c 3181 2008-04-01 15:19:44Z karstenm $");


/*
 *  VGaussianKernel
 *
 *  Generate a 1D Gaussian kernel.
 */
VBoolean VGaussianKernel(int ncoeffs, VDouble coeffs[], double s) 
{
  double sum, x, b = -1.0 / (2.0 * s * s);
  VDouble *p1, *p2, *pend = coeffs + ncoeffs;

  /* Setup depends on whether the number of coefficients asked for is odd or even: */
  if(ncoeffs & 1) {
    p1 = & coeffs[ncoeffs / 2];
    *p1-- = 1.0;
    p2 = p1 + 2;
    x = 1.0;
    sum = 0.5;
  } else {
    p1 = & coeffs[ncoeffs / 2];
    p2 = p1 + 1;
    x = 0.5;
    sum = 0.0;
  }

  /* Fill the vector with coefficients, then make them sum to 1.0: */
  while(p2 < pend) {
    sum += *p1-- = *p2++ = exp(x * x * b);
    x += 1.0;
  }
  sum += sum;
  for(p1 = coeffs; p1 < pend; p1++)
    *p1 /= sum;
  return TRUE;
}


/*
 *  VConvolveImageSep
 *
 *  Convolve an image with a separable filter.
 *  The filter is specified as separate band-wise, row-wise, and column-
 *  filters. Each is applied in sucession in accordance with its respective
 *  pad method and shift. If a specific band is specified, a one-band
 *  image is produced by centering the kernel on that band of src.
 */
VImage VConvolveImageSep(VImage src, VImage dest, VBand band, VImage masks[3],
                         VConvolvePadMethod pad_methods[3], int shifts[3]) 
{
  int last_filter, i;
  VImage result, prev_image;

  /* Identify the last in the series of filters: */
  for(last_filter = 2; last_filter >= 0; last_filter--)
    if(masks[last_filter])
      break;
  if(last_filter == -1) {
    VWarning("VConvolveSep: No mask specified");
    return NULL;
  }
  result = src;

  for(i = 0; i < 3; i++) {
    prev_image = result;
    if(! masks[i]) continue;
    result = VConvolveImage(prev_image, (i == last_filter) ? dest : NULL,
			    band, masks[i], pad_methods[i], shifts[i]);
    if(prev_image != src)
      VDestroyImage(prev_image);
    band = VAllBands;
  }
  return result;
}


/*
 *  VGaussianConvolveImage
 *
 *  Convolve an image with a Gaussian filter, or a filter that is the
 *  partial derivate of a Gaussian in the X or Y direction.
 *  The 2D filter is implemented as two 1D filters, each with the specified
 *  standard deviation and length.
 */
VImage VGaussianConv(VImage src, VImage dest, VBand band,double sigma, int filter_size) 
{
  VImage src_float, result, tmp, filters[3];
  VRepnKind dest_pixel_repn, conv_pixel_repn;
  static VConvolvePadMethod pad_methods[3] =
    { VConvolvePadBorder, VConvolvePadBorder, VConvolvePadBorder };
  static int shifts[3] = { 0, 0, 0 };

  /* Ensure that the standard deviation and filter size are positive: */
  if(sigma <= 0.0) {
    VWarning("VGaussianConvolveImage: Sigma (%g) is not positive", sigma);
    return NULL;
  }
  if(filter_size <= 0) {
    VWarning("VGaussianConvolveImage: Filter size (%d) is not positive",
	     filter_size);
    return NULL;
  }

  /* Determine the appropriate destination pixel type. Convolving a bit
     or ubyte image with a derivative of a Gaussian produces an sbyte image;
     otherwise the destination has the same type as src. */
  dest_pixel_repn = VPixelRepn(src);

  /* Convert the source image to Float or Double for the sake of efficiency: */
  if(! VIsFloatPtRepn(VPixelRepn(src))) {
    if(!(src_float = VConvertImageRange(src, NULL, band, VFloatRepn)))
      return NULL;
    band = VAllBands;
  } else
    src_float = src;
  conv_pixel_repn = VPixelRepn(src_float);

  /* Make separate row and column filters: */
  filters[0] = NULL;
  filters[1] = VCreateImage(1, filter_size, 1, VDoubleRepn);
  filters[2] = VCreateImage(1, 1, filter_size, VDoubleRepn);
  VGaussianKernel(filter_size, VImageData(filters[1]), sigma);
  VGaussianKernel(filter_size, VImageData(filters[2]), sigma);

  /* Convolve src with the two filters: */
  result =
    VConvolveImageSep(src_float, (conv_pixel_repn == dest_pixel_repn) ?
		      dest : NULL, band, filters, pad_methods, shifts);

  /* Convert the result to the apprporate destination type, if necessary: */
  if(result && conv_pixel_repn != dest_pixel_repn) {
    tmp = result;
    result = VConvertImageRange(tmp, dest, VAllBands, dest_pixel_repn);
    VDestroyImage(tmp);
  }

  if(src_float != src) VDestroyImage(src_float);
  VDestroyImage(filters[1]);
  VDestroyImage(filters[2]);
  return result;
}
