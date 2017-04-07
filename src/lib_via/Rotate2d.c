/*! \file
  2D image rotation.


\par Reference:
"A Fast Algorithm for General Raster Rotation"
by Alan Paeth, Graphics Interface '86, pp. 77-81.


\par Author:
 Daniel Ko, UBC Laboratory for Computational Intelligence
*/


/*
 *  This file contains routines for rotating images.
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

/* From the standard C library: */
#include <math.h>


/* From the Vista library: */
#include "viaio/Vlib.h"
#include "viaio/os.h"
#include "viaio/VImage.h"
#include "viaio/mu.h"

#define PI 3.14159265

extern VImage RotateImage (VImage src, VImage dest, VBand band, double angle);



/*!
\fn VImage VRotateImage2d (VImage src, VImage dest, VBand band, double angle)
\param src   input image
\param dest  output image 
\param band  slice(s) to be processed.
\param angle totation angle.
*/
VImage VRotateImage2d (VImage src, VImage dest, VBand band, double angle)
{
#define Tolerance 0.01

  double angle_remained;
  VImage tmp1 = NULL, tmp2;

  /* Normalize angle to the range [-Pi, +Pi]: */
  angle = fmod (angle, 2.0 * PI);
  if (angle > PI)
    angle -= 2.0 * PI;

  /* 
   * Two stage rotation:
   * 1. Use flipping to rotate the image by a multiple of 90 degrees.
   * 2. Use RotateImage() to rotate the remaining angle.
   */
  if (angle >= -0.25 * PI && angle <= 0.25 * PI) {

    /* Do nothing: */
    angle_remained = angle;
    tmp2 = src;

  } else if (angle >= 0.25 * PI && angle <= 0.75 * PI) {

    /* Rotate by 90 degrees: */
    tmp1 = VTransposeImage (src, NULL, band);
    tmp2 = VFlipImage (tmp1, NULL, VAllBands, TRUE);
    angle_remained = angle - 0.5 * PI;

  } else if (angle >= 0.75 * PI || angle <= -0.75 * PI) {

    /* Rotate by 180 degrees: */
    tmp1 = VFlipImage (src, NULL, band, TRUE);
    tmp2 = VFlipImage (tmp1, NULL, VAllBands, FALSE);
    angle_remained = angle > 0.0 ? angle - PI : angle + PI;

  } else {

    /* Rotate by -90 degress: */
    tmp1 = VTransposeImage (src, NULL, band);
    tmp2 = VFlipImage (tmp1, NULL, VAllBands, FALSE);
    angle_remained = angle + 0.5 * PI;
  }

  if (VMax (angle_remained, -angle_remained) > Tolerance) {

    /* Go on to stage 2: */
    dest = RotateImage (tmp2, dest, tmp2 == src ? band : VAllBands, 
			angle_remained);

  } else {

    /* Stage 1 only is good enough: */
    dest = VCopyImage (tmp2, dest, tmp2 == src ? band : VAllBands);
  }
    
  /* Clean up: */
  if (tmp2 != src) {
    VDestroyImage (tmp1);
    VDestroyImage (tmp2);
  }
    
  return dest;
    
#undef Tolerance
}


/*
 *  RotateImage
 *
 *  Rotate an image by an angle in the interval [-Pi/2, +Pi/2].
 *
 *  The rotation algorithm is adapted from
 *  "A Fast Algorithm for General Raster Rotation"
 *  by Alan Paeth, Graphics Interface '86, pp. 77-81.
 */

VImage RotateImage (VImage src, VImage dest, VBand band, double angle)
{
  double shearx, sheary;
  int margin_nrows, margin_ncols;
  int src_nbands;
  VImage sx_image, sysx_image, sxsysx_image, cropped_image;
    
  /* Ensure that range of "angle" is valid: */
  if (angle <= -0.5 * PI || angle >= 0.5 * PI) {
    VWarning ("VRotateImage: Internal error in RotateImage");
    return NULL;
  }

  /* Determine "shearx" and "sheary": */
  shearx = - tan (angle / 2.0);
  sheary = sin (angle);
    
  /* Read properties of "src": */
  src_nbands = VImageNBands (src);

    
  /* Check to ensure that "band" exists: */
  if (band != VAllBands && (band < 0 || band >= src_nbands)) {
    VWarning ("VRotateImage: Band %d referenced in image of %d bands",
	      band, VImageNBands (src));
    return NULL;
  }

  /* First shear in x direction: */
  sx_image = VShearImageX (src, NULL, band, shearx);
  if (sx_image == NULL)
    return NULL;

  /* Then shear in y direction: */
  sysx_image = VShearImageY (sx_image, NULL, VAllBands, sheary);
  if (sysx_image == NULL)
    return NULL;

  /* Finally, shear in x direction again: */
  sxsysx_image = VShearImageX (sysx_image, NULL, VAllBands, shearx);
  if (sxsysx_image == NULL)
    return NULL;

  /* Calculate margin (unwanted blank space) sizes: */
  margin_nrows = VImageNRows (src) * fabs (shearx) * fabs (sheary);
  margin_ncols =
    (VImageNRows (sysx_image) - VImageNRows (src)) * fabs (shearx);

  /* Remove margins: */
  cropped_image = VCropImage (sxsysx_image, dest, VAllBands,
			      margin_nrows, margin_ncols,
			      VImageNRows (sxsysx_image) -
			      2 * margin_nrows,
			      VImageNColumns (sxsysx_image) -
			      2 * margin_ncols);

  /* Clean up: */
  VDestroyImage (sx_image);
  VDestroyImage (sysx_image);
  VDestroyImage (sxsysx_image);

  VCopyImageAttrs (src, cropped_image);
    
  return cropped_image;
}
