/*
 *  $Id: Shear.c 3177 2008-04-01 14:47:24Z karstenm $
 *
 *  This file contains routines for shearing images.
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

/* From the Vista library: */
#include "viaio/Vlib.h"
#include "viaio/os.h"
#include "viaio/VImage.h"
#include "viaio/mu.h"

/* From the standard C libaray: */
#include <math.h>

/* File identification string: */
VRcsId ("$Id: Shear.c 3177 2008-04-01 14:47:24Z karstenm $");

/* Macro(s): */
#define XYToRow(X,Y,DX,DY) ((DY)-(Y))   /* From xy-coord. to row-coord. */
#define XYToCol(X,Y,DX,DY) ((DX)+(X))   /* From xy-coord. to column-coord. */


/*
 * VShearImageX
 *
 * Shear an image in the x direction (i.e. the direction in which
 * column increases).
 *
 * The parameter "shear" is the lower left element of a shear matrix:
 *
 *   1        0
 *   shear    1
 *
 * The shearing algorithm is adapted from
 * "A Fast Algorithm for General Raster Rotation"
 * by Alan Paeth, Graphics Interface '86, pp. 77-81.
 */

VImage VShearImageX (VImage src, VImage dest, VBand band, double shear)
{
  int src_nrows, src_ncols, src_nbands;
  int dest_nrows, dest_ncols, dest_nbands;
  VRepnKind src_repn, dest_repn;
  double skew, skewi, skewf, left, oleft;
  int x, y, b, src_curband, dest_curband;
  int extra_ncols, sdx, sdy, ddx, ddy;

  /* Read properties of "src": */
  src_nrows  = VImageNRows (src);
  src_ncols  = VImageNColumns (src);
  src_nbands = VImageNBands (src);
  src_repn   = VPixelRepn (src);

  /* Check to ensure that "band" exists: */
  if (band != VAllBands && (band < 0 || band >= src_nbands)) {
    VWarning ("VShearImageY: Band %d referenced in image of %d bands",
	      band, VImageNBands (src));
    return NULL;
  }
    
  /* Determine properties of "dest": */
  dest_nrows  = src_nrows;
  extra_ncols = (int) floor (shear * (src_nrows - 0.5));
  dest_ncols  = src_ncols + VMax (extra_ncols, -extra_ncols) + 1;
  dest_nbands = (band == VAllBands) ? src_nbands : 1;
  dest_repn   = src_repn;

  /* Calculate displacements in converting from xy-coord. to
     row/column-coord. : */
  sdx = 0;
  sdy = src_nrows - 1;
  ddx = (extra_ncols < 0) ? -extra_ncols : 0;
  ddy = dest_nrows - 1;
    
  /* Create dest image */
  dest = VSelectDestImage ("VShearImageX", dest,
			   dest_nbands, dest_nrows, dest_ncols, dest_repn);
    
  if (dest == NULL)
    return NULL;

  /* Set all dest pixels to zero: */
  VFillImage (dest, VAllBands, 0.0);
    
  /* Treat lower left-hand corner of image as origin, with
     y up and x to the right. */
    
  /*
   * Shear a row of an image of a particular type: 
   */
#define ShearRow(type)                                         \
    {                                                              \
	type pixel;                                                \
	for (x = 0; x < src_ncols; x++) {                          \
	    pixel = VPixel (src,                                   \
			    src_curband,                           \
			    XYToRow (x, y, sdx, sdy),              \
			    XYToCol (x, y, sdx, sdy),              \
			    type);                                 \
	    left = pixel * skewf;                                  \
	    pixel = pixel - left + oleft;                          \
	    VPixel (dest,                                          \
		    dest_curband,                                  \
		    XYToRow ((int) (x + skewi), y, ddx, ddy),      \
		    XYToCol ((int) (x + skewi), y, ddx, ddy),      \
		    type)                                          \
		= pixel;                                           \
	    oleft = left;                                          \
	}                                                          \
	VPixel (dest,                                              \
		dest_curband,                                      \
		XYToRow ((int) (src_ncols + skewi), y, ddx, ddy),  \
		XYToCol ((int) (src_ncols + skewi), y, ddx, ddy),  \
		type)                                              \
	    = (type) oleft;                                        \
    }

  /* For each band in the dest image do: */
  for (b = 0; b < dest_nbands; b++) {

    src_curband = (band == VAllBands) ? b : band;
    dest_curband = (band == VAllBands) ? b : 0;
	
    /* For each row in the source image do: */
    for (y = 0; y < src_nrows; y++) {
	    
      skew = shear * (y + 0.5);
      skewi = floor (skew);
      skewf = skew - skewi;

      oleft = 0.0;

      /* Shear a row according to pixel representation: */
      switch (src_repn) {
      case VBitRepn:      ShearRow (VBit);    break;
      case VUByteRepn:	ShearRow (VUByte);  break;
      case VSByteRepn:	ShearRow (VSByte);  break;
      case VShortRepn:	ShearRow (VShort);  break;
      case VLongRepn:	ShearRow (VLong);   break;
      case VFloatRepn:	ShearRow (VFloat);  break;
      case VDoubleRepn:	ShearRow (VDouble); break;
      default: break;
      }
    }
  }

  VCopyImageAttrs (src, dest);
  return dest;

#undef ShearRow
}


/*
 * VShearImageY
 *
 * Shear an image in the y direction (i.e. the direction in which
 * row decreases).
 *
 * The parameter "shear" is the upper right element of a shear matrix:
 *
 *   1    shear
 *   o    1
 *
 * The shearing algorithm is adapted from
 * "A Fast Algorithm for General Raster Rotation"
 * by Alan Paeth, Graphics Interface '86, pp. 77-81.
 */

VImage VShearImageY (VImage src, VImage dest, VBand band, double shear)
{
  int src_nrows, src_ncols, src_nbands;
  int dest_nrows, dest_ncols, dest_nbands;
  VRepnKind src_repn, dest_repn;
  double skew, skewi, skewf, left, oleft;
  int x, y, b, src_curband, dest_curband;
  int extra_nrows, sdx, sdy, ddx, ddy;

  /* Read properties of "src": */
  src_nrows  = VImageNRows (src);
  src_ncols  = VImageNColumns (src);
  src_nbands = VImageNBands (src);
  src_repn   = VPixelRepn (src);

  /* Check to ensure that "band" exists: */
  if (band != VAllBands && (band < 0 || band >= src_nbands)) {
    VWarning ("VShearImageY: Band %d referenced in image of %d bands",
	      band, VImageNBands (src));
    return NULL;
  }
    
  /* Determine properties of "dest": */
  extra_nrows = (int) floor (shear * (src_ncols - 0.5));
  dest_nrows  = src_nrows + VMax (extra_nrows, -extra_nrows) + 1;
  dest_ncols  = src_ncols;
  dest_nbands = (band == VAllBands) ? src_nbands : 1;
  dest_repn   = src_repn;

  /* Calculate displacements in converting from xy-coord. to
     row/column-coord. : */
  sdx = 0;
  sdy = src_nrows - 1;
  ddx = 0;
  ddy = dest_nrows - 1 - ((extra_nrows < 0) ? -extra_nrows : 0);
    
  /* Create dest image */
  dest = VSelectDestImage ("VShearImageY", dest, 
			   dest_nbands, dest_nrows, dest_ncols, dest_repn);
  if (dest == NULL)
    return NULL;

  /* Set all dest pixels to zero: */
  VFillImage (dest, VAllBands, 0.0);
    
  /* Treat lower left-hand corner of image as origin, with
     y up and x to the right. */


  /*
   * Shear a column of an image of a particular type: 
   */
#define ShearColumn(type)                                           \
    {                                                                   \
	type pixel;                                                     \
	for (y = 0; y < src_nrows; y++) {                               \
	    pixel = VPixel (src,                                        \
			    src_curband,                                \
			    XYToRow (x, y, sdx, sdy),                   \
			    XYToCol (x, y, sdx, sdy),                   \
			    type);                                      \
	    left = pixel * skewf;                                       \
	    pixel = pixel - left + oleft;                               \
	    VPixel (dest,                                               \
		    dest_curband,                                       \
		    XYToRow (x, (int) (y + skewi),                      \
			     ddx, ddy),                                 \
		    XYToCol (x, (int) (y + skewi),                      \
			     ddx, ddy),                                 \
		    type)                                               \
		= pixel;                                                \
	    oleft = left;                                               \
	}                                                               \
	VPixel (dest,                                                   \
		dest_curband,                                           \
		XYToRow (x, (int) (src_nrows + skewi), ddx, ddy),       \
		XYToCol (x, (int) (src_nrows + skewi), ddx, ddy),       \
		type)                                                   \
	    = (type) oleft;                                             \
    }

  /* For each band in the dest image do: */
  for (b = 0; b < dest_nbands; b++) {

    src_curband = (band == VAllBands) ? b : band;
    dest_curband = (band == VAllBands) ? b : 0;
	
    /* For each column in the source image do: */
    for (x = 0; x < src_ncols; x++) {
	    
      skew = shear * (x + 0.5);
      skewi = floor (skew);
      skewf = skew - skewi;
	    
      oleft = 0.0;

      /* Shear a column according to pixel representation: */
      switch (src_repn) {
      case VBitRepn:      ShearColumn (VBit);    break;
      case VUByteRepn:	ShearColumn (VUByte);  break;
      case VSByteRepn:	ShearColumn (VSByte);  break;
      case VShortRepn:	ShearColumn (VShort);  break;
      case VLongRepn:	ShearColumn (VLong);   break;
      case VFloatRepn:	ShearColumn (VFloat);  break;
      case VDoubleRepn:	ShearColumn (VDouble); break;
      default: break;
      }
	    
    }

  }

  VCopyImageAttrs (src, dest);
  return dest;

#undef ShearColumn
}
