/*
 *  $Id: Image.c 3177 2008-04-01 14:47:24Z karstenm $
 *
 *  This file provides basic support for images (the VImage class).
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
#include "viaio/file.h"
#include "viaio/os.h"
#include "viaio/VImage.h"

/* From the standard C library: */
#include <math.h>


/*
 *  Dictionaries for VImage attributes and their values.
 */

/* Keywords for representing band interpretation values: */
VDictEntry VBandInterpDict[] = {
  { "complex",	VBandInterpComplex },
  { "gradient",	VBandInterpGradient },
  { "intensity",	VBandInterpIntensity },
  { "orientation",	VBandInterpOrientation },
  { "rgb",		VBandInterpRGB },
  { "stereo_pair",	VBandInterpStereoPair },
  { NULL }
};


/*
 *  VCreateImage
 *
 *  Allocates memory for a new image with specified properties.
 *  Returns a pointer to the image if successful, zero otherwise.
 */

VImage VCreateImage (int nbands, int nrows, int ncolumns, VRepnKind pixel_repn)
{
  size_t row_size = ncolumns * VRepnSize (pixel_repn);
  size_t data_size = nbands * nrows * row_size;
  size_t row_index_size = nbands * nrows * sizeof (char *);
  size_t band_index_size = nbands * sizeof (char **);
  size_t pixel_size;
  char *p;
  VImage image;
  int band, row;

#define AlignUp(v, b) ((((v) + (b) - 1) / (b)) * (b))

  /* Check parameters: */
  if (nbands < 1) {
    VWarning ("VCreateImage: Invalid number of bands: %d", (int) nbands);
    return NULL;
  }
  if (nrows < 1) {
    VWarning ("VCreateImage: Invalid number of rows: %d", (int) nrows);
    return NULL;
  }
  if (ncolumns < 1) {
    VWarning ("VCreateImage: Invalid number of columns: %d",
	      (int) ncolumns);
    return NULL;
  }
  if (pixel_repn != VBitRepn && pixel_repn != VUByteRepn &&
      pixel_repn != VSByteRepn && pixel_repn != VShortRepn &&
      pixel_repn != VLongRepn && pixel_repn != VFloatRepn &&
      pixel_repn != VIntegerRepn && pixel_repn != VUIntegerRepn &&
      pixel_repn != VUShortRepn && pixel_repn != VULongRepn &&
      pixel_repn != VDoubleRepn) {
    VWarning ("VCreateImage: Invalid pixel representation: %d",
	      (int) pixel_repn);
    return NULL;
  }

  /* Allocate memory for the VImage, its indices, and pixel values, while
     padding enough to ensure pixel values are appropriately aligned: */
  pixel_size = VRepnSize (pixel_repn);
  p = VMalloc (AlignUp (sizeof (VImageRec) + row_index_size +
			band_index_size, pixel_size) + data_size);

  /* Initialize the VImage: */
  image = (VImage) p;
  image->nbands = nbands;
  image->nrows = nrows;
  image->ncolumns = ncolumns;
  image->flags = VImageSingleAlloc;
  image->pixel_repn = pixel_repn;
  image->attributes = VCreateAttrList ();
  image->band_index = (VPointer **) (p += sizeof (VImageRec));
  image->row_index = (VPointer *) (p += band_index_size);
  image->data = (VPointer) AlignUp ((long) p + row_index_size, pixel_size);
  image->nframes = nbands;
  image->nviewpoints = image->ncolors = image->ncomponents = 1;

  /* Initialize the indices: */
  for (band = 0; band < nbands; band++)
    image->band_index[band] = image->row_index + band * nrows;
  for (row = 0, p = image->data; row < nbands * nrows; row++, p += row_size)
    image->row_index[row] = p;

  return image;

#undef AlignUp
}


/*
 *  VCreateImageLike
 *
 *  Create an image with the same properties as an existing one.
 */

VImage VCreateImageLike (VImage src)
{
  return VCopyImageAttrs (src, NULL);
}	


/*
 *  VDestroyImage
 *
 *  Frees memory occupied by an image.
 */

void VDestroyImage (VImage image)
{
  if (! image)
    return;
  if (! (image->flags & VImageSingleAlloc)) {
    VFree (image->data);
    VFree ((VPointer) image->row_index);
    VFree ((VPointer) image->band_index);
  }
  VDestroyAttrList (VImageAttrList (image));
  VFree ((VPointer) image);
}


/*
 *  VGetPixel
 *
 *  Fetch a pixel value, regardless of type, and return it as a Double.
 */

VDouble VGetPixel (VImage image, int band, int row, int column)
{
  VPointer p = VPixelPtr (image, band, row, column);

  switch (VPixelRepn (image)) {

  case VBitRepn:
    return (VDouble) * (VBit *) p;

  case VUByteRepn:
    return (VDouble) * (VUByte *) p;

  case VSByteRepn:
    return (VDouble) * (VSByte *) p;

  case VShortRepn:
    return (VDouble) * (VShort *) p;

  case VUShortRepn:
    return (VDouble) * (VUShort *) p;

  case VIntegerRepn:
    return (VDouble) * (VInteger *) p;

  case VUIntegerRepn:
    return (VDouble) * (VUInteger *) p;

  case VLongRepn:
    return (VDouble) * (VLong *) p;

  case VULongRepn:
    return (VDouble) * (VULong *) p;


  case VFloatRepn:
    return (VDouble) * (VFloat *) p;

  case VDoubleRepn:
    return (VDouble) * (VDouble *) p;

  default:
    VError ("VGetPixel: %s images not supported", VPixelRepnName (image));
  }
  return 0.0;		/* to make lint happy */
}


/*
 *  VSetPixel
 *
 *  Set a pixel, regardless of type, and to a value passed as a Double.
 */

void VSetPixel (VImage image, int band, int row, int column,
		VDoublePromoted value)
{
  VPointer p = VPixelPtr (image, band, row, column);

  switch (VPixelRepn (image)) {

  case VBitRepn:
    * (VBit *) p = value;
  break;

  case VUByteRepn:
    * (VUByte *) p = value;
  break;

  case VSByteRepn:
    * (VSByte *) p = value;
  break;

  case VShortRepn:
    * (VShort *) p = value;
  break;  

  case VUShortRepn:
    * (VUShort *) p = value;
  break;

  case VIntegerRepn:
    * (VInteger *) p = value;
  break;

  case VUIntegerRepn:
    * (VUInteger *) p = value;
  break;

  case VLongRepn:
    * (VLong *) p = value;
  break;

  case VFloatRepn:
    * (VFloat *) p = value;
  break;

  case VDoubleRepn:
    * (VDouble *) p = value;
  break;

  default:
    VError ("VSetPixel: %s images not supported", VPixelRepnName (image));
  }
}


/*
 *  VCopyImage
 *
 *  Copy the pixels and attributes of one image to another.
 *  Returns a pointer to the destination image if successful, zero otherwise.
 *  The band parameter may be VAllBands, in which case all bands of pixel
 *  values are copied, or a particular band number, in which case only a
 *  single band is copied to a 1-band destination image.
 */

VImage VCopyImage (VImage src, VImage dest, VBand band)
{
  VImage result;

  if (src == dest &&
      (band == VAllBands || (band == 0 && VImageNBands (src) == 1)))
    return src;

  if ((result = VCopyImagePixels (src, dest, band)) != 0)
    VCopyImageAttrs (src, result);
  return result;
}


/*
 *  VCopyImageAttrs
 *
 *  Give a destination image the same attributes as a source image.
 *  However if the destination image doesn't have the same number of bands
 *  as the source image, any band interpretation attributes are deleted.
 */

VImage VCopyImageAttrs (VImage src, VImage dest)
{
  VAttrList list;

  if (src == dest)
    return dest;
  if (! dest) {
    dest = VCreateImage (VImageNBands (src), VImageNRows (src),
			 VImageNColumns (src), VPixelRepn (src));
    if (! dest)
      return NULL;
  }

  /* Clone the source image's attribute list if it isn't empty: */
  if (! VAttrListEmpty (VImageAttrList (src))) {
    list = VImageAttrList (dest);
    VImageAttrList (dest) = VCopyAttrList (VImageAttrList (src));
  } else if (! VAttrListEmpty (VImageAttrList (dest))) {
    list = VImageAttrList (dest);
    VImageAttrList (dest) = VCreateAttrList ();
  } else list = NULL;
  if (list)
    VDestroyAttrList (list);

  /* Preserve band interpretation attributes only if the source and
     destination images have the same number of bands: */
  if (VImageNBands (src) > 1 && VImageNBands (dest) == VImageNBands (src)) {
    VImageNFrames (dest) = VImageNFrames (src);
    VImageNViewpoints (dest) = VImageNViewpoints (src);
    VImageNColors (dest) = VImageNColors (src);
    VImageNComponents (dest) = VImageNComponents (src);
  } else {
    VExtractAttr (VImageAttrList (dest), VFrameInterpAttr, NULL,
		  VBitRepn, NULL, FALSE);
    VExtractAttr (VImageAttrList (dest), VViewpointInterpAttr, NULL,
		  VBitRepn, NULL, FALSE);
    VExtractAttr (VImageAttrList (dest), VColorInterpAttr, NULL,
		  VBitRepn, NULL, FALSE);
    VExtractAttr (VImageAttrList (dest), VComponentInterpAttr, NULL,
		  VBitRepn, NULL, FALSE);
    VImageNComponents (dest) = VImageNColors (dest) = 
      VImageNViewpoints (dest) = 1;
    VImageNFrames (dest) = VImageNBands (dest);
  }
  return dest;
}


/*
 *  VCopyImagePixels
 *
 *  Copy the pixels of one image to another.
 *  Returns a pointer to the destination image if successful, zero otherwise.
 *  The band parameter may be VAllBands, in which case all bands of pixel
 *  values are copied, or a particular band number, in which case only a
 *  single band is copied to a 1-band destination image.
 */

VImage VCopyImagePixels (VImage src, VImage dest, VBand band)
{
  int npixels;
  VPointer src_pixels;
  VImage result;

  /* Locate the source and destination of the copy: */
  if (! VSelectBand ("VCopyImagePixels", src, band, & npixels, & src_pixels))
    return NULL;
  result = VSelectDestImage ("VCopyImagePixels", dest,
			     band == VAllBands ? VImageNBands (src) : 1,
			     VImageNRows (src), VImageNColumns (src),
			     VPixelRepn (src));
  if (! result)
    return NULL;

  /* Copy pixel values from src to dest: */
  memcpy (VImageData (result), src_pixels, npixels * VPixelSize (src));

  return result;
}


/*
 *  VCopyBand
 *
 *  Copy a band of pixel data from one image to another.
 *  Band src_band of image src is copied to band dest_band of image dest.
 *  The destination image must exist, having the same pixel representation
 *  and size as the source image. Either src_band or dst_band may be
 *  VAllBands, provided they both describe the same number of bands.
 */

VBoolean VCopyBand (VImage src, VBand src_band, VImage dest, VBand dest_band)
{
  int nbands, src_npixels, dest_npixels;
  VPointer src_pixels, dest_pixels;

  /* The destination image must exist: */
  if (! dest) {
    VWarning ("VCopyBand: No destination specified");
    return FALSE;
  }

  /* VAllBands not accepted for destination band: */
  if (dest_band < 0 || dest_band >= VImageNBands (dest)) {
    VWarning ("VCopyBand: Band %d referenced in image of %d bands",
	      (int) dest_band, (int) VImageNBands (dest));
    return FALSE;
  }

  /* Ensure that the destination image has the appropriate dimensions
     and pixel representation: */
  nbands = dest_band;
  if (src_band == VAllBands)
    nbands += VImageNBands (src) - 1;
  if (nbands < VImageNBands (dest))
    nbands = VImageNBands (dest);
  if (! VSelectDestImage ("VCopyBand", dest, nbands, VImageNRows (src),
			  VImageNColumns (src), VPixelRepn (src)))
    return FALSE;

  /* Locate the specified source and destination bands: */
  if (! VSelectBand ("VCopyBand", src, src_band,
		     & src_npixels, & src_pixels))
    return FALSE;
  if (! VSelectBand ("VCopyBand", dest, dest_band,
		     & dest_npixels, & dest_pixels))
    return FALSE;

  /* Copy from the source band to the destination band: */
  memcpy (dest_pixels, src_pixels, src_npixels * VPixelSize (src));

  return TRUE;
}


/*
 *  VCombineBands
 *
 *  Copy a series of bands from various images into a destination image.
 *  nbands is the number of bands to be copied, and the number of elements
 *  in src_images (a list of images from which to obtain the bands) and
 *  src_bands (a list of their respective band numbers). The bands are
 *  copied into the destination image in sequence.
 */

VImage VCombineBands (int nels, VImage src_images[], VBand src_bands[],
		      VImage dest)
{
  int n, i;
  VImage result, src = src_images[0];

  /* Count the number of bands needed in the destination image: */
  for (i = n = 0; i < nels; i++)
    n += (src_bands[i] == VAllBands) ? VImageNBands (src_images[i]) : 1;

  /* Check or allocate the destination image: */
  result = VSelectDestImage ("VCombineBands", dest, n,
			     VImageNRows (src), VImageNColumns (src), 
			     VPixelRepn (src));
  if (! result)
    return NULL;

  /* Copy each source band into the destination image: */
  for (i = n = 0; i < nels; i++) {
    if (! VCopyBand (src_images[i], src_bands[i], result, n)) {
      if (result != dest)
	VDestroyImage (result);
      return NULL;
    }
    n += (src_bands[i] == VAllBands) ? VImageNBands (src_images[i]) : 1;
  }
  return result;
}


/*
 *  VCombineBandsVa
 *
 *  A varargs version of VCombineBands. It is called by:
 *
 *	dest = VCombineBandsVa (dest, src_image1, src_band1, ...,
 *				(VImage) NULL);
 */

VImage VCombineBandsVa (VImage dest, ...)
{
  va_list args;
  VImage src, result;
  int nbands;
  VBand src_band, dest_band;

  /* Count the number of bands to be combined: */
  va_start (args, dest);
  for (nbands = 0; (src = va_arg (args, VImage)); nbands +=
	 (va_arg (args, VBand) == VAllBands) ? VImageNBands (src) : 1) ;
  va_end (args);

  /* Check or allocate the destination image: */
  va_start (args, dest);
  src = va_arg (args, VImage);
  va_end (args);
  result = VSelectDestImage ("VCombineBandsVa", dest, nbands,
			     VImageNRows (src), VImageNColumns (src),
			     VPixelRepn (src));
  if (! result)
    return NULL;

  /* Copy each source band into the destination image: */
  va_start (args, dest);
  for (dest_band = 0; (src = va_arg (args, VImage)); ) {
    src_band = va_arg (args, VBand);
    if (! VCopyBand (src, src_band, result, dest_band)) {
      if (result != dest)
	VDestroyImage (result);
      return NULL;
    }
    dest_band += (src_band == VAllBands) ? VImageNBands (src) : 1;
  }
  va_end (args);
  return result;
}


/*
 *  VSelectDestImage
 *
 *  Check that a destination image provided for an operation has the
 *  appropriate number of rows, columns and bands, and a pixel representation.
 *  Or, if no destination image is provided, create one.
 */

/*  There are two ways to use this routine. If an operation is such that
 *  it can be carried out with a destination image that is the same as
 *  the source image, one follows this procedure:
 *
 *	VImage result;
 *
 *	result = VSelectDestImage (...);
 *	if (! result)
 *	    return NULL;
 *
 *	On successful completion:
 *      VCopyImageAttrs (src, result);
 *	return result;
 *
 *	On failure:
 *	if (result != dest)
 *	    VDestroyImage (result);
 *	return NULL;
 *
 *  And if an operation *cannot* be carried out with a destination image
 *  that is the same as the source image, one does:
 *
 *	VImage result;
 *
 *	result = VSelectDestImage (...);	use or create dest image
 *	if (! result)
 *	    return NULL;
 *	if (src == result)
 *	    result = VCreateImage (...);	allocate a work image
 *
 *	On successful completion:
 *	if (src == dest) {
 *	    VCopyImagePixels (result, dest, VAllBands); move work to dest
 *	    VDestroyImage (result);
 *	    return dest;
 *	} else {
 *	    VCopyImageAttrs (src, result);
 *	    return result;
 *      }
 *
 *	On failure:
 *	if (result != dest)
 *	    VDestroyImage (result);
 *	return NULL;
 */

VImage VSelectDestImage (VStringConst routine, VImage dest,
			 int nbands, int nrows, int ncolumns,
			 VRepnKind pixel_repn)
{
  /* If no destination image was specified, allocate one: */
  if (! dest)
    return VCreateImage (nbands, nrows, ncolumns, pixel_repn);

    /* Otherwise check that the destination provided has the appropriate
       characteristics: */
  if (VImageNBands (dest) != nbands) {
    VWarning ("%s: Destination image has %d bands; %d expected",
	      routine, VImageNBands (dest), nbands);
    return NULL;
  }
  if (VImageNRows (dest) != nrows) {
    VWarning ("%s: Destination image has %d rows; %d expected",
	      routine, VImageNRows (dest), nrows);
    return NULL;
  }
  if (VImageNColumns (dest) != ncolumns) {
    VWarning ("%s: Destination image has %d columns; %d expected",
	      routine, VImageNColumns (dest), ncolumns);
    return NULL;
  }
  if (VPixelRepn (dest) != pixel_repn) {
    VWarning ("%s: Destination image has %s pixels; %s expected", routine, 
	      VPixelRepnName (dest), VRepnName (pixel_repn));
    return NULL;
  }
  return dest;
}


/*
 *  VSelectBand
 *
 *  Check a band specification and use it to determine the number and
 *  address of a block of pixels.
 */

VBoolean VSelectBand (VStringConst routine, VImage image, VBand band,
		      int *npixels, VPointer *first_pixel)
{
  if (band == VAllBands) {
    if (npixels)
      *npixels = VImageNPixels (image);
    if (first_pixel)
      *first_pixel = VImageData (image);
  } else if (band >= 0 && band < VImageNBands (image)) {
    if (npixels)
      *npixels = VImageNRows (image) * VImageNColumns (image);
    if (first_pixel)
      *first_pixel = image->band_index[band][0];
  } else {
    VWarning ("%s: Band %d referenced in image of %d band(s)",
	      routine, band, VImageNBands (image));
    return FALSE;
  }
  return TRUE;
}


/*
 *  VImageFrameInterp, VImageViewpointInterp,
 *  VImageColorInterp, VImageComponentInterp
 *
 *  Routines for accessing an image's band interpretation information.
 *
 *  Each routine returns a VBandInterpXxx constant describing how
 *  image bands are to be interpreted at a particular level of the band
 *  hierarchy.
 *
 *  If that level's dimension is 1 and there is no band interpretation
 *  attribute for the level, VBandInterpNone is returned.
 *
 *  If the dimension is >1 and there is no attribute, VBandInterpOther is
 *  returned.
 *
 *  If an image's band interpretation information is inconsistent
 *  (e.g., the color_interp attribute says RGB but ncolors is 2) then
 *  VWarning is called and VBandInterpOther is returned.
 */

VBandInterp VImageFrameInterp (VImage image)
{
  VLong interp;
  VGetAttrResult result;

  if (VImageNBands (image) !=
      (VImageNFrames (image) * VImageNViewpoints (image) *
       VImageNColors (image) * VImageNComponents (image)))
    VWarning ("VImageFrameInterp: No. bands (%d) conflicts with no. "
	      "of frames, etc. (%d %d %d %d)",
	      VImageNBands (image), VImageNFrames (image),
	      VImageNViewpoints (image), VImageNColors (image),
	      VImageNComponents (image));

  if (! VImageAttrList (image) ||
      (result =
       VGetAttr (VImageAttrList (image), VFrameInterpAttr,
		 VBandInterpDict, VLongRepn, & interp)) == VAttrMissing)
    return VImageNFrames (image) > 1 ? VBandInterpOther : VBandInterpNone;

  if (result == VAttrBadValue)
    return VBandInterpOther;

  switch (interp) {

  }
  return VBandInterpOther;
}

VBandInterp VImageViewpointInterp (VImage image)
{
  VLong interp;
  VGetAttrResult result;

  if (VImageNBands (image) !=
      (VImageNFrames (image) * VImageNViewpoints (image) *
       VImageNColors (image) * VImageNComponents (image)))
    VWarning ("VImageViewpointInterp: No. bands (%d) conflicts with no. "
	      "of frames, etc. (%d %d %d %d)",
	      VImageNBands (image), VImageNFrames (image),
	      VImageNViewpoints (image), VImageNColors (image),
	      VImageNComponents (image));
    
  if (! VImageAttrList (image) ||
      (result =
       VGetAttr (VImageAttrList (image), VViewpointInterpAttr,
		 VBandInterpDict, VLongRepn, & interp)) == VAttrMissing)
    return VImageNViewpoints (image) > 1 ?
      VBandInterpOther : VBandInterpNone;

  if (result == VAttrBadValue)
    return VBandInterpOther;

  switch (interp) {

  case VBandInterpStereoPair:
    if (VImageNViewpoints (image) != 2) {
      VWarning ("VBandViewpointInterp: "
		"Stereo-pair image has %d viewpoint dimension(s)",
		VImageNViewpoints (image));
      return VBandInterpOther;
    }
    return VBandInterpStereoPair;
  }
  return VBandInterpOther;
}

VBandInterp VImageColorInterp (VImage image)
{
  VLong interp;
  VGetAttrResult result;

  if (VImageNBands (image) !=
      (VImageNFrames (image) * VImageNViewpoints (image) *
       VImageNColors (image) * VImageNComponents (image)))
    VWarning ("VImageColorInterp: No. bands (%d) conflicts with no. "
	      "of frames, etc. (%d %d %d %d)",
	      VImageNBands (image), VImageNFrames (image),
	      VImageNViewpoints (image), VImageNColors (image),
	      VImageNComponents (image));
    
  if (! VImageAttrList (image) ||
      (result =
       VGetAttr (VImageAttrList (image), VColorInterpAttr,
		 VBandInterpDict, VLongRepn, & interp)) == VAttrMissing)
    return VImageNColors (image) > 1 ? VBandInterpOther : VBandInterpNone;

  if (result == VAttrBadValue)
    return VBandInterpOther;

  switch (interp) {

  case VBandInterpRGB:
    if (VImageNColors (image) != 3) {
      VWarning ("VBandColorInterp: RGB image has %d color dimension(s)",
		VImageNColors (image));
      return VBandInterpOther;
    }
    return VBandInterpRGB;
  }
  return VBandInterpOther;
}

VBandInterp VImageComponentInterp (VImage image)
{
  VLong interp;
  VGetAttrResult result;

  if (VImageNBands (image) !=
      (VImageNFrames (image) * VImageNViewpoints (image) *
       VImageNColors (image) * VImageNComponents (image)))
    VWarning ("VImageComponentInterp: No. bands (%d) conflicts with no. "
	      "of frames, etc. (%d %d %d %d)",
	      VImageNBands (image), VImageNFrames (image),
	      VImageNViewpoints (image), VImageNColors (image),
	      VImageNComponents (image));
    
  if (! VImageAttrList (image) ||
      (result =
       VGetAttr (VImageAttrList (image), VComponentInterpAttr,
		 VBandInterpDict, VLongRepn, & interp)) == VAttrMissing)
    return VImageNComponents (image) > 1 ?
      VBandInterpOther : VBandInterpNone;

  if (result == VAttrBadValue)
    return VBandInterpOther;

  switch (interp) {

  case VBandInterpComplex:
    if (VImageNComponents (image) != 2) {
      VWarning ("VBandColorInterp: Complex image has %d component(s)",
		VImageNComponents (image));
      return VBandInterpOther;
    }
    return VBandInterpComplex;

  case VBandInterpGradient:
    if (VImageNComponents (image) > 3) {
      VWarning ("VBandColorInterp: Gradient image has %d component(s)",
		VImageNComponents (image));
      return VBandInterpOther;
    }
    return VBandInterpGradient;

  case VBandInterpIntensity:
    if (VImageNComponents (image) > 1) {
      VWarning ("VBandColorInterp: Intensity image has %d component(s)",
		VImageNComponents (image));
      return VBandInterpOther;
    }
    return VBandInterpIntensity;

  case VBandInterpOrientation:
    if (VImageNComponents (image) > 1) {
      VWarning ("VBandColorInterp: "
		"Orientation image has %d component(s)",
		VImageNComponents (image));
      return VBandInterpOther;
    }
    return VBandInterpOrientation;
  }
  return VBandInterpOther;
}


/*
 *  VSetBandInterp
 *
 *  Set all of an image's standard band interpretation attributes.
 */

VBoolean VSetBandInterp (VImage image,
			 VBandInterp frame_interp, int nframes,
			 VBandInterp viewpoint_interp, int nviewpoints,
			 VBandInterp color_interp, int ncolors,
			 VBandInterp component_interp, int ncomponents)
{
  VBoolean result = TRUE;
  VString str;

  if (VImageNBands (image) !=
      nframes * nviewpoints * ncolors * ncomponents) {
    VWarning ("VSetBandInterp: No. bands (%d) conflicts with no. "
	      "of frames, etc. (%d %d %d %d)", VImageNBands (image),
	      nframes, nviewpoints, ncolors, ncomponents);
    result = FALSE;
  }

  if (frame_interp == VBandInterpNone)
    result &= VExtractAttr (VImageAttrList (image), VFrameInterpAttr,
			    NULL, VStringRepn, & str, FALSE);
  else VSetAttr (VImageAttrList (image), VFrameInterpAttr,
		 VBandInterpDict, VLongRepn, (VLong) frame_interp);
  VImageNFrames (image) = nframes;

  if (viewpoint_interp == VBandInterpNone)
    result &= VExtractAttr (VImageAttrList (image), VViewpointInterpAttr,
			    NULL, VStringRepn, & str, FALSE);
  else VSetAttr (VImageAttrList (image), VViewpointInterpAttr,
		 VBandInterpDict, VLongRepn, (VLong) viewpoint_interp);
  VImageNViewpoints (image) = nviewpoints;

  if (color_interp == VBandInterpNone)
    result &= VExtractAttr (VImageAttrList (image), VColorInterpAttr,
			    NULL, VStringRepn, & str, FALSE);
  else VSetAttr (VImageAttrList (image), VColorInterpAttr,
		 VBandInterpDict, VLongRepn, (VLong) color_interp);
  VImageNColors (image) = ncolors;

  if (component_interp == VBandInterpNone)
    result &= VExtractAttr (VImageAttrList (image), VComponentInterpAttr,
			    NULL, VStringRepn, & str, FALSE);
  else VSetAttr (VImageAttrList (image), VComponentInterpAttr,
		 VBandInterpDict, VLongRepn, (VLong) component_interp);
  VImageNComponents (image) = ncomponents;

  return result;
}


/*
 *  VReadImages
 *
 *  Read a Vista data file, extract the images from it, and return a list
 *  of them.
 */

int VReadImages (FILE *file, VAttrList *attributes, VImage **images)
{
  return VReadObjects (file, VImageRepn, attributes, (VPointer **) images);
}


/*
 *  VWriteImages
 *
 *  Write a list of images to a Vista data file.
 */

VBoolean VWriteImages (FILE *file, VAttrList attributes,
		       int nimages, VImage images[])
{
  return VWriteObjects (file, VImageRepn, attributes, nimages,
			(VPointer *) images);
}
