/*
 *  $Id: ImageType.c 3177 2008-04-01 14:47:24Z karstenm $
 *
 *  This file contains methods for the image (VImage) type.
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


/*
 *  Table of methods.
 */

/* Later in this file: */
static VCopyMethod VImageCopyMethod;
static VDecodeMethod VImageDecodeMethod;
static VEncodeAttrMethod VImageEncodeAttrMethod;
static VEncodeDataMethod VImageEncodeDataMethod;

/* Used in Type.c to register this type: */
VTypeMethods VImageMethods = {
  VImageCopyMethod,			/* copy a VImage */
  (VDestroyMethod *) VDestroyImage,	/* destroy a VImage */
  VImageDecodeMethod,			/* decode a VImage's value */
  VImageEncodeAttrMethod,		/* encode a VImage's attr list */
  VImageEncodeDataMethod		/* encode a VImage's binary data */
};


/*
 *  VImageCopyMethod
 *
 *  The "copy" method registered for the "image" type.
 *  Copy a VImage object.
 */

static VPointer VImageCopyMethod (VPointer value)
{
  return VCopyImage ((VImage) value, NULL, VAllBands);
}


/*
 *  VImageDecodeMethod
 *
 *  The "decode" method registered for the "image" type.
 *  Convert an attribute list plus binary data to a VImage object.
 */

static VPointer VImageDecodeMethod (VStringConst name, VBundle b)
{
  VImage image;
  VLong nbands, nrows, ncolumns, pixel_repn;
  VLong nframes, nviewpoints, ncolors, ncomponents;
  VAttrList list;
  size_t length;

#define Extract(name, dict, locn, required)	\
	VExtractAttr (b->list, name, dict, VLongRepn, & locn, required)

  /* Extract the number of bands, rows, columns, pixel repn, etc.: */
  nbands = nframes = nviewpoints = ncolors = ncomponents = 1;	/* defaults */
  if (! Extract (VNBandsAttr, NULL, nbands, FALSE) ||
      ! Extract (VNRowsAttr, NULL, nrows, TRUE) ||
      ! Extract (VNColumnsAttr, NULL, ncolumns, TRUE) ||
      ! Extract (VRepnAttr, VNumericRepnDict, pixel_repn, TRUE) ||
      ! Extract (VNFramesAttr, NULL, nframes, FALSE) ||
      ! Extract (VNViewpointsAttr, NULL, nviewpoints, FALSE) ||
      ! Extract (VNColorsAttr, NULL, ncolors, FALSE) ||
      ! Extract (VNComponentsAttr, NULL, ncomponents, FALSE))
    return NULL;

  /* Ensure that nbands == nframes * nviewpoints * ncolors * ncomponents.
     For backwards compatibility, set ncomponents to nbands if nbands != 1
     but nframes == nviewpoints == ncolors == ncomponents == 1. */
  if (nbands != nframes * nviewpoints * ncolors * ncomponents) {
    if (nbands != 1 && nframes == 1 && nviewpoints == 1 &&
	ncolors == 1 && ncomponents == 1)
      ncomponents = nbands;
    else {
      VWarning ("VImageDecodeMethod: %s image has inconsistent nbands",
		name);
      return NULL;
    }
  }

  /* Create an image with the specified properties: */
  if (! (image = VCreateImage ((int) nbands, (int) nrows, (int) ncolumns,
			       (VRepnKind) pixel_repn)))
    return NULL;
  VImageNFrames (image) = nframes;
  VImageNViewpoints (image) = nviewpoints;
  VImageNColors (image) = ncolors;
  VImageNComponents (image) = ncomponents;

  /* Give it whatever attributes remain: */
  list = VImageAttrList (image);
  VImageAttrList (image) = b->list;
  b->list = list;

  /* Check that the expected amount of binary data was read: */
  length = VImageNPixels (image);
  if (VPixelRepn (image) == VBitRepn)
    length = (length + 7) / 8;
  else length *= VPixelPrecision (image) / 8;
  if (length != b->length) {
    VWarning ("VImageDecodeMethod: %s image has wrong data length", name);
  Fail:   VDestroyImage (image);
    return NULL;
  }

  /* Unpack the binary pixel data: */
  length = VImageSize (image);
  if (! VUnpackData (VPixelRepn (image), VImageNPixels (image),
		     b->data, VMsbFirst, & length, & VImageData (image),
		     NULL))
    goto Fail;
  return image;

#undef Extract
}


/*
 *  VImageEncodeAttrMethod
 *
 *  The "encode_attrs" method registered for the "image" type.
 *  Encode an attribute list value for a VImage object.
 */

static VAttrList VImageEncodeAttrMethod (VPointer value, size_t *lengthp)
{
  VImage image = value;
  VAttrList list;
  size_t length;

#define OptionallyPrepend(value, name)				\
	if (value != 1)							\
	    VPrependAttr (list, name, NULL, VLongRepn, (VLong) value)

  /* Temporarily prepend several attributes to the image's attribute list: */
  if ((list = VImageAttrList (image)) == NULL)
    list = VImageAttrList (image) = VCreateAttrList ();
  VPrependAttr (list, VRepnAttr, VNumericRepnDict,
		VLongRepn, (VLong) VPixelRepn (image));
  VPrependAttr (list, VNColumnsAttr, NULL,
		VLongRepn, (VLong) VImageNColumns (image));
  VPrependAttr (list, VNRowsAttr, NULL,
		VLongRepn, (VLong) VImageNRows (image));
  OptionallyPrepend (VImageNComponents (image), VNComponentsAttr);
  OptionallyPrepend (VImageNColors (image), VNColorsAttr);
  OptionallyPrepend (VImageNViewpoints (image), VNViewpointsAttr);
  OptionallyPrepend (VImageNFrames (image), VNFramesAttr);
  OptionallyPrepend (VImageNBands (image), VNBandsAttr);

  /* Compute the file space needed for the image's binary data: */
  length = VImageNPixels (image);
  if (VPixelRepn (image) == VBitRepn)
    length = (length + 7) / 8;
  else length *= VPixelPrecision (image) / 8;
  *lengthp = length;

  return list;

#undef OptionallyPrepend
}


/*
 *  VImageEncodeDataMethod
 *
 *  The "encode_data" method registered for the "image" type.
 *  Encode the pixel values for a VImage object.
 */

static VPointer VImageEncodeDataMethod (VPointer value, VAttrList list,
					size_t length, VBoolean *free_itp)
{
  VImage image = value;
  VAttrListPosn posn;
  size_t len;
  VPointer ptr;

  /* Remove the attributes prepended by the VImageEncodeAttrsMethod: */
  for (VFirstAttr (list, & posn);
       strcmp (VGetAttrName (& posn), VRepnAttr) != 0;
       VDeleteAttr (& posn)) ;
  VDeleteAttr (& posn);

  /* Pack and return pixel data: */
  if (! VPackData (VPixelRepn (image), VImageNPixels (image),
		   VImageData (image), VMsbFirst, & len, & ptr, free_itp))
    return NULL;
  if (len != length)
    VError ("VImageEncodeDataMethod: Encoded data has unexpected length");
  return ptr;
}
