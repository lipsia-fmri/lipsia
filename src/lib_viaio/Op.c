/*
 *  This file contains routines for performing pixel-wise arithmetic and
 *  logical operations on images.
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

/* From the standard C library: */
#include <math.h>

extern double rint(double);

/* Later in this file: */
static void SetClippingRange (VStringConst, VImage,
			      VDouble *, VDouble *, double *, double *);

static VStringConst bad_op_msg = "%s: Invalid operation (%d) for %s pixels";
static VStringConst bound_msg = "%s: Clipping bound out of range; set to %g";
static VStringConst clip_msg = "%s: Destination pixel value(s) clipped";


#define DoBinaryOps(type, arg2, logical_ops)				\
    {									\
	type t1, t2;							\
	type local_min = min, local_max = max;				\
									\
	switch (op) {							\
									\
	case VImageOpAdd:						\
	    DoBinaryOp (*src1_pp++ + arg2);				\
	    break;							\
									\
	case VImageOpDist:						\
	    DoBinaryOp ((t1 = *src1_pp++ - arg2) > 0 ? t1 : -t1);	\
	    break;							\
									\
	case VImageOpDiv:						\
	    for (i = 0; i < npixels; i++) {				\
		t1 = (type) *src1_pp++ / (type) arg2;			\
		Clip (t1, local_min, local_max);			\
		*res_pp++ = t1;						\
	    }								\
	    break;							\
									\
	case VImageOpMax:						\
	    DoBinaryOp (((t1 = *src1_pp++) > (t2 = arg2)) ? t1 : t2);	\
	    break;							\
									\
	case VImageOpMin:						\
	    DoBinaryOp (((t1 = *src1_pp++) < (t2 = arg2)) ? t1 : t2);	\
	    break;							\
									\
	case VImageOpMult:						\
	    DoBinaryOp (*src1_pp++ * arg2);				\
	    break;							\
									\
	case VImageOpSub:						\
	    DoBinaryOp (*src1_pp++ - arg2);				\
	    break;							\
									\
	logical_ops							\
									\
	default:							\
	    VWarning (bad_op_msg, routine, op, VPixelRepnName (src1));	\
	    goto Error;							\
	}								\
    }

#define LogicalOps(arg2)						\
	case VImageOpAnd: DoBinaryOp (*src1_pp++ & arg2); break;	\
	case VImageOpOr: DoBinaryOp (*src1_pp++ | arg2); break;		\
	case VImageOpXor: DoBinaryOp (*src1_pp++ ^ arg2); break;

#define DoBinaryOp(expr)						\
    for (i = 0; i < npixels; i++) {					\
	t1 = expr;							\
	Clip (t1, local_min, local_max);				\
	*res_pp++ = t1;							\
    }

#define DoUnaryOps(type, round)						\
    {									\
	type t;								\
	type local_min = min, local_max = max;				\
									\
	switch (op) {							\
									\
	case VImageOpAbs:						\
	    for (i = 0; i < npixels; i++) {				\
		t = *src_pp++;						\
		if (t < 0)						\
		    t = -t;						\
		Clip (t, local_min, local_max);				\
		*res_pp++ = t;						\
	    }								\
	    break;							\
									\
	case VImageOpExp: DoUnaryOp (exp, Nothing round); break;	\
									\
	case VImageOpLog: DoUnaryOp (log, Nothing round); break;	\
									\
	case VImageOpNot:						\
	    for (i = 0; i < npixels; i++) {				\
		t = ! *src_pp++;					\
		Clip (t, local_min, local_max);				\
		*res_pp++ = t;						\
	    }								\
	    break;							\
									\
	case VImageOpSqrt: DoUnaryOp (sqrt, Nothing round); break;	\
									\
	case VImageOpSquare:						\
	    for (i = 0; i < npixels; i++) {				\
		t = *src_pp++;						\
		t *= t;							\
		Clip (t, local_min, local_max);				\
		*res_pp++ = t;						\
	    }								\
	    break;							\
									\
	default:							\
	    VWarning (bad_op_msg, routine, op, VPixelRepnName (src1));	\
	    goto Error;							\
	}								\
    }

#define DoUnaryOp(fnc, round)						\
    for (i = 0; i < npixels; i++) {					\
	t = round (fnc ((double) *src_pp++));				\
	Clip (t, local_min, local_max);					\
	*res_pp++ = t;							\
    }

#define Clip(v, min, max)						\
    if (v < min) {							\
	v = min;							\
	clipped = TRUE;							\
    } else if (v > max) {						\
	v = max;							\
	clipped = TRUE;							\
    }

#define Nothing


/*
 *  VImageOpI
 *
 *  Perform a binary arithmetic or logical operation with two images.
 */

VImage VImageOpI (VImage src1, VImage dest, VBand band1, VImageOpKind op,
		  VImage src2, VBand band2, VDouble *minp, VDouble *maxp)
{
    int i, npixels, npixels1, npixels2;
    VPointer src1_pixels, src2_pixels, res_pixels;
    VImage result;
    double min, max;
    VBoolean clipped = FALSE;
    static VStringConst routine = "VImageOpI";

    /* Require the two images to have similar properties: */
    if (VImageNRows (src1) != VImageNRows (src2) ||
	VImageNColumns (src1) != VImageNColumns (src2) ||
	VPixelRepn (src1) != VPixelRepn (src2)) {
	VWarning ("%s: Source images have dissimilar properties", routine);
	return NULL;
    }

    /* We can operate on one band from each image, all bands from each
       image, or all bands from the first and one from the second: */
    if (band1 == VAllBands) {
	if (band2 == VAllBands && VImageNBands (src1) != VImageNBands (src2)) {
	    VWarning ("%s: Source images have differing numbers of bands",
		      routine);
	    return NULL;
	}
    } else {
	if (band2 != VAllBands) {
	    VWarning ("%s: Incorrect combination of band1, band2", routine);
	    return NULL;
	}
    }

    /* Locate the source pixels: */
    if (! VSelectBand (routine, src1, band1, & npixels1, & src1_pixels))
	return NULL;
    if (! VSelectBand (routine, src2, band2, & npixels2, & src2_pixels))
	return NULL;

    /* Locate the destination pixels: */
    result = VSelectDestImage (routine, dest,
			       band1 == VAllBands ? VImageNBands (src1) : 1,
			       VImageNRows (src1), VImageNColumns (src1),
			       VPixelRepn (src1));
    if (! result)
	return NULL;
    res_pixels = VPixelPtr (result, 0, 0, 0);

    /* Establish the clipping range: */
    SetClippingRange (routine, src1, minp, maxp, & min, & max);

    /* Do the operation npixels2 pixels at a time. Npixels2 represents
       either a single band of src2 pixels or all src2 pixels, depending
       on band2. */
    while (npixels1 > 0) {
	npixels = npixels2;
	switch (VPixelRepn (src1)) {

	case VBitRepn:
	    {
		VBit *src1_pp = src1_pixels, *res_pp = res_pixels;
		VBit *src2_pp = src2_pixels;
		DoBinaryOps (long, *src2_pp++, LogicalOps (*src2_pp++));
	    }
	    break;

	case VUByteRepn:
	    {
		VUByte *src1_pp = src1_pixels, *res_pp = res_pixels;
		VUByte *src2_pp = src2_pixels;
		DoBinaryOps (long, *src2_pp++, LogicalOps (*src2_pp++));
	    }
	    break;

	case VSByteRepn:
	    {
		VSByte *src1_pp = src1_pixels, *res_pp = res_pixels;
		VSByte *src2_pp = src2_pixels;
		DoBinaryOps (long, *src2_pp++, LogicalOps (*src2_pp++));
	    }
	    break;

	case VShortRepn:
	    {
		VShort *src1_pp = src1_pixels, *res_pp = res_pixels;
		VShort *src2_pp = src2_pixels;
		DoBinaryOps (long, *src2_pp++, LogicalOps (*src2_pp++));
	    }
	    break;

	case VUShortRepn:
	    {
		VUShort *src1_pp = src1_pixels, *res_pp = res_pixels;
		VUShort *src2_pp = src2_pixels;
		DoBinaryOps (long, *src2_pp++, LogicalOps (*src2_pp++));
	    }
	    break;

	case VIntegerRepn:
	    {
		VInteger *src1_pp = src1_pixels, *res_pp = res_pixels;
		VInteger *src2_pp = src2_pixels;
		DoBinaryOps (long, *src2_pp++, LogicalOps (*src2_pp++));
	    }
	    break;

	case VUIntegerRepn:
	    {
		VUInteger *src1_pp = src1_pixels, *res_pp = res_pixels;
		VUInteger *src2_pp = src2_pixels;
		DoBinaryOps (long, *src2_pp++, LogicalOps (*src2_pp++));
	    }
	    break;

	case VLongRepn:
	    {
		VLong *src1_pp = src1_pixels, *res_pp = res_pixels;
		VLong *src2_pp = src2_pixels;
		DoBinaryOps (long, *src2_pp++, LogicalOps (*src2_pp++));
	    }
	    break;

	case VULongRepn:
	    {
		VULong *src1_pp = src1_pixels, *res_pp = res_pixels;
		VULong *src2_pp = src2_pixels;
		DoBinaryOps (long, *src2_pp++, LogicalOps (*src2_pp++));
	    }
	    break;

	case VFloatRepn:
	    {
		VFloat *src1_pp = src1_pixels, *res_pp = res_pixels;
		VFloat *src2_pp = src2_pixels;
		DoBinaryOps (double, *src2_pp++, Nothing);
	    }
	    break;

	case VDoubleRepn:
	    {
		VDouble *src1_pp = src1_pixels, *res_pp = res_pixels;
		VDouble *src2_pp = src2_pixels;
		DoBinaryOps (double, *src2_pp++, Nothing);
	    }
	    break;

	default:
	    break;
	}
	src1_pixels = (char *) src1_pixels + npixels2 * VPixelSize (src1);
	res_pixels = (char *) res_pixels + npixels2 * VPixelSize (result);
	npixels1 -= npixels2;
    }

    if (clipped)
	VWarning (clip_msg, routine);

    VCopyImageAttrs (src1, result);
    return result;

Error:
    if (result != dest)
	VDestroyImage (result);
    return NULL;
}


/*
 *  VImageOpU
 *
 *  Perform a unary arithmetic or logical operation on an image.
 */

VImage VImageOpU (VImage src1, VImage dest, VBand band1, VImageOpKind op,
		  VDouble *minp, VDouble *maxp)
{
    int i, npixels;
    VPointer src_pixels, res_pixels;
    VImage result;
    double min, max;
    VBoolean clipped = FALSE;
    static VStringConst routine = "VImageOpU";

    /* Locate the source pixels: */
    if (! VSelectBand (routine, src1, band1, & npixels, & src_pixels))
	return NULL;

    /* Locate the destination pixels: */
    result = VSelectDestImage (routine, dest,
			       band1 == VAllBands ? VImageNBands (src1) : 1,
			       VImageNRows (src1), VImageNColumns (src1),
			       VPixelRepn (src1));
    if (! result)
	return NULL;
    res_pixels = VPixelPtr (result, 0, 0, 0);

    /* Establish the clipping range: */
    SetClippingRange (routine, src1, minp, maxp, & min, & max);

    /* Do the operation: */
    switch (VPixelRepn (src1)) {

    case VBitRepn:
	{
	    VBit *src_pp = src_pixels, *res_pp = res_pixels;
	    DoUnaryOps (long, rint);
	}
	break;

    case VUByteRepn:
	{
	    VUByte *src_pp = src_pixels, *res_pp = res_pixels;
	    DoUnaryOps (long, rint);
	}
	break;

    case VSByteRepn:
	{
	    VSByte *src_pp = src_pixels, *res_pp = res_pixels;
	    DoUnaryOps (long, rint);
	}
	break;

    case VShortRepn:
	{
	    VShort *src_pp = src_pixels, *res_pp = res_pixels;
	    DoUnaryOps (long, rint);
	}
	break;

    case VUShortRepn:
	{
	    VUShort *src_pp = src_pixels, *res_pp = res_pixels;
	    DoUnaryOps (long, rint);
	}
	break;

    case VIntegerRepn:
	{
	    VInteger *src_pp = src_pixels, *res_pp = res_pixels;
	    DoUnaryOps (long, rint);
	}
	break;

    case VUIntegerRepn:
	{
	    VUInteger *src_pp = src_pixels, *res_pp = res_pixels;
	    DoUnaryOps (long, rint);
	}
	break;

    case VLongRepn:
	{
	    VLong *src_pp = src_pixels, *res_pp = res_pixels;
	    DoUnaryOps (long, rint);
	}
	break;

    case VULongRepn:
	{
	    VULong *src_pp = src_pixels, *res_pp = res_pixels;
	    DoUnaryOps (long, rint);
	}
	break;

    case VFloatRepn:
	{
	    VFloat *src_pp = src_pixels, *res_pp = res_pixels;
	    DoUnaryOps (double,	Nothing);
	}
	break;

    case VDoubleRepn:
	{
	    VDouble *src_pp = src_pixels, *res_pp = res_pixels;
	    DoUnaryOps (double,	Nothing);
	}
	break;

    default:
	break;
    }

    if (clipped)
	VWarning (clip_msg, routine);

    VCopyImageAttrs (src1, result);
    return result;

Error:
    if (result != dest)
	VDestroyImage (result);
    return NULL;
}


/*
 *  VImageOpV
 *
 *  Perform a binary arithmetic or logical operation with an image and a
 *  scalar.
 */

VImage VImageOpV (VImage src1, VImage dest, VBand band1, VImageOpKind op,
		  VDouble value, VDouble *minp, VDouble *maxp)
{
    int i, npixels;
    VPointer src_pixels, res_pixels;
    VImage result;
    long lvalue = 0;
    double min, max, fvalue = 0.0;
    VBoolean clipped = FALSE;
    static VStringConst routine = "VImageOpV";

    /* Locate the source pixels: */
    if (! VSelectBand (routine, src1, band1, & npixels, & src_pixels))
	return NULL;

    /* Locate the destination pixels: */
    result = VSelectDestImage (routine, dest,
			       band1 == VAllBands ? VImageNBands (src1) : 1,
			       VImageNRows (src1), VImageNColumns (src1),
			       VPixelRepn (src1));
    if (! result)
	return NULL;
    res_pixels = VPixelPtr (result, 0, 0, 0);

    /* Establish the clipping range: */
    SetClippingRange (routine, src1, minp, maxp, & min, & max);

    /* Convert value to the representation we'll need for computation: */
    if (VIsFloatPtRepn (VPixelRepn (src1)))
	fvalue = value;
    else lvalue = value;

    /* Do the operation: */
    switch (VPixelRepn (src1)) {

    case VBitRepn:
	{
	    VBit *src1_pp = src_pixels, *res_pp = res_pixels;
	    DoBinaryOps (long, lvalue, LogicalOps (lvalue));
	}
	break;

    case VUByteRepn:
	{
	    VUByte *src1_pp = src_pixels, *res_pp = res_pixels;
	    DoBinaryOps (long, lvalue, LogicalOps (lvalue));
	}
	break;

    case VSByteRepn:
	{
	    VSByte *src1_pp = src_pixels, *res_pp = res_pixels;
	    DoBinaryOps (long, lvalue, LogicalOps (lvalue));
	}
	break;

    case VShortRepn:
	{
	    VShort *src1_pp = src_pixels, *res_pp = res_pixels;
	    DoBinaryOps (long, lvalue, LogicalOps (lvalue));
	}
	break;

    case VUShortRepn:
	{
	    VUShort *src1_pp = src_pixels, *res_pp = res_pixels;
	    DoBinaryOps (long, lvalue, LogicalOps (lvalue));
	}
	break;

    case VIntegerRepn:
	{
	    VInteger *src1_pp = src_pixels, *res_pp = res_pixels;
	    DoBinaryOps (long, lvalue, LogicalOps (lvalue));
	}
	break;

    case VUIntegerRepn:
	{
	    VUInteger *src1_pp = src_pixels, *res_pp = res_pixels;
	    DoBinaryOps (long, lvalue, LogicalOps (lvalue));
	}
	break;

    case VLongRepn:
	{
	    VLong *src1_pp = src_pixels, *res_pp = res_pixels;
	    DoBinaryOps (long, lvalue, LogicalOps (lvalue));
	}
	break;

    case VULongRepn:
	{
	    VULong *src1_pp = src_pixels, *res_pp = res_pixels;
	    DoBinaryOps (long, lvalue, LogicalOps (lvalue));
	}
	break;

    case VFloatRepn:
	{
	    VFloat *src1_pp = src_pixels, *res_pp = res_pixels;
	    DoBinaryOps (double, fvalue, Nothing);
	}
	break;

    case VDoubleRepn:
	{
	    VDouble *src1_pp = src_pixels, *res_pp = res_pixels;
	    DoBinaryOps (double, fvalue, Nothing);
	}
	break;

    default:
	break;
    }

    if (clipped)
	VWarning (clip_msg, routine);

    VCopyImageAttrs (src1, result);
    return result;

Error:
    if (result != dest)
	VDestroyImage (result);
    return NULL;
}


/*
 *  SetClippingRange
 *
 *  Examine src, minp, and maxp parameters, and use them to establish
 *  actual clipping range.
 */

static void SetClippingRange (VStringConst routine, VImage src,
			      VDouble *minp, VDouble *maxp,
			      double *min, double *max)
{
    if (VIsFloatPtRepn (VPixelRepn (src))) {
	*min = -HUGE_VAL;		/* - Inf */
	*max = HUGE_VAL;		/* + Inf */
    } else {
	*min = VPixelMinValue (src);
	*max = VPixelMaxValue (src);
    }

    if (minp) {
	if (*minp < *min)
	    VWarning (bound_msg, routine, *min);
	else *min = *minp;
    }
    if (maxp) {
	if (*maxp > *max)
	    VWarning (bound_msg, routine, *max);
	else *max = *maxp;
    }
    if (*min > *max)
	VWarning ("%s: Clipping range is empty", routine);
}
