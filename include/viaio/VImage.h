/*
 *  $Id: VImage.h 726 2004-03-08 13:12:45Z lohmann $
 *
 *  Definitions associated with images: their representation in files and
 *  in memory, and operations that can be performed with them.
 */

#ifndef V_VImage_h
#define V_VImage_h 1

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

/* From the standard C library: */
#include <stdio.h>

/* For definition of pid_t: */
#include <sys/types.h>


#ifdef __cplusplus
extern "C" {
#endif


/*
 *  The representation of an image in memory.
 */

/* Description of an image: */
typedef struct V_ImageRec {
    int nbands;				/* number of bands */
    int nrows;				/* number of rows */
    int ncolumns;			/* number of columns */
    VRepnKind pixel_repn;		/* representation of pixel values */
    unsigned long flags;		/* various flags */
    VAttrList attributes;		/* list of other image attributes */
    VPointer data;			/* array of image pixel values */
    VPointer *row_index;		/* ptr to first pixel of each row */
    VPointer **band_index;		/* ptr to first row of each band */
    int nframes;			/* number of motion frames */
    int nviewpoints;			/* number of camera viewpoints */
    int ncolors;			/* number of color channels */
    int ncomponents;			/* number of vector components */
} VImageRec;

/* Codes for flags: */
enum {
    VImageSingleAlloc = 0x01		/* one free() releases everything */
};


/*
 *  Macros for accessing image in memory.
 */

/* The number of bands, rows, columns, etc. in an image: */
#define VImageNBands(image)	((image)->nbands)
#define VImageNRows(image)	((image)->nrows)
#define VImageNColumns(image)	((image)->ncolumns)
#define VImageNFrames(image)	((image)->nframes)
#define VImageNViewpoints(image) ((image)->nviewpoints)
#define VImageNColors(image)	((image)->ncolors)
#define VImageNComponents(image) ((image)->ncomponents)

/* An image's pixel representation: */
#define VPixelRepn(image)	((image)->pixel_repn)

/* An image's pixel data: */
#define VImageData(image)	((image)->data)

/* An image's list of attributes: */
#define VImageAttrList(image)	((image)->attributes)

/* The number of pixels in an image: */
#define VImageNPixels(image) \
    ((image)->nbands * (image)->nrows * (image)->ncolumns)

/* The size, in bytes, of a pixel on the present architecture: */
#define VPixelSize(image)	(VRepnSize ((image)->pixel_repn))

/* The minimum size, in bits, needed to represent a pixel: */
#define VPixelPrecision(image)  (VRepnPrecision ((image)->pixel_repn))

/* The name of an image's pixel representation: */
#define VPixelRepnName(image)	(VRepnName ((image)->pixel_repn))

/* The minimum and maximum legal pixel values: */
#define VPixelMinValue(image)	(VRepnMinValue ((image)->pixel_repn))
#define VPixelMaxValue(image)	(VRepnMaxValue ((image)->pixel_repn))

/* The size, in bytes, of an entire image on the present architecture: */
#define VImageSize(image)	(VImageNPixels(image) * VPixelSize(image))

/* A pointer to a specified pixel: */
#define VPixelPtr(image, band, row, column) \
    ((VPointer) ((char *) ((image)->band_index[band][row]) + \
		 (column) * VPixelSize (image)))

/* A specified pixel of a specified type: */
#define VPixel(image, band, row, column, type) \
    (* ((type *) (image)->band_index[band][row] + (column)))

/* A pointer to a three-dimensional array of pixels: */
#define VPixelArray(image, type) \
    ((type ***) (image)->band_index)

/* The band corresponding to a particular frame, viewpoint, etc.: */
#define VBandIndex(image, frame, viewpoint, color, component) \
    (((((frame) * (image)->nviewpoints + (viewpoint)) * (image)->ncolors + \
       (color)) * (image)->ncomponents) + (component))

/* Test whether two images have the same representational range: */
#define VSameImageRange(image1, image2)	\
    ((image1)->nbands == (image2)->nbands && \
     (image1)->nrows == (image2)->nrows && \
     (image1)->ncolumns == (image2)->ncolumns && \
     (image1)->pixel_repn == (image2)->pixel_repn)

/* Test whether two images have the same size: */
#define VSameImageSize(image1, image2) \
    ((image1)->nbands == (image2)->nbands && \
     (image1)->nrows == (image2)->nrows && \
     (image1)->ncolumns == (image2)->ncolumns)


/*
 *  Attributes used to represent an image.
 */

#define VColorInterpAttr	"color_interp"
#define VComponentInterpAttr	"component_interp"
#define VFrameInterpAttr	"frame_interp"
#define VNBandsAttr		"nbands"
#define VNColorsAttr		"ncolors"
#define VNComponentsAttr	"ncomponents"
#define VNFramesAttr		"nframes"
#define VNViewpointsAttr	"nviewpoints"
#define VPixelAspectRatioAttr	"pixel_aspect_ratio"
#define VViewpointInterpAttr	"viewpoint_interp"

/* Values of band interpretation attributes: */
typedef enum {

    /* Used at all levels of the hierarchy (frame, viewpoint, color, ...) */
    VBandInterpNone,		/* no interpretation specified */
    VBandInterpOther,		/* unknown interpretation specified */

    /* Describing motion frames: */

    /* Describing camera viewpoints: */
    VBandInterpStereoPair,

    /* Describing color channels: */
    VBandInterpRGB,

    /* Describing vector components: */
    VBandInterpComplex,
    VBandInterpGradient,
    VBandInterpIntensity,
    VBandInterpOrientation
} VBandInterp;

extern VDictEntry VBandInterpDict[];


/*
 *  Miscellaneous things.
 */

/* For specifying a band: */
typedef int VBand;
#define VAllBands	-1		/* all bands */

/* For specifying how to pad when convolving: */
typedef enum {
    VConvolvePadNone,			/* set border pixels to zero */
    VConvolvePadZero,			/* pad with zero */
    VConvolvePadBorder,			/* pad with border pixel values */
    VConvolvePadWrap,			/* pad with wrapped around values */
    VConvolvePadTrim			/* don't pad, trim destination image */
} VConvolvePadMethod;

/* For specifying which type of image filtering is to be performed: */
typedef enum {
    VFilterGaussian,			/* 2D Gaussian smoothing */
    VFilterGaussianDx,			/* d / dX */
    VFilterGaussianDy			/* d / dY */
} VFilterKind;

/* Keywords for representing kind of operation: */
typedef enum {

    /* Unary operations: */
    VImageOpAbs,
    VImageOpExp,
    VImageOpLog,
    VImageOpNot,
    VImageOpSqrt,
    VImageOpSquare,

    /* Binary (dyadic) operations: */
    VImageOpAdd,
    VImageOpAnd,
    VImageOpDist,
    VImageOpDiv,
    VImageOpMax,
    VImageOpMin,
    VImageOpMult,
    VImageOpOr,
    VImageOpSub,
    VImageOpXor
} VImageOpKind;


/*
 *  Declarations of library routines.
 */

/* From Adjust.c: */

extern VImage VAdjustImage (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */,
    double		/* brightness */,
    double		/* contrast */
#endif
);

/* From Canny.c: */

extern VImage VCanny (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */,
    double		/* sigma */,
    VImage *		/* orientation */,
    double *		/* noise */
#endif
);

/* From Complex.c: */

extern VImage VBuildComplexImage (
#if NeedFunctionPrototypes
    VImage		/* r_src */,
    VBand		/* r_band */,
    VImage		/* i_src */,
    VBand		/* i_band */,
    VImage		/* dest */,
    VRepnKind		/* pixel_repn */
#endif
);

extern VImage VImageMagnitude (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */
#endif
);

extern VImage VImagePhase (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */
#endif
);

/* From ConvertC.c: */

extern VImage VConvertImageCopy (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */,
    VRepnKind		/* pixel_repn */
#endif
);

/* From ConvertL.c: */

extern VImage VConvertImageLinear (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */,
    VRepnKind		/* pixel_repn */,
    double		/* a */,
    double		/* b */
#endif
);

/* From ConvertO.c: */

extern VImage VConvertImageOpt (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */,
    VRepnKind		/* pixel_repn */,
    int			/* method */
#endif
);

/* From ConvertR.c: */

extern VImage VConvertImageRange (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */,
    VRepnKind		/* pixel_repn */
#endif
);

/* From Convolve.c: */

extern VImage VConvolveImage (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */,
    VImage		/* mask */,
    VConvolvePadMethod	/* pad_method */,
    int			/* shift */
#endif
);

/* From ConvolveSep.c: */

extern VImage VConvolveImageSep (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */,
    VImage [3]		/* masks */,
    VConvolvePadMethod [3] /* pad_methods */,
    int [3]		/* shifts */
#endif
);

/* From Crop.c: */

extern VImage VCropImage (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */,
    int			/* top */,
    int			/* left */,
    int			/* height */,
    int			/* width */
#endif
);

/* From Dither.c: */

extern VBoolean VDither (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */,
    int			/* top */,
    int			/* left */,
    int			/* height */,
    int			/* width */,
    int	[]		/* nvalues */,
    VBooleanPromoted	/* absolute */
#endif
);

/* From Fft.c: */

extern VImage VImageFFT (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBooleanPromoted	/* inverse */,
    double		/* norm */
#endif
);

/* In Fill.c: */

extern VBoolean VFillImage (
#if NeedFunctionPrototypes
    VImage 		/* image */,
    VBand		/* band */,
    VDoublePromoted	/* value */
#endif
);

/* From Flip.c: */

extern VImage VFlipImage (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */,
    VBooleanPromoted	/* vertical */
#endif
);

/* From ForkImageDpy.c: */

extern pid_t VForkImageDisplay (
#if NeedFunctionPrototypes
    VImage		/* image */,
    VStringConst	/* program */,
    VStringConst	/* title */
#endif
);

/* From GaussConv.c: */

extern VImage VGaussianConvolveImage (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */,
    double		/* sigma */,
    int			/* filter_size */,
    VFilterKind		/* filter_kind */
#endif
);

/* From Gradient.c: */

extern VImage VImageGradient (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */
#endif
);

/* From GradientSeq.c: */

extern VImage VImageGradientSeq (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */,
    int			/* kernel_size */,
    VDouble []		/* kernel */
#endif
);

/* From Image.c: */

extern VImage VCreateImage (
#if NeedFunctionPrototypes
    int			/* nbands */,
    int			/* nrows */,
    int			/* ncols */,
    VRepnKind		/* pixel_repn */
#endif
);

extern VImage VCreateImageLike (
#if NeedFunctionPrototypes
    VImage		/* src  */
#endif
);

extern void VDestroyImage (
#if NeedFunctionPrototypes
    VImage		/* image */
#endif
);

extern VDouble VGetPixel (
#if NeedFunctionPrototypes
    VImage		/* image */,
    int			/* band */,
    int			/* row */,
    int			/* column */
#endif
);

extern void VSetPixel (
#if NeedFunctionPrototypes
    VImage		/* image */,
    int			/* band */,
    int			/* row */,
    int			/* column */,
    VDoublePromoted	/* value */
#endif
);

extern VImage VCopyImage (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */
#endif
);

extern VImage VCopyImageAttrs (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */
#endif
);

extern VImage VCopyImagePixels (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */
#endif
);

extern VBoolean VCopyBand (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VBand		/* src_band */,
    VImage		/* dest */,
    VBand		/* dest_band */
#endif
);

extern VImage VCombineBands (
#if NeedFunctionPrototypes
    int			/* nels */,
    VImage []		/* src_images */,
    VBand []		/* src_bands */,
    VImage		/* dest */
#endif
);

extern VImage VCombineBandsVa (
#if NeedVarargsPrototypes
    VImage		/* dest */,
    ...
#endif
);

extern VImage VSelectDestImage (
#if NeedFunctionPrototypes
    VStringConst	/* routine */,
    VImage		/* dest */,
    int			/* nbands */,
    int			/* nrows */,
    int			/* ncolumns */,
    VRepnKind		/* pixel_repn */
#endif
);

extern VBoolean VSelectBand (
#if NeedFunctionPrototypes
    VStringConst	/* routine */,
    VImage		/* image */,
    VBand		/* band */,
    int *		/* npixels */,
    VPointer *		/* first_pixel */
#endif
);

extern VBandInterp VImageFrameInterp (
#if NeedFunctionPrototypes
    VImage		/* image */
#endif
);

extern VBandInterp VImageViewpointInterp (
#if NeedFunctionPrototypes
    VImage		/* image */
#endif
);

extern VBandInterp VImageColorInterp (
#if NeedFunctionPrototypes
    VImage		/* image */
#endif
);

extern VBandInterp VImageComponentInterp (
#if NeedFunctionPrototypes
    VImage		/* image */
#endif
);

extern VBoolean VSetBandInterp (
#if NeedFunctionPrototypes
    VImage		/* image */,
    VBandInterp		/* frame_interp */,
    int			/* nframes */,
    VBandInterp		/* viewpoint_interp */,
    int			/* nviewpoints */,
    VBandInterp		/* color_interp */,
    int			/* ncolors */,
    VBandInterp		/* component_interp */,
    int			/* ncomponents */
#endif
);

extern int VReadImages (
#if NeedFunctionPrototypes
    FILE *		/* file */,
    VAttrList *		/* attributes */,
    VImage **		/* images */
#endif
);

extern VBoolean VWriteImages (
#if NeedFunctionPrototypes
    FILE *		/* file */,
    VAttrList		/* attributes */,
    int			/* nimages */,
    VImage []		/* images */
#endif
);

/* From ImageDpy.c: */

extern void VImageWindowSize (
#if NeedFunctionPrototypes
    VImage		/* image */,
    int			/* max_width */,
    int			/* max_height */,
    int *		/* width */,
    int *		/* height */
#endif
);

/* From ImageToPS.c: */

extern VBoolean VImageBandToPS (
#if NeedFunctionPrototypes
    FILE *		/* f */,
    VImage		/* image */,
    VBand		/* band */
#endif
);

extern VBoolean VRGBImageToPS (
#if NeedFunctionPrototypes
    FILE *		/* f */,
    VImage		/* image */,
    VBand		/* band */
#endif
);

extern VBoolean VGradientImageToPS (
#if NeedFunctionPrototypes
    FILE *		/* f */,
    VImage		/* image */,
    VBand		/* band */,
    double		/* arrow_width */
#endif
);

/* From Invert.c: */

extern VImage VInvertImage (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */
#endif
);

extern VImage VNegateImage (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */
#endif
);

/* From Kernel.c: */

extern VBoolean VGaussianKernel (
#if NeedFunctionPrototypes
    int			/* ncoeffs */,
    VDouble [] 		/* coeffs */,
    double		/* sigma */
#endif
);

extern VBoolean VGaussianD1Kernel (
#if NeedFunctionPrototypes
    int			/* ncoeffs */,
    VDouble [] 		/* coeffs */,
    double		/* sigma */
#endif
);

extern VBoolean VGaussianSplineKernel (
#if NeedFunctionPrototypes
    int			/* ncoeffs */,
    VDouble [] 		/* coeffs */,
    double		/* sigma */
#endif
);

/* From Op.c: */

extern VImage VImageOpI (
#if NeedFunctionPrototypes
    VImage		/* src1 */,
    VImage		/* dest */,
    VBand		/* band1 */,
    VImageOpKind	/* op */,
    VImage		/* src2 */,
    VBand		/* band2 */,
    VDouble *		/* minptr */,
    VDouble *		/* maxptr */
#endif
);

extern VImage VImageOpU (
#if NeedFunctionPrototypes
    VImage		/* src1 */,
    VImage		/* dest */,
    VBand		/* band1 */,
    VImageOpKind	/* op */,
    VDouble *		/* minptr */,
    VDouble *		/* maxptr */
#endif
);

extern VImage VImageOpV (
#if NeedFunctionPrototypes
    VImage		/* src1 */,
    VImage		/* dest */,
    VBand		/* band1 */,
    VImageOpKind	/* op */,
    double		/* value */,
    VDouble *		/* minptr */,
    VDouble *		/* maxptr */
#endif
);

/* From OptFlowWLS.c: */

VImage VOptFlowWLS (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    int			/* method */,
    double		/* smoothing */,
    double		/* threshold */,
    double [2]		/* noise */,
    double		/* prior_variance */
#endif
);

/* From ReadPlain.c: */

VImage VReadPlain (
#if NeedFunctionPrototypes
    FILE *		/* plain format file */
#endif
);

/* From ReadPnm.c: */

VImage VReadPnm (
#if NeedFunctionPrototypes
    FILE *		/* pnm file */
#endif
);

/* From RGBToGray.c: */

extern VImage VRGBImageToGray (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */
#endif
);

/* From Rotate.c: */

extern VImage VRotateImage (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */,
    double		/* angle */
#endif
);

/* From Sample.c: */

extern VImage VSampleImage (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */,
    int			/* row reduction factor */,
    int			/* column reduction factor */
#endif
);

/* From Scale.c: */

VImage VScaleImage (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */,
    double		/* row scaling factor */,
    double		/* column scaling factor */
#endif
);

VImage VReduceImage (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */,
    int			/* row reduction factor */,
    int			/* column reduction factor */
#endif
);

VImage VExpandImage (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */,
    int			/* row expansion factor */,
    int			/* column expansion factor */
#endif
);

/* From Shear.c: */

extern VImage VShearImageX (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */,
    double		/* shear */
#endif
);

extern VImage VShearImageY (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */,
    double		/* shear */
#endif
);

/* From Stats.c: */

extern VBoolean VImageStats (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VBand		/* band */,
    VDouble *		/* pmin */,
    VDouble *		/* pmax */,
    VDouble *		/* pmean */,
    VDouble *		/* pvar */
#endif
);

/* From Synth.c: */

/* Later in this file: */
VBoolean VRampImage (
#if NeedFunctionPrototypes
    VImage		/* image */,
    VBand		/* band */,
    double [2]		/* origin */,
    double		/* orientation */,
    double 		/* base */,
    double 		/* slope */
#endif
);

VBoolean VSineGratingImage (
#if NeedFunctionPrototypes
    VImage		/* image */,
    VBand		/* band */,
    double [2]		/* origin */,
    double		/* amplitude */,
    double		/* base */,
    double		/* frequency */,
    double		/* orientation */,
    double		/* phase */
#endif
);

VBoolean VZonePlateImage (
#if NeedFunctionPrototypes
    VImage		/* image */,
    VBand		/* band */,
    double [2]		/* origin */,
    double		/* amplitude */,
    double		/* base */,
    double		/* frequency */,
    double		/* phase */
#endif
);

VBoolean VBinomialNoiseImage (
#if NeedFunctionPrototypes
    VImage		/* image */,
    VBand		/* band */,
    VDoublePromoted	/* v0 */,
    VDoublePromoted	/* v1 */,
    double		/* p */
#endif
);

VBoolean VNormalNoiseImage (
#if NeedFunctionPrototypes
    VImage		/* image */,
    VBand		/* band */,
    double		/* mean */,
    double 		/* std_dev */
#endif
);

VBoolean VUniformNoiseImage (
#if NeedFunctionPrototypes
    VImage		/* image */,
    VBand		/* band */,
    VDoublePromoted	/* min_value */,
    VDoublePromoted	/* max_value */
#endif
);

/* From Transpose.c: */

extern VImage VTransposeImage (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */
#endif
);

/* From UbcIff.c: */

extern VImage VReadUbcIff (
#if NeedFunctionPrototypes
    FILE *		/* f */
#endif
);

extern VBoolean VWriteUbcIff (
#if NeedFunctionPrototypes
    FILE *		/* f */,
    VImage		/* image */,
    VBand		/* band */
#endif
);

/* From ZeroC.c: */

extern VImage VZeroCrossings (
#if NeedFunctionPrototypes
    VImage		/* src */,
    VImage		/* dest */,
    VBand		/* band */
#endif
);

#ifdef __cplusplus
}
#endif

#endif /* V_VImage_h */
