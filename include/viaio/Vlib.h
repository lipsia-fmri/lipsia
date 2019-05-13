/*
 *  $Id: Vlib.h,v 1.1.1.1 2004/03/08 13:12:46 lohmann Exp $
 *
 *  Definitions associated with the Vista library of vision software.
 */

#ifndef V_Vlib_h
#define V_Vlib_h 1

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


/* For portability: */
#include "viaio/Xfuncproto.h"

/* Get definition of size_t and NULL: */
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __APPLE__
  #include <stdio.h>
  FILE *fmemopen(void *buf, size_t size, const char *mode);
#endif

/*
 *  Basic constants.
 */

#ifndef FALSE
#define FALSE 0
#define TRUE 1
#endif


#define TINY 1.0e-6
#define NO_MINVAL -1.0e+6


/* function definition for writing history entry */
extern char *GetLipsiaName(char *prgname);


/*
 *  Basic types.
 */

/* Numeric types: */
/* (These definitions may be platform-specific.) */
typedef char VBit;			             /* 0 or 1 */
typedef double VDouble;		        	/* >= 64-bit IEEE floating point */
typedef float VFloat;			          /* >= 32-bit IEEE floating point */

typedef int16_t VShort;			 /* >= 16-bit signed integer */
typedef uint16_t VUShort;		 /* >= 16-bit unsigned integer */
typedef int32_t VInteger;		 /* >= 32-bit signed integer */
typedef uint32_t VUInteger;		 /* >= 32-bit unsigned integer */
typedef int64_t VLong;	                 /* >= 64-bit signed integer */
typedef uint64_t VULong;		 /* >= 64-bit unsigned integer */


#if __STDC__ 
typedef signed char VSByte;		/* integer in [-128,127] */
#else
typedef char VSByte;			/* integer in [-128,127] */
#endif
typedef unsigned char VUByte;		/* integer in [0,255] */




/* Macros for generating constants of particular numeric types: */
/* (These definitions may be platform-specific.) */
#define VBitConst(c)	(c)
#define VUByteConst(c)	(c)
#define VSByteConst(c)	(c)
#define VShortConst(c)	(c)
#define VLongConst(c)	(c ## l)
#if __STDC__
#define VFloatConst(c)	(c ## f)
#define VDoubleConst(c) (c)
#else
#define VFloatConst(c)	((VFloat) c)
#define VDoubleConst(c)	((VDouble) c)
#define VIntegerConst(c)  ((VInteger) c)
#define VUIntegerConst(c)  ((VUInteger) c)
#define VUShortConst(c)	((VUShort) c)
#define VULongConst(c)	((VULong) c)
#endif

/* Miscellaneous types: */
/* (These definitions may be platform-specific.) */
typedef char VBoolean;			/* TRUE or FALSE */
#if __STDC__ || defined(__cplusplus) || defined(c_plusplus) 
typedef void *VPointer;			/* generic pointer */
typedef const char *VStringConst;	/* null-terminated string constant */
#else
typedef char *VPointer;			/* generic pointer */
typedef char *VStringConst;		/* null-terminated string constant */
#endif
typedef char *VString;			/* null-terminated string */

/* What they are promoted to, as arguments, on the various platforms: */
/* (These definitions may be platform-specific.) */
typedef int VBitPromoted;
typedef int VBooleanPromoted;
typedef double VDoublePromoted;
typedef double VFloatPromoted;
typedef long VLongPromoted;
typedef long VULongPromoted;
typedef int VSBytePromoted;
typedef int VShortPromoted;
typedef int VUShortPromoted;
typedef int VIntegerPromoted;
typedef int VUIntegerPromoted;
typedef int VUBytePromoted;

/* Standard object types: */
typedef struct V_EdgesRec *VEdges;	/* edge set */
typedef struct V_ImageRec *VImage;	/* image    */
typedef struct V_VolumesRec *Volumes;   /* volumes  */

/* Codes for referring to representations: */
typedef enum {
    VUnknownRepn,

    /* Integer numbers: */
    VBitRepn,				/* 1-bit integer, [0, 1] */
    VUByteRepn,				/* 8-bit integer, [0, 255] */
    VSByteRepn,				/* 8-bit integer, [-128, 127] */
    VShortRepn,				/* 16-bit integer, [-32768, 32767] */
    VLongRepn,				/* 64-bit integer, [-2**63, 2**63-1] */

    /* Floating point numbers: */
    VFloatRepn,				/* 32-bit IEEE floating point */
    VDoubleRepn,			/* 64-bit IEEE floating point */

    /* Miscellaneous representations: */
    VAttrListRepn,			/* attribute list */
    VBooleanRepn,			/* TRUE or FALSE */
    VBundleRepn,			/* object of named type */
    VListRepn,				/* list of opaque objects */
    VPointerRepn,			/* pointer to opaque object */
    VStringRepn,			/* null-terminated string */

    /* Standard object types: */
    VEdgesRepn,				/* edge set */
    VImageRepn,				/* image */
    
    /* New object types */
    VolumesRepn,                        /* volumes */
    VUShortRepn,                        /* unsigned short */
    VIntegerRepn,                       /* 32-bit integer */
    VUIntegerRepn,                      /* unsigned 32-bit integer */
    VULongRepn,
    VNRepnKinds				/* number of predefined types */
} VRepnKind;


/*
 *  Attribute list representation.
 */

/* Each attribute name/value pair is represented by: */
typedef struct V_AttrRec {
    struct V_AttrRec *next;		/* next in list */
    struct V_AttrRec *prev;		/* previous in list */
    VRepnKind repn;			/* rep'n of attribute value */
    VPointer value;			/* pointer to attribute value */
    char name[1];			/* beginning of name string */
} VAttrRec;

/* A list of attributes is represented by a header node: */
typedef VAttrRec *VAttrList;

/* Macros for testing an attribute list: */
#define VAttrListEmpty(l) ((l) == NULL || (l)->next == NULL)

/* Position within a list of attributes: */
typedef struct {
    VAttrList list;			/* the list */
    struct V_AttrRec *ptr;		/* position within the list */
} VAttrListPosn;

/* Macros for moving up and down an attribute list: */
#define VFirstAttr(l,p) ((void) ((p)->list = (l), (p)->ptr = (l)->next))
#define VLastAttr(l,p)  ((void) ((p)->list = (l), (p)->ptr = (l)->prev))
#define VAttrExists(p)  ((p)->ptr != NULL)
#define VNextAttr(p)	((void) ((p)->ptr = (p)->ptr ? (p)->ptr->next : NULL))
#define VPrevAttr(p)	((void) ((p)->ptr = (p)->ptr ? (p)->ptr->prev : NULL))

/* Macros for accessing the attribute at a particular position in an list: */
#define VGetAttrName(p) ((p)->ptr->name)
#define VGetAttrRepn(p) ((p)->ptr->repn)

/* Result of trying to retrieve an attribute's value: */
typedef enum {
    VAttrFound,				/* successfully retrieved value */
    VAttrMissing,			/* didn't find attribute */
    VAttrBadValue			/* incompatible value */
} VGetAttrResult;

/* Names of generic attributes: */
#define VCommentAttr	"comment"
#define VDataAttr	"data"
#define VHistoryAttr	"history"
#define VLengthAttr	"length"
#define VNameAttr	"name"
#define VNColumnsAttr	"ncolumns"
#define VNRowsAttr	"nrows"
#define VRepnAttr	"repn"

/* An object whose type is named but not registered: */
typedef struct {
    VAttrList list;			/* object's attribute list value */
    size_t length;			/* length of binary data */
    VPointer data;			/* pointer to binary data */
    char type_name[1];			/* beginning of object's type's name */
} VBundleRec, *VBundle;


/*
 *  Object type registry.
 */

/* Type of procedure for copying object's value: */
typedef VPointer VCopyMethod (
#if NeedFunctionPrototypes
    VPointer		/* value */
#endif
);

/* Type of procedure for destroying object's value: */
typedef void VDestroyMethod (
#if NeedFunctionPrototypes
    VPointer		/* value */
#endif
);

/* Type of procedure for decoding object's binary data: */
typedef VPointer VDecodeMethod (
#if NeedFunctionPrototypes
    VStringConst	/* name */,
    VBundle		/* bundle */
#endif
);

/* Type of procedure for encoding object's attr list: */
typedef VAttrList VEncodeAttrMethod (
#if NeedFunctionPrototypes
    VPointer		/* value */,
    size_t *		/* length */
#endif
);

/* Type of procedure for encoding object's binary data: */
typedef VPointer VEncodeDataMethod (
#if NeedFunctionPrototypes
    VPointer		/* value */,
    VAttrList		/* list */,
    size_t		/* length */,
    VBoolean *		/* free_it */
#endif
);

/* Set of methods supporting an object type: */
typedef struct {
    VCopyMethod *copy;
    VDestroyMethod *destroy;
    VDecodeMethod *decode;
    VEncodeAttrMethod *encode_attr;
    VEncodeDataMethod *encode_data;
} VTypeMethods;

/* Information about a representation: */
typedef struct {
    VStringConst name;			/* name string */
    size_t size;			/* size, in bytes */
    int precision;			/* precision, in bits */
    VDouble min_value;			/* min and max representable values */
    VDouble max_value;
    VTypeMethods *methods;		/* associated methods */
} VRepnInfoRec;

/* Table indexed by representation: */
extern VRepnInfoRec *VRepnInfo;

/* Macros for accessing information about representations: */
#define VRepnSize(repn)		(VRepnInfo[repn].size)
#define VRepnPrecision(repn)	(VRepnInfo[repn].precision)
#define VRepnName(repn)		(VRepnInfo[repn].name)
#define VRepnMinValue(repn)	(VRepnInfo[repn].min_value)
#define VRepnMaxValue(repn)	(VRepnInfo[repn].max_value)
#define VRepnMethods(repn)	(VRepnInfo[repn].methods)

#define VIsIntegerRepn(repn)	(((repn) >= VBitRepn && (repn) <= VShortRepn) || ((repn) >= VUShortRepn && (repn) <= VULongRepn))
#define VIsFloatPtRepn(repn)	((repn) == VFloatRepn || (repn) == VDoubleRepn)


/*
 *  Dictionary describing keyword/value mapping.
 */

/* Dictionary entry: */
typedef struct {

    /* The following are initialized by the dictionary provider: */
    VStringConst keyword;		/* keyword string */
    VLong ivalue;			/* value, if an integer */
    VStringConst svalue;		/* value, if a string */

    /* The following are used only by code in VLookupDictValue: */
    VBoolean icached;			/* whether integer value cached */
    VBoolean fcached;			/* whether float value cached */
    VDouble fvalue;			/* cached floating-point value */
} VDictEntry;

/* Dictionaries of keywords for attribute values: */
extern VDictEntry VBooleanDict[];	/* boolean values */
extern VDictEntry VNumericRepnDict[];	/* numeric representation kinds */


/*
 *  Other definitions.
 */

/* Order in which to pack bits or bytes: */
typedef enum { VLsbFirst, VMsbFirst } VPackOrder;


/*
 *  Declarations of library routines.
 */

/* From Alloc.c: */

extern VPointer VCalloc (
#if NeedFunctionPrototypes
    size_t		/* n */,
    size_t		/* size */
#endif
);

extern void VFree (
#if NeedFunctionPrototypes
    VPointer		/* p */
#endif
);

extern VPointer VMalloc (
#if NeedFunctionPrototypes
    size_t		/* size */
#endif
);

extern VPointer VRealloc (
#if NeedFunctionPrototypes
    VPointer		/* p */,
    size_t		/* size */
#endif
);

/* From Attr.c: */

extern void VAppendAttr (
#if NeedVarargsPrototypes
    VAttrList		/* list */,
    VStringConst	/* name */,
    VDictEntry *	/* dict */,
    VRepnKind		/* repn */,
    ...			/* value */
#endif
);

extern VAttrList VCopyAttrList (
#if NeedFunctionPrototypes
    VAttrList		/* list */
#endif
);

extern VAttrList VCreateAttrList (
#if NeedFunctionPrototypes
    void
#endif
);

extern VBundle VCreateBundle (
#if NeedFunctionPrototypes
    VStringConst	/* type_name */,
    VAttrList		/* list */,
    size_t		/* length */,
    VPointer		/* data */
#endif
);

extern VBoolean VDecodeAttrValue (
#if NeedFunctionPrototypes
    VStringConst	/* str */,
    VDictEntry *	/* dict */,
    VRepnKind		/* repn */,
    VPointer		/* value */
#endif
);

extern void VDeleteAttr (
#if NeedFunctionPrototypes
    VAttrListPosn *	/* posn */
#endif
);

extern void VDestroyAttrList (
#if NeedFunctionPrototypes
    VAttrList		/* list */
#endif
);

extern void VDestroyBundle (
#if NeedFunctionPrototypes
    VBundle		/* bundle */
#endif
);

extern VStringConst VEncodeAttrValue (
#if NeedVarargsPrototypes
    VDictEntry *	/* dict */,
    VRepnKind		/* repn */,
    ...			/* value */
#endif
);

extern VBoolean VExtractAttr (
#if NeedFunctionPrototypes
    VAttrList		/* list */,
    VStringConst	/* name */,
    VDictEntry *	/* dict */,
    VRepnKind		/* repn */,
    VPointer		/* value */,
    VBooleanPromoted	/* required */
#endif
);

extern VGetAttrResult VGetAttr (
#if NeedFunctionPrototypes
    VAttrList		/* list */,
    VStringConst	/* name */,
    VDictEntry *	/* dict */,
    VRepnKind		/* repn */,
    VPointer		/* value */
#endif
);

extern VBoolean VGetAttrValue (
#if NeedFunctionPrototypes
    VAttrListPosn *	/* posn */,
    VDictEntry *	/* dict */,
    VRepnKind		/* repn */,
    VPointer		/* value */
#endif
);

extern void VInsertAttr (
#if NeedVarargsPrototypes
    VAttrListPosn *	/* posn */,
    VBooleanPromoted	/* after */,
    VStringConst	/* name */,
    VDictEntry *	/* dict */,
    VRepnKind		/* repn */,
    ...			/* value */
#endif
);

extern VBoolean VLookupAttr (
#if NeedFunctionPrototypes
    VAttrList		/* list */,
    VStringConst	/* name */,
    VAttrListPosn *	/* posn */
#endif
);

extern void VPrependAttr (
#if NeedVarargsPrototypes
    VAttrList		/* list */,
    VStringConst	/* name */,
    VDictEntry *	/* dict */,
    VRepnKind		/* repn */,
    ...			/* value */
#endif
);

extern void VSetAttr (
#if NeedVarargsPrototypes
    VAttrList		/* list */,
    VStringConst	/* name */,
    VDictEntry *	/* dict */,
    VRepnKind		/* repn */,
    ...			/* value */
#endif
);

extern void VSetAttrValue (
#if NeedVarargsPrototypes
    VAttrListPosn *	/* posn */,
    VDictEntry *	/* dict */,
    VRepnKind		/* repn */,
    ...			/* value */
#endif
);


/* new: geometry info */

extern void VCopyGeoInfo (
    VAttrList 	         /* src list */,
    VAttrList 	         /* dst list */
);

extern void VSetGeoInfo (
    VAttrList 	         /* geo list */,
    VAttrList 	         /* dst list */
);

extern VAttrList VGetGeoInfo(VAttrList);

extern double *VGetGeoPixdim(VAttrList,double *);
extern double *VGetGeoDim(VAttrList,double *);
extern double *VGetGeoSlicetime(VAttrList,double *,int);

extern void VSetGeoPixdim(VAttrList,double *);
extern void VSetGeoDim(VAttrList,double *);
extern void VGetSlicetimes(int,int,int,double,double *);


/* new: file IO, May 2017, updated Feb 2018 */
extern int CheckGzip(char *filename);
extern VAttrList VReadAttrList(VString filename,VLong tr,VBoolean attrtype,VBoolean do_scaling);
extern VAttrList VReadAttrListZ(FILE *,VString,VLong,VBoolean,VBoolean);

extern VImage VReadImage(VAttrList list);
extern int VAttrListNumImages(VAttrList list);
extern VImage *VAttrListGetImages(VAttrList list,int n);
extern void VImageDimensions(VImage *src,int nimages,int *bands,int *rows,int *cols);
extern void VDimensions(VImage *src,int nimages,int *nslices,int *nrows,int *ncols,int *nvolumes);
extern double VGetVoxel(VImage *src,int nimages,int b,int r,int c,int j);
extern int VSetVoxel(VImage *src,int nvolumes,int b,int r,int c,int j,double value);

  
/* automatic minval, July 2017, update May 2019 */
extern int Foreground(double *data,size_t n,int norm);
extern void VMaskMinval(VAttrList list,VImage mask,double minval);
extern void VMinval(VAttrList list,VString mask_filename,double minval);
extern void VMultMinval(VAttrList *list,int nlists,VString mask_filename,double minval);
extern int VSetVoxel(VImage *src,int nvolumes,int b,int r,int c,int j,double value);
extern double VGetVoxel(VImage *src,int nimages,int b,int r,int c,int j);
extern VImage VReadImageFile(VString filename);
  
/* From Dictionary.c: */

extern VDictEntry *VLookupDictKeyword (
#if NeedFunctionPrototypes
    VDictEntry *	/* dict */,
    VStringConst	/* keyword */
#endif
);

extern VDictEntry *VLookupDictValue (
#if NeedVarargsPrototypes
    VDictEntry *	/* dict */,
    VRepnKind		/* repn */,
    ...			/* value */
#endif
);

/* From Error.c: */

typedef void VErrorHandler (
#if NeedFunctionPrototypes
    VStringConst	/* message */
#endif
);

typedef void VWarningHandler (
#if NeedFunctionPrototypes
    VStringConst	/* message */
#endif
);

extern void VSetErrorHandler (
#if NeedFunctionPrototypes
    VErrorHandler *	/* fnc */
#endif
);

extern void VSetWarningHandler (
#if NeedFunctionPrototypes
    VWarningHandler *	/* fnc */
#endif
);

extern void VSetProgramName (
#if NeedFunctionPrototypes
    VStringConst	/* name */
#endif
);

extern void VError (
#if NeedVarargsPrototypes
    VStringConst	/* format */,
    ...			/* args */
#endif
);

extern void VWarning (
#if NeedVarargsPrototypes
    VStringConst	/* format */,
    ...			/* args */
#endif
);

extern void VSystemError (
#if NeedVarargsPrototypes
    VStringConst	/* format */,
    ...			/* args */
#endif
);

extern void VSystemWarning (
#if NeedVarargsPrototypes
    VStringConst	/* format */,
    ...			/* args */
#endif
);

extern void VDefaultError (
#if NeedFunctionPrototypes
    VStringConst	/* message */
#endif
);

extern void VDefaultWarning (
#if NeedFunctionPrototypes
    VStringConst	/* message */
#endif
);

/* From PackData.c: */

VBoolean VPackData (
#if NeedFunctionPrototypes
    VRepnKind		/* repn */,
    size_t		/* nels */,
    VPointer		/* unpacked */,
    VPackOrder		/* packed_order */,
    size_t *		/* length */,
    VPointer *		/* packed */,
    VBoolean *		/* alloced */
#endif
);

VBoolean VUnpackData (
#if NeedFunctionPrototypes
    VRepnKind		/* repn */,
    size_t		/* nels */,
    VPointer		/* packed */,
    VPackOrder		/* packed_order */,
    size_t *		/* length */,
    VPointer *		/* unpacked */,
    VBoolean *		/* alloced */
#endif
);

void VPackBits (
#if NeedFunctionPrototypes
    size_t		/* nels */,
    VPackOrder		/* packed_order */,
    VBit *		/* unpacked */,
    char *		/* packed */
#endif
);

void VUnpackBits (
#if NeedFunctionPrototypes
    size_t		/* nels */,
    VPackOrder		/* packed_order */,
    char *		/* packed */,
    VBit *		/* unpacked */
#endif
);

/* From Type.c: */

VRepnKind VRegisterType (
#if NeedFunctionPrototypes
    VStringConst	/* name */,
    VTypeMethods *	/* methods */
#endif
);

VRepnKind VLookupType (
#if NeedFunctionPrototypes
    VStringConst	/* name */
#endif
);

#ifdef __cplusplus
}
#endif

#endif /* V_Vlib_h */
