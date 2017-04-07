/*
 *  $Id: file.h 726 2004-03-08 13:12:45Z lohmann $
 *
 *  This file contains definitions related to the Vista data file format.
 */

#ifndef V_file_h
#define V_file_h 1

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
#include <viaio/Vlib.h>

/* From the standard C library: */
#ifdef NULL	/* because some stdio's blindly defined NULL */
#undef NULL
#endif
#include <stdio.h>
#ifndef NULL
#define NULL 0
#endif


#ifdef __cplusplus
extern "C" {
#endif

/*
 *  Components of the Vista data file header, and the delimeter separating
 *  attributes from binary data.
 */

#define VFileHeader		"V-data"
#define VFileVersion		2
#define VFileDelimiter		"\f\n"


/*
 *  Limit on the length of an attribute name.
 */

#define VMaxAttrNameLength	256


/*
 *  Type of function supplied as a filter to VReadFile.
 */

typedef VBoolean VReadFileFilterProc (
#if NeedFunctionPrototypes
    VBundle	    	    /* bundle */,
    VRepnKind		    /* repn */
#endif
);


/*
 *  Declarations of library routines.

 */

/* Open a file for input: */
extern FILE *VOpenInputFile (
#if NeedFunctionPrototypes
    VStringConst	/* filename */,
    VBoolean		/* nofail */
#endif
);

/* Open a file for output: */
extern FILE *VOpenOutputFile (
#if NeedFunctionPrototypes
    VStringConst	/* filename */,
    VBoolean		/* nofail */
#endif
);

/* Read objects of a certain type: */
extern int VReadObjects (
#if NeedFunctionPrototypes
    FILE *		/* file */,
    VRepnKind		/* repn */,
    VAttrList *		/* attributes */,
    VPointer **		/* objects */
#endif
);

/* Read a Vista data file: */
extern VAttrList VReadFile (
#if NeedFunctionPrototypes
    FILE *		/* f */,
    VReadFileFilterProc	* /* filter */
#endif
);

/* Write objects of a certain type: */
extern VBoolean VWriteObjects (
#if NeedFunctionPrototypes
    FILE *		/* file */,
    VRepnKind		/* repn */,
    VAttrList		/* attributes */,
    int			/* nobjects */,
    VPointer []		/* objects */
#endif
);

/* Write a Vista data file: */
extern VBoolean VWriteFile (
#if NeedFunctionPrototypes
    FILE *		/* f */,
    VAttrList		/* list */
#endif
);

#ifdef __cplusplus
}
#endif

#endif /* V_file_h */
