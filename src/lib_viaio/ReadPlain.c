/*
 *  $Id: ReadPlain.c 3177 2008-04-01 14:47:24Z karstenm $
 *
 *  This file contains routines for reading Vista plain format files.
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
#include "viaio/file.h"
#include "viaio/os.h"
#include "viaio/VImage.h"

/* From the standard C library: */
#include <ctype.h>

/* File identification string: */
VRcsId ("$Id: ReadPlain.c 3177 2008-04-01 14:47:24Z karstenm $");

/* Later in this file: */
static VBoolean ParseHeader (FILE *, VRepnKind *, long *, long *, long *);
static VBoolean NextToken (FILE *);

/* Buffer larger enough to contain any input token: */
static char token[32];
static int lineno;


/* 
 *  VReadPlain
 * 
 *  Read an image from a Vista plain format file.
 */

VImage VReadPlain (FILE *f)
{
    VRepnKind repn;
    long ncolumns, nrows, nbands;
    VImage work = NULL, result;
    VDouble *pp;
    int i;
    char *cp;
    
    lineno = 1;

    /* Parse file header: */
    if (! ParseHeader (f, & repn, & nbands, & nrows, & ncolumns))
	goto Error;
    
    /* Fill a Double image with pixel values from the file: */
    work = VCreateImage (nbands, nrows, ncolumns, VDoubleRepn);
    if (! work)
	return NULL;
    for (i = VImageNPixels (work), pp = VImageData (work); i > 0; i--) {
	if (! NextToken (f)) {
	    VWarning ("VReadPlain: Unexpected EOF");
	    return FALSE;
	}
	*pp++ = strtod (token, & cp);
   if (*cp) {
	    VWarning ("VReadPlain: Bad pixel value: %s", token);
	    goto Error;
	}
    }

    /* There should be no more tokens: */
    if (NextToken (f)) {
	VWarning ("VReadPlain: File continues beyond expected EOF");
	goto Error;
    }

    /* Convert the image to the representation specified in the file: */
    if (repn == VDoubleRepn)
	return work;
    result = VConvertImageCopy (work, NULL, VAllBands, repn);
    VDestroyImage (work);
    return result;

Error:
    VWarning ("VReadPlain: Error in Vista plain format file near line %d",
	      lineno);
    if (work)
	VDestroyImage (work);
    return NULL;
}


/*
 *  ParseHeader
 *
 *  Parse the header of a Vista plain format file.
 */

static VBoolean ParseHeader (FILE *f, VRepnKind *repn,
			     long *nbands, long *nrows, long *ncolumns)
{
    char *cp;

    /* Read repn: */
    if (! NextToken (f))
	return FALSE;
    switch (*repn = VLookupType (token)) {

    case VBitRepn:
    case VUByteRepn:
    case VSByteRepn:
    case VShortRepn:
    case VLongRepn:
    case VFloatRepn:
    case VDoubleRepn:
	break;

    default:
	VWarning ("VReadPlain: Representation \"%s\" not recognized", token);
	return FALSE;
    }
    
    /* Read nbands: */
    if (! NextToken (f)) {
	VWarning ("VReadPlain: Unexpected EOF");
	return FALSE;
    }
    *nbands = strtol (token, & cp, 10);
    if (*cp || *nbands < 1) {
	VWarning ("VReadPlain: Bad number of bands: %s", token);
	return FALSE;
    }
    
    /* Read nrows: */
    if (! NextToken (f)) {
	VWarning ("VReadPlain: Unexpected EOF");
	return FALSE;
    }
    *nrows = strtol (token, & cp, 10);
    if (*cp || *nrows < 1) {
	VWarning ("VReadPlain: Bad number of rows: %s", token);
	return FALSE;
    }

    /* Read ncolumns: */
    if (! NextToken (f)) {
	VWarning ("VReadPlain: Unexpected EOF");
	return FALSE;
    }
    *ncolumns = strtol (token, & cp, 10);
    if (*cp || *ncolumns < 1) {
	VWarning ("VReadPlain: Bad number of columns: %s", token);
	return FALSE;
    }
    
    return TRUE;
}


/*
 *  NextToken
 *
 *  Copies the next input token into `token', returning TRUE if successful
 *  or FALSE if encountering an EOF. Also tracks the current line number
 *  in the input file.
 */

static VBoolean NextToken (FILE *f)
{
    char *cp = token;
    int ch;

    /* Ignore leading whitespace and comments: */
    while (1) {

	/* Ignore leading spaces: */
	do {
	    ch = getc (f);
	    if (ch == '\n')
		lineno++;
	} while (isspace (ch));
	
	/* Ignore leading comments: */
	while (ch == '#') {
	    do {
		ch = getc (f);
		if (ch == EOF)
		    return FALSE;
	    } while (ch != '\n');
	    lineno++;
	    ch = getc (f);
	    if (ch == '\n')
		lineno++;
	}
	
	/* Return FALSE if EOF encountered: */
	if (ch == EOF)
	    return FALSE;
	
	/* Exit loop if "ch" is start of token: */
	if (! isspace (ch))
	    break;
    }

    /* Read until whitespace is encountered: */
    while (ch != EOF && ! isspace (ch)) {
	if (cp - token < sizeof (token) - 1)
	    *cp++ = ch;
	ch = getc (f);
    }
    if (ch == '\n')
	lineno++;
    *cp = 0;
    return TRUE;
}
