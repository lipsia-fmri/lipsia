/*
 *  $Id: ReadPnm.c 3177 2008-04-01 14:47:24Z karstenm $
 *
 *  This file contains routines for reading portable anymap (PNM) files.
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
VRcsId ("$Id: ReadPnm.c 3177 2008-04-01 14:47:24Z karstenm $");

/* Buffer larger enough to contain any input token: */
static char token[32];
static int lineno;

/* Types of file format. ParsePnmHeader() depends on the order of these: */
typedef enum { Pbm, Pgm, Ppm, PbmRaw, PgmRaw, PpmRaw } PnmFormat;

/* Later in this file: */
static VBoolean ParsePnmHeader (FILE *, PnmFormat *, long *, long *, long *);
static VImage ReadPbm (FILE *, VBoolean, long, long);
static VImage ReadPgm (FILE *, VBoolean, long, long, long);
static VImage ReadPpm (FILE *, VBoolean, long, long, long);
static int ReadByte (FILE *, VBoolean);
static VBoolean NextToken (FILE *);
static char NextDataChar (FILE *);


/* 
 *  VReadPnm
 * 
 *  Read an image from a Portable Bitmap, Graymap, or Pixmap (PBM, PGM, or PPM)
 *  file and return it as a VImage.
 *  
 */

VImage VReadPnm (FILE *f)
{
    PnmFormat format;
    long width, height, maxc;
    VImage image;

    lineno = 1;

    /* Parse PNM file header: */
    if (ParsePnmHeader (f, & format, & width, & height, & maxc)) {
    
	/* Perform conversion: */
	switch (format) {
	case Pbm:    image = ReadPbm (f, FALSE, width, height);         break;
	case PbmRaw: image = ReadPbm (f, TRUE, width, height);       	break;
	case Pgm:    image = ReadPgm (f, FALSE, width, height, maxc);   break;
	case PgmRaw: image = ReadPgm (f, TRUE, width, height, maxc); 	break;
	case Ppm:    image = ReadPpm (f, FALSE, width, height, maxc);   break;
	case PpmRaw: image = ReadPpm (f, TRUE, width, height, maxc); 	break;
	}
    } else image = NULL;

    if (! image)
	VWarning ("VReadPnm: Error in PNM file near line %d", lineno);
    else if (NextToken (f)) {
	VWarning ("VReadPnm: File continues beyond expected EOF");
	VDestroyImage (image);
	image = NULL;
    }
    return image;
}


/*
 *  ParsePnmHeader
 *
 *  Parse the header of a PNM file.
 */

static VBoolean ParsePnmHeader (FILE *f, PnmFormat *format,
				long *width, long *height, long *maxc)
{
    char *cp;

    /* Read magic number to determine PNM format: */
    if (! NextToken (f))
	return FALSE;
    if (strlen (token) != 2 ||
	token[0] != 'P' || token[1] < '1' || token[1] > '6') {
	VWarning ("VReadPnm: File format \"%s\" not recognized", token);
	return FALSE;
    }
    *format = token[1] - '1';

    /* Read width: */
    if (! NextToken (f)) {
	VWarning ("VReadPnm: Unexpected EOF");
	return FALSE;
    }
    *width = strtol (token, & cp, 10);
    if (*cp || *width < 1) {
	VWarning ("VReadPnm: Bad width: %s", token);
	return FALSE;
    }

    /* Read height: */
    if (! NextToken (f)) {
	VWarning ("VReadPnm: Unexpected EOF");
	return FALSE;
    }
    *height = strtol (token, & cp, 10);
    if (*cp || *height < 1) {
	VWarning ("VReadPnm: Bad height: %s", token);
	return FALSE;
    }

    /* Read maxc (max color value) for PGM and PPM: */
    if (*format != Pbm && *format != PbmRaw) {
	if (! NextToken (f)) {
	    VWarning ("VReadPnm: Unexpected EOF");
	    return FALSE;
	}
	*maxc = strtol (token, & cp, 10);
	if (*cp || *maxc < 1) {
	    VWarning ("VReadPnm: Bad max color value: %s", token);
	    return FALSE;
	}
    }

    return TRUE;
}


/*
 *  ReadPbm
 *
 *  Read a PBM format file.
 */

static VImage ReadPbm (FILE *f, VBoolean raw, long width, long height)
{
    VImage image;
    VBit *pp;
    int i, j, m, bit;
    char ch;

    if (! (image = VCreateImage (1, height, width, VBitRepn)))
	return NULL;
    pp = VImageData (image);
    if (raw)
	for (i = 0; i < height; i++) {
	    bit = -1;
	    for (j = 0; j < width; j++, bit--) {
		if (bit < 0) {
		    bit = 7;
		    if ((m = getc (f)) == EOF)
			goto Eof;
		}
		*pp++ = ! ((m >> bit) & 1);
	    }
	}
    else
	for (i = width * height; i > 0; i--) {
	    if ((ch = NextDataChar (f)) == EOF)
		goto Eof;
	    *pp++ = (ch == '0');
	}

    return image;

Eof:
    VWarning ("VReadPnm: Unexpected EOF");
    VDestroyImage (image);
    return NULL;
}


/*
 *  ReadPgm
 *
 *  Read a PGM format file.
 */

static VImage ReadPgm (FILE *f, VBoolean raw, long width, long height,
		       long maxc)
{
    VImage image;
    VUByte *pp;
    int i, m;
    double scaling_factor;

    if (! (image = VCreateImage (1, height, width, VUByteRepn)))
	return NULL;
    pp = VImageData (image);
    scaling_factor = 255.0 / maxc;
    for (i = width * height; i > 0; i--) {
	if ((m = ReadByte (f, raw)) == EOF) {
	    VDestroyImage (image);
	    return NULL;
	}
	*pp++ = m * scaling_factor;
    }
    return image;
}


/*
 *  ReadPpm
 *
 *  Read a PPM format file.
 */

static VImage ReadPpm (FILE *f, VBoolean raw, long width, long height,
		       long maxc)
{
    VImage image;
    VUByte *pp[3];
    int i, j, m;
    double scaling_factor;
 
    if (! (image = VCreateImage (3, height, width, VUByteRepn)))
	return NULL;    
    VImageNFrames (image) = 1;
    VImageNColors (image) = 3;
    VAppendAttr (VImageAttrList (image), VColorInterpAttr, VBandInterpDict,
		 VLongRepn, (VLong) VBandInterpRGB);
    for (j = 0; j < 3; j++)
	pp[j] = VPixelPtr (image, j, 0, 0);
    /* scaling_factor = 255.0 / maxc; */
    for (i = width * height; i > 0; i--)
	for (j = 0; j < 3; j++) {
	    if ((m = ReadByte (f, raw)) == EOF) {
		VDestroyImage (image);
		return NULL;
	    }
	    *pp[j]++ = m;
	}
    return image;
}


/*
 *  ReadByte
 *
 *  Return the next byte value from the file, which may be represented
 *  as a decimal integer or as a single, binary-encoded byte.
 */

static int ReadByte (FILE *f, VBoolean raw)
{
    int ch;
    char *cp;

    if (raw) {
	if ((ch = getc (f)) != EOF)
	    return ch;
    } else {
	if (NextToken (f)) {
	    ch = strtol (token, & cp, 10);
	    return *cp ? -1 : ch;
	}
    }
    VWarning ("VReadPnm: Unexpected EOF");
    return EOF;
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
    int ch = NextDataChar (f);

    if (ch == EOF)
	return FALSE;

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


/*
 *  NextDataChar
 *
 *  Reads the next non-comment, non-whitespace character. Also tracks
 *  the current line number in the input file.
 */

static char NextDataChar (FILE *f)
{
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
		    return EOF;
	    } while (ch != '\n');
	    lineno++;
	    ch = getc (f);
	    if (ch == '\n')
		lineno++;
	}
	
	/* Return "ch" if it's EOF or not whitespace: */
	if (ch == EOF || ! isspace (ch))
	    return ch;
    }
}
