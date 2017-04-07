/*
 *  $Id: option.h 3001 2007-11-30 12:42:11Z karstenm $
 *
 *  Definitions associated with command option and parameter parsing.
 */

#ifndef V_option_h
#define V_option_h 1

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

/* From the standard C library: */
#include <stdio.h>


#ifdef __cplusplus
extern "C" {
#endif


/*
 *  Declarations of data structures.
 */

/* Accepted command options are described by a table of these entries: */
typedef struct {
    VStringConst keyword;		/* keyword signalling option */
    VRepnKind repn;			/* type of value supplied by option */
    int number;				/* number of values supplied */
    VPointer value;			/* location for storing value(s) */
    VBoolean *found;			/* whether optional arg found */
    VDictEntry *dict;			/* optional dict of value keywords */
    VStringConst blurb;			/* on-line help blurb */
} VOptionDescRec;

/* Some possible values for the found field: */
#define VRequiredOpt (& V_RequiredOpt)
#define VOptionalOpt (& V_OptionalOpt)
extern VBoolean V_RequiredOpt;
extern VBoolean V_OptionalOpt;

/* If an option takes multiple values, they are represented by a VArgVector: */
typedef struct {
    int number;				/* number of arguments */
    VPointer vector;			/* vector of arguments */
} VArgVector;


/*
 *  Declarations of library routines.
 */

/* history management: */
void VHistory (
#if NeedFunctionPrototypes
    int,
    VOptionDescRec*,
    char*,
    VAttrList*,
    VAttrList*
#endif
);

VAttrList VReadHistory(
#if NeedFunctionPrototypes
    VAttrList*
#endif
);

void VPrependHistory(
#if NeedFunctionPrototypes
int,
VOptionDescRec *,
char*,
VAttrList*
#endif
);

char * getVersion();


/* From Option.c: */

VBoolean VIdentifyFiles (
#if NeedFunctionPrototypes
    int			/* noptions */,
    VOptionDescRec []	/* options */,
    VStringConst	/* keyword */,
    int *		/* argc */,
    char **		/* argv */,
    int			/* fd */
#endif
);

VBoolean VParseCommand (
#if NeedFunctionPrototypes
    int			/* noptions */,
    VOptionDescRec []	/* options */,
    int *		/* argc */,
    char **		/* argv */
#endif
);

VBoolean VParseCommand_nl (
#if NeedFunctionPrototypes
    int			/* noptions */,
    VOptionDescRec []	/* options */,
    int *		/* argc */,
    char **		/* argv */
#endif
);

void VParseFilterCmd (
#if NeedFunctionPrototypes
    int			/* noptions */,
    VOptionDescRec []	/* opts */,
    int			/* argc */,
    char **		/* argv */,
    FILE **		/* inp */,
    FILE **		/* outp */
#endif
);

void VPrintOptions (
#if NeedFunctionPrototypes
    FILE *		/* f */,
    int			/* noptions */,
    VOptionDescRec []	/* options */
#endif
);

int VPrintOptionValue (
#if NeedFunctionPrototypes
    FILE *		/* f */,
    VOptionDescRec *	/* option */
#endif
);

void VReportBadArgs (
#if NeedFunctionPrototypes
    int			/* argc */,
    char **		/* argv */
#endif
);

void VReportUsage (
#if NeedFunctionPrototypes
    VStringConst	/* program */,
    int			/* noptions */,
    VOptionDescRec []	/* options */,
    VStringConst	/* other_args */
#endif
);

void VReportUsage_nl (
#if NeedFunctionPrototypes
    VStringConst	/* program */,
    int			/* noptions */,
    VOptionDescRec []	/* options */,
    VStringConst	/* other_args */
#endif
);

void VReportValidOptions (
#if NeedFunctionPrototypes
    int			/* noptions */,
    VOptionDescRec []	/* options */
#endif
);

/* From Param.c: */

VBoolean VLoadParameters (
#if NeedFunctionPrototypes
    int			/* nparams */,
    VOptionDescRec []	/* params */,
    VStringConst	/* filename */,
    VStringConst	/* object */,
    VPointer		/* values */,
    VBoolean		/* ignore_unrecog */
#endif
);

VBoolean VParseParamDefn (
#if NeedFunctionPrototypes
    VStringConst	/* buf */,
    VString		/* name */,
    VRepnKind *		/* repn */,
    VString		/* value */,
    VString		/* comment */
#endif
);

VBoolean VParseParamOptions (
#if NeedFunctionPrototypes
    int			/* nparams */,
    VOptionDescRec []	/* params */,
    int *		/* argc */,
    char **		/* argv */,
    VPointer		/* values */
#endif
);

void VPrintParameters (
#if NeedFunctionPrototypes
    FILE *		/* f */,
    int			/* nparams */,
    VOptionDescRec []	/* params */,
    VPointer		/* values */
#endif
);

void VReportValidParamOptions (
#if NeedFunctionPrototypes
    int			/* nparams */,
    VOptionDescRec []	/* params */,
    VPointer		/* values */
#endif
);

#ifdef __cplusplus
}
#endif

#endif /* V_option_h */
