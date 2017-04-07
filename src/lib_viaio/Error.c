/*
 *  $Id: Error.c 3177 2008-04-01 14:47:24Z karstenm $
 *
 *  This file contains routines for reporting errors.
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

/* From the standard C library: */
#include <errno.h>

/* File identification string: */
VRcsId ("$Id: Error.c 3177 2008-04-01 14:47:24Z karstenm $");

/* Name of the program, for use in error messages: */
static VStringConst programName;

/* Routines to call to report errors: */
static VErrorHandler *errorHandler = VDefaultError;
static VWarningHandler *warningHandler = VDefaultWarning;

/* Maximum length of an error message: */
#define maxErrorMessageLength	500

/* Later in this file: */
static void FormatMsg (char *, VStringConst, VStringConst, va_list *,
		       VStringConst);


/*
 *  VSetErrorHandler
 *
 *  Establish a caller-supplied routine as the handler for error messages.
 */

void VSetErrorHandler (VErrorHandler *fnc)
{
    errorHandler = fnc ? fnc : VDefaultError;
}


/*
 *  VSetWarningHandler
 *
 *  Establish a caller-supplied routine as the handler for warning messages.
 */

void VSetWarningHandler (VWarningHandler *fnc)
{
    warningHandler = fnc ? fnc : VDefaultWarning;
}


/*
 *  VSetProgramName
 *
 *  Note the program name for use in error messages.
 */

void VSetProgramName (VStringConst name)
{
    programName = name;
}


/*
 *  VError
 *
 *  Report a fatal program error.
 */

void VError (VStringConst format, ...)
{
    va_list args;
    char buf[maxErrorMessageLength];

    va_start (args, format);
    FormatMsg (buf, "Fatal", format, & args, 0);
    va_end (args);
    (*errorHandler) (buf);
    fprintf (stderr, "%s: Fatal: Unexpected return from error handler.\n",
	     programName);
    exit (EXIT_FAILURE);
}


/*
 *  VWarning
 * 
 *  Report a non-fatal program error.
 */

void VWarning (VStringConst format, ...)
{
    va_list args;
    char buf[maxErrorMessageLength];

    va_start (args, format);
    FormatMsg (buf, "Warning", format, & args, 0);
    va_end (args);
    (*warningHandler) (buf);
}


/*
 *  VSystemError
 *
 *  Report a fatal error incurred by a system call.
 */

void VSystemError (VStringConst format, ...)
{
    va_list args;
    char buf[maxErrorMessageLength];

    va_start (args, format);
    FormatMsg (buf, "Fatal", format, & args, strerror (errno));
    va_end (args);
    (*errorHandler) (buf);
    fprintf (stderr, "%s: Fatal: Unexpected return from error handler.\n",
	     programName);
    exit (EXIT_FAILURE);
}


/*
 *  VSystemWarning
 *
 *  Report non-fatal error incurred by a system call.
 */

void VSystemWarning (VStringConst format, ...)
{
    va_list args;
    char buf[maxErrorMessageLength];

    va_start (args, format);
    FormatMsg (buf, "Warning", format, & args, strerror (errno));
    va_end (args);
    (*warningHandler) (buf);
}


/*
 *  FormatMsg -- format an error message for output.
 */

static void FormatMsg (char *buf, VStringConst severity, VStringConst format,
		       va_list *args, VStringConst extra)
{
    if (programName)
	sprintf (buf, "%s: %s: ", programName, severity);
    else sprintf (buf, "%s: ", severity);
    vsprintf (buf + strlen (buf), format, *args);
    if (extra)
	sprintf (buf + strlen (buf), ": %s", extra);
    strcat (buf, ".\n");
}


/*
 *  VDefaultError
 *
 *  Default error handler.
 */

void VDefaultError (VStringConst msg)
{
    fputs (msg, stderr);
    exit (EXIT_FAILURE);
}


/*
 *  VDefaultWarning
 *
 *  Default warning handler.
 */

void VDefaultWarning (VStringConst msg)
{
    fputs (msg, stderr);
}
