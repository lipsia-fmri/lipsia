/*
 *  $Id: Alloc.c 3177 2008-04-01 14:47:24Z karstenm $
 *
 *  This file contains routines for allocating and freeing memory, with
 *  error checking.
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

/*
 *  VMalloc
 *
 *  Perform error checking on malloc() call.
 */

VPointer VMalloc (size_t size)
{
    VPointer p;

    if (size == 0)
	return NULL;
    if (! (p = (VPointer) malloc (size)))
	VSystemError ("VMalloc: Memory allocation failure");
    return p;
}


/*
 *  VRealloc
 *
 *  Perform error checking on realloc() call.
 */

VPointer VRealloc (VPointer p, size_t size)
{
    if (size == 0) {
	VFree (p);
	return NULL;
    }
    if (! p)
	return VMalloc (size);
    if (! (p = (VPointer) realloc (p, size)))
	VSystemError ("VRealloc: Memory allocation failure");
    return p;
}


/*
 *  VCalloc
 *
 *  Perform error checking on calloc() call.
 */

VPointer VCalloc (size_t n, size_t size)
{
    VPointer p;

    if (n == 0 || size == 0)
	return NULL;
    if (! (p = (VPointer) calloc (n, size)))
	VSystemError ("VCalloc: Memory allocation failure");
    return p;
}


/*
 *  VFree
 *
 *  Perform error checking on free() call.
 */

void VFree (VPointer p)
{
    if (p)
	free ((char *) p);
}
