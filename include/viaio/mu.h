/*
 *  $Id: mu.h 726 2004-03-08 13:12:45Z lohmann $
 *
 *  Miscellaneous utility macros provided by Vista.
 */

#ifndef V_mu_h
#define V_mu_h 1

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
#include <string.h>

/* Under SunOS 4.x, definition of memxxx() aren't in <string.h>. */
#ifdef sun
#if NeedFunctionPrototypes
extern void *memcpy (void *, const void *, size_t);
extern void *memset (void *, int, size_t);
#else
extern void *memcpy ();
extern void *memset ();
#endif
#endif


/*
 *  Convenience macros.
 */

/* The minimum and maximum of two values: */
/* #if __GNUC__ && ! __STRICT_ANSI__ */
/* #define VMin(a,b)			\ */
/*     ( { typedef _ta = (a), _tb = (b);	\ */
/* 	_ta _a = (a); _tb _b = (b);	\ */
/* 	_a < _b ? _a : _b; } ) */
/* #define VMax(a,b)			\ */
/*     ( { typedef _ta = (a), _tb = (b);	\ */
/* 	_ta _a = (a); _tb _b = (b);	\ */
/* 	_a > _b ? _a : _b; } ) */
/* #else */
#define VMax(a,b)	((a) > (b) ? (a) : (b))
#define VMin(a,b)	((a) < (b) ? (a) : (b))
/* #endif */

/* The offset of a field within a particular structure type: */
#define VOffset(type, field) \
    ((size_t) (((char *) & ((type) 0)->field) - (char *) 0))

#define VOffsetOf(type, field)	VOffset(type *, field)

/* The length of a statically-allocated one-dimensional array: */
#define VNumber(array)	((size_t) (sizeof (array) / sizeof ((array)[0])))

/* Zero out a one-dimensional array: */
#define VZero(array, nels) \
    ((void) memset ((void *) array, 0, (size_t) (nels) * sizeof ((array)[0])))

/* Copy one vector to another: */
#define VCopy(from, to, nels) \
    ((void) memcpy ((void *) (to), (void *) (from), \
		    (size_t) (nels) * sizeof ((from)[0])))

/* Allocate storage for a particular type of object: */
#define VNew(type)	((type *) VMalloc (sizeof (type)))

/* Copy a string into a new block of storage: */
#define VNewString(str) \
    ((VString) ((str) ? strcpy ((char *) VMalloc (strlen (str) + 1), str) : 0))

#endif /* V_mu_h */
