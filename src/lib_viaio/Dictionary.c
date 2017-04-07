/*
 *  $Id: Dictionary.c 3177 2008-04-01 14:47:24Z karstenm $
 *
 *  This file contains dictionaries that symbols used in V data files to
 *  values used internally.
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
 *  Dictionaries for generic attributes and their values.
 */

/* Keywords for representing TRUE or FALSE: */
VDictEntry VBooleanDict[] = {
    { "false",		FALSE },
    { "true",		TRUE },
    { "no",		FALSE },
    { "yes",		TRUE },
    { "off",		FALSE },
    { "on",		TRUE },
    { NULL }
};

/* Keywords for representing kinds of numeric representation: */
VDictEntry VNumericRepnDict[] = {
    { "bit",		VBitRepn },
    { "double",		VDoubleRepn },
    { "float",		VFloatRepn },
    { "long",		VLongRepn },
    { "ulong",		VULongRepn },
    { "sbyte",		VSByteRepn },
    { "short",		VShortRepn },
    { "ushort",		VUShortRepn },
    { "int",		VIntegerRepn },
    { "uint",		VUIntegerRepn },
    { "ubyte",		VUByteRepn },
    { NULL }
};


/*
 *  VLookupDictKeyword
 *
 *  Look up an entry in an attribute value dictionary, by keyword.
 *  (It's assumed that dictionaries are pretty small -- a linear search
 *  is done.)
 */

VDictEntry *VLookupDictKeyword (VDictEntry *dict, VStringConst keyword)
{
    if (dict)
	for ( ; dict->keyword; dict++)
	    if (strcmp (keyword, dict->keyword) == 0)
		return dict;
    return NULL;
}


/*
 *  VLookupDictValue
 *
 *  Look up an entry in an attribute dictionary, by value.
 *  Calling sequence:
 *
 *	VLookupDictValue (VDictEntry *dict, VRepnKind repn, xxx value)
 *
 *  where xxx is a type that corresponds to repn.
 */

VDictEntry *VLookupDictValue (VDictEntry *dict, VRepnKind repn, ...)
{
    va_list args;
    VLong i_value = 0;
    VDouble f_value = 0.0;
    VString s_value = NULL;
    VBoolean i_valid;

    /* Unravel the arguments passed: */
    if (! dict)
	return NULL;
    va_start (args, repn);
    switch (repn) {
    case VBitRepn: i_value = va_arg (args, VBitPromoted); break;
    case VUByteRepn: i_value = va_arg (args, VUBytePromoted); break;
    case VSByteRepn: i_value = va_arg (args, VSBytePromoted); break;
    case VShortRepn: i_value = va_arg (args, VShortPromoted); break;
    case VUShortRepn: i_value = va_arg (args, VUShortPromoted); break;
    case VIntegerRepn: i_value = va_arg (args, VIntegerPromoted); break;
    case VUIntegerRepn: i_value = va_arg (args, VUIntegerPromoted); break;
    case VLongRepn: i_value = va_arg (args, VLongPromoted); break;
    case VULongRepn: i_value = va_arg (args, VULongPromoted); break;
    case VFloatRepn: f_value = va_arg (args, VFloatPromoted); break;
    case VDoubleRepn: f_value = va_arg (args, VDoublePromoted); break;
    case VBooleanRepn: i_value = va_arg (args, VBooleanPromoted); break;
    case VStringRepn: s_value = va_arg (args, VString); break;
    default:
	VError ("VLookupDictValue: Can't lookup %s value", VRepnName (repn));
    }
    va_end (args);

    /* Search the dictionary by value: */
    switch (repn) {

    case VBitRepn:
    case VUByteRepn:
    case VSByteRepn:
    case VShortRepn:
    case VUShortRepn:
    case VIntegerRepn:
    case VUIntegerRepn:
    case VLongRepn:
    case VULongRepn:
    case VBooleanRepn:
	for ( ; dict->keyword; dict++) {

	    /* Is the entry's value only stored as a string? */
	    if (dict->svalue && ! dict->icached) {

		/* Yes -- try to convert the string to an integer, and
		   cache that value: */
		if (! VDecodeAttrValue (dict->svalue, NULL,
					VLongRepn, & dict->ivalue))
		    break;
		dict->icached = TRUE;
	    }

	    /* Test against the integer value stored in the entry: */
	    if (i_value == dict->ivalue)
		return dict;
	}
	break;

    case VFloatRepn:
    case VDoubleRepn:
	for ( ; dict->keyword; dict++) {

	    /* Does the entry include a cached floating point value? */
	    if (! dict->fcached) {

		/* No -- obtain it from an integer or string value: */
		if (dict->svalue) {
		    if (! VDecodeAttrValue (dict->svalue, NULL,
					    VDoubleRepn, & dict->fvalue))
			break;
		} else dict->fvalue = dict->ivalue;
		dict->fcached = TRUE;
	    }

	    /* Test against the cached float value now stored in the entry: */
	    if (f_value == dict->fvalue)
		return dict;
	}
	break;

    case VStringRepn:

	/* In case we're searching a dictionary with only integer values
	   stored, try to convert the supplied string value to an integer: */
	i_valid = VDecodeAttrValue (s_value, NULL, VLongRepn, & i_value);

	for ( ; dict->keyword; dict++) {

	    /* If the entry includes a string value, compare with it: */
	    if (dict->svalue) {
		if (strcmp (s_value, dict->svalue) == 0)
		    return dict;
	    }

	    /* Otherwise, compare with its integer value: */
	    else if (i_valid && i_value == dict->ivalue)
		return dict;
	}
	break;

    default:
	break;
    }
    return NULL;
}
