/*
 *  $Id: Type.c 3177 2008-04-01 14:47:24Z karstenm $
 *
 *  Routines in this file allow one to register types (e.g., "image") along
 *  with sets of routines for operating on objects of those types.
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
#include "viaio/VList.h"
#include "viaio/mu.h"
#include <limits.h>

/*
 *  Information about built-in types.
 */

extern VTypeMethods VEdgesMethods;	/* in EdgesType.c */
extern VTypeMethods VImageMethods;	/* in ImageType.c */
extern VTypeMethods VolumesMethods;     /* in VolumesType.c */

static VRepnInfoRec builtin_repn_info[] = {
    { "unknown" },

    /* Integer and floating-point numbers: */
    { "bit",	   sizeof (VBit),         1, 0, 1 },
    { "ubyte",	   sizeof (VUByte),	  8, 0, UCHAR_MAX },
    { "sbyte",	   sizeof (VSByte),	  8, CHAR_MIN, CHAR_MAX },
    { "short",	   sizeof (VShort),      16, -32768.0, 32767.0 },
    { "long",	   sizeof (VLong),       64, LONG_MIN,LONG_MAX }, 
    { "float",	   sizeof (VFloat),      32, 
	  -3.40282346638528860e+38, 3.40282346638528860e+38 },
    { "double",	   sizeof (VDouble),   64,
	  -1.797693134862315708e+308, 1.797693134862315708e+308 },

    /* Miscellaneous representations: */
    { "attr-list", sizeof (VAttrList),  0, 0.0, 0.0 },
    { "boolean",   sizeof (VBoolean),	1, 0.0, 0.0 },
    { "bundle",	   sizeof (VPointer),	0, 0.0, 0.0 },
    { "list",	   sizeof (VList),	0, 0.0, 0.0 },
    { "pointer",   sizeof (VPointer),	0, 0.0, 0.0 },
    { "string",	   sizeof (VString),	0, 0.0, 0.0 },

    /* Standard object types: */
    { "edges",	   sizeof (VPointer),	0, 0.0, 0.0, & VEdgesMethods },
    { "image",	   sizeof (VPointer),	0, 0.0, 0.0, & VImageMethods },
    
    /* new object types */
    { "volumes",   sizeof (VPointer),   0, 0.0, 0.0, & VolumesMethods },
    { "ushort",	   sizeof (VUShort),     16, 0, USHRT_MAX },
    { "int",	   sizeof (VInteger),    32, INT_MIN, INT_MAX },
    { "uint",	   sizeof (VUInteger),   32, 0, UINT_MAX },
    { "ulong",	   sizeof (VULong),      64, 0, ULONG_MAX }
};

VRepnInfoRec *VRepnInfo = builtin_repn_info;
static VRepnKind nRepnKinds = VNRepnKinds;


/*
 *  VRegisterType
 *
 *  Register some handlers for dealing with objects of a particular type.
 *  Returns the VRepnKind code assigned the new type.
 */

VRepnKind VRegisterType (VStringConst name, VTypeMethods *methods)
{
    VRepnInfoRec *p;

    /* Move the existing type information into a bigger table: */
    if (VRepnInfo == builtin_repn_info) {
	VRepnInfo = VMalloc ((VNRepnKinds + 1) * sizeof (VRepnInfoRec));
	VCopy (builtin_repn_info, VRepnInfo, VNRepnKinds);
    } else
	VRepnInfo =
	    VRealloc (VRepnInfo, (nRepnKinds + 1) * sizeof (VRepnInfoRec));

    /* Write the new type's info into the last table entry: */
    p = VRepnInfo + nRepnKinds;
    p->name = VNewString (name);
    p->size = p->precision = p->min_value = p->max_value = 0.0;
    p->methods = methods;

    return nRepnKinds++;
}


/*
 *  VLookupType
 *
 *  Locate information about a named type.
 */

VRepnKind VLookupType (VStringConst name)
{
    VRepnKind repn;

    for (repn = VUnknownRepn; repn < nRepnKinds; repn++)
	if (strcmp (VRepnInfo[repn].name, name) == 0)
	    return repn;
    return VUnknownRepn;
}
