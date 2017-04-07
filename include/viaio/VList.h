/*
 *  $Id: VList.h 726 2004-03-08 13:12:45Z lohmann $
 *
 *  Definitions associated with VList.
 */

#ifndef V_VList_h
#define V_VList_h 1

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

/* From the standard C library: */
#include <stdio.h>



#ifdef __cplusplus
extern "C" {
#endif


/*
 *  Declarations of data structure.
 */

/* List element: */
typedef struct V_Node *VNodePtrType;
struct V_Node {
    VPointer item;		/* pointer to data item */
    VNodePtrType prev;		/* pointer to previous node */
    VNodePtrType next;		/* pointer to next node */
};

/* List head: */
typedef struct V_List {
    VNodePtrType current;	/* pointer to current node */
    VNodePtrType head;		/* pointer to head node */
    VNodePtrType tail;		/* pointer to tail node */
    int count;			/* number of nodes in VList */
} *VList;


/*
 *  Definitions of macros.
 */

#define VListCount(vlist)	((vlist)->count)
#define VListCurr(vlist)	((vlist)->current->item)
#define VListGetCurr(vlist)	((vlist)->current)
#define VListSetCurr(vlist,curr) ((void) ((vlist)->current = (curr)))


/*
 *  Declarations of library routines.
 */

/* From List.c: */

extern VList VListCreate (
#if NeedFunctionPrototypes
    void
#endif
);

extern VPointer VListFirst (
#if NeedFunctionPrototypes
    VList		/* vlist */
#endif
);

extern VPointer VListLast (
#if NeedFunctionPrototypes
    VList		/* vlist */
#endif
);

extern VPointer VListNext (
#if NeedFunctionPrototypes
    VList		/* vlist */
#endif
);

extern VPointer VListPrev (
#if NeedFunctionPrototypes
    VList		/* vlist */
#endif
);

extern void VListAdd (
#if NeedFunctionPrototypes
    VList		/* vlist */,
    VPointer		/* item */
#endif
);

extern void VListInsert (
#if NeedFunctionPrototypes
    VList		/* vlist */,
    VPointer		/* item */
#endif
);

extern void VListAppend (
#if NeedFunctionPrototypes
    VList		/* vlist */,
    VPointer		/* item */
#endif
);

extern void VListPrepend (
#if NeedFunctionPrototypes
    VList		/* vlist */,
    VPointer		/* item */
#endif
);

extern VPointer VListRemove (
#if NeedFunctionPrototypes
    VList		/* vlist */
#endif
);

extern void VListConcat (
#if NeedFunctionPrototypes
    VList		/* vlist1 */,
    VList		/* vlist2 */
#endif
);

extern void VListDestroy (
#if NeedFunctionPrototypes
    VList		/* vlist */,
    void (*) (
#if NeedNestedPrototypes
	VPointer	    /* opaque_object */
#endif
	      )		/* item_free */
#endif
);

extern VPointer VListTrim (
#if NeedFunctionPrototypes
    VList		/* vlist */
#endif
);

extern VPointer VListSearch (
#if NeedFunctionPrototypes
    VList		/* vlist */,
    int (*) ()		/* comp */,
    VPointer		/* comp_arg */
#endif
);

#ifdef __cplusplus
}
#endif

#endif /* V_VList_h */
