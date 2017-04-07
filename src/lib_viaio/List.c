/*
 *  $Id: List.c 3177 2008-04-01 14:47:24Z karstenm $
 *
 *  This file contains the routines for manipulating VList instances.
 *
 *  [Thanks to Murray Goldberg (goldberg@cs.ubc.ca)]
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
#include "viaio/os.h"
#include "viaio/VList.h"



/* Later in this file: */
static VNodePtrType MakeNode (VPointer item,
			      VNodePtrType prev, VNodePtrType next);


/*
 * MakeNode
 *
 * Make a node
 */

static VNodePtrType MakeNode (VPointer item, VNodePtrType prev,
			     VNodePtrType next)
{
    VNodePtrType result = VMalloc (sizeof (struct V_Node));

    result->item = item;
    result->prev = prev;
    result->next = next;

    return result;
}


/*
 * VListCreate
 *
 * Make a new, empty list, and returns its reference.
 */

VList VListCreate (void)
{
    VList vlist = VMalloc (sizeof (struct V_List));
    VNodePtrType dummy_head, dummy_tail;

    dummy_head = VMalloc (sizeof (struct V_Node));

    dummy_tail = VMalloc (sizeof (struct V_Node));

    dummy_head->item = NULL;
    dummy_head->prev = NULL;
    dummy_head->next = dummy_tail;

    dummy_tail->item = NULL;
    dummy_tail->prev = dummy_head;
    dummy_tail->next = NULL;

    vlist->head	   = dummy_head;
    vlist->tail	   = dummy_tail;
    vlist->current = dummy_head;
    vlist->count   = 0;

    return vlist;
}


/*
 * VListFirst
 *
 * Return a pointer to the first item in vlist, and
 * make the first item the current item.
 */

VPointer VListFirst (VList vlist)
{
    if ( vlist->count == 0 )   /* empty vist, move beyond beginning */
	vlist->current = vlist->head;
    else		       /* vlist not empty, move to beginning */
	vlist->current = vlist->head->next;

    return vlist->current->item;
}


/*
 * VListLast
 *
 * Return a pointer to the last item in vlist, and
 * make the last item the current item.
 */

VPointer VListLast (VList vlist)
{
    if ( vlist->count == 0 )   /* empty vlist, move beyond end */
	vlist->current = vlist->tail;
    else		       /* vlist not empty, move to end */
	vlist->current = vlist->tail->prev;

    return vlist->current->item;
}


/*
 * VListNext
 *
 * Advance vlist's current item by one, return the
 * new current item. Return NULL if the new current
 * item is beyond the end of vlist.
 */

VPointer VListNext (VList vlist)
{
    if ( vlist->current == vlist->tail )
	/* already beyond end, no action */
	;
    else   /* move to next node */
	vlist->current = vlist->current->next;

    return vlist->current->item;
}


/*
 * VListPrev
 *
 * Back up vlist's current item by one, return the
 * new current item. Return NULL if the new current
 * item is before the beginning of vlist.
 */

VPointer VListPrev (VList vlist)
{
    if ( vlist->current == vlist->head )
	/* already before beginning, no action */
	;
    else   /* move to previous node */
	vlist->current = vlist->current->prev;

    return vlist->current->item;
}


/*
 * VListAdd
 *
 * Add item to vlist immediately after the current
 * item, and make item the current item. If the
 * current pointer is before the beginning of vlist,
 * item is added at the beginning. If the current
 * pointer is beyond the end of vlist, item is
 * added at the end.
 */

void VListAdd (VList vlist, VPointer item)
{
    VNodePtrType add_me;

    if ( vlist->current == vlist->tail )
	/* current pointer beyond end, add to end */
	vlist->current = vlist->tail->prev;

    add_me = MakeNode (item, vlist->current, vlist->current->next);

    add_me->prev->next = add_me;
    add_me->next->prev = add_me;

    vlist->current = add_me;
    vlist->count++;
}


/*
 * VListInsert
 *
 * Add item to vlist immediately before the current
 * item, and make item the current item. If the
 * current pointer is before the beginning of vlist,
 * item is added at the beginning. If the current
 * pointer is beyond the end of vlist, item is
 * added at the end.
 */

void VListInsert (VList vlist, VPointer item)
{
    VNodePtrType add_me;

    if ( vlist->current == vlist->head )
	/* current pointer before beginning, add to beginning */
	vlist->current = vlist->head->next;

    add_me = MakeNode (item, vlist->current->prev, vlist->current);

    add_me->prev->next = add_me;
    add_me->next->prev = add_me;

    vlist->current = add_me;
    vlist->count++;
}


/*
 * VListAppend
 *
 * Add item to the end of vlist, and make item
 * the current item.
 */

void VListAppend (VList vlist, VPointer item)
{
    vlist->current = vlist->tail;   /* move beyond end */
    VListAdd (vlist, item);
}


/*
 * VListPrepend
 *
 * Add item to the beginning of vlist, and make
 * item the current item.
 */

void VListPrepend (VList vlist, VPointer item)
{
    vlist->current = vlist->head;   /* move before beginning */
    VListAdd (vlist, item);
}


/*
 * VListRemove
 *
 * Return current item and take it out of vlist. Make
 * the next item the current one.
 */

VPointer VListRemove (VList vlist)
{
    VPointer return_me;
    VNodePtrType free_me;

    return_me = vlist->current->item;

    if ((vlist->current == vlist->tail)
	|| (vlist->current == vlist->head) )
	/* current pointer before beginning or beyond end, no action */
	;
    else {  /* free current node */

	vlist->current->prev->next = vlist->current->next;
	vlist->current->next->prev = vlist->current->prev;
	free_me = vlist->current;
	vlist->current = vlist->current->next;

	VFree (free_me);
	vlist->count--;
    }

    return return_me;
}


/*
 * VListConcat
 *
 * Add vlist2 to the end of vlist1. The current
 * pointer is set to the current pointer of vlist1.
 * vlist2 no longer exists after the operation.
 */

void VListConcat (VList vlist1, VList vlist2)
{
    VNodePtrType free_me, free_me_too;

    free_me = vlist1->tail;
    free_me_too = vlist2->head;

    vlist1->tail->prev->next = vlist2->head->next;
    vlist2->head->next->prev = vlist1->tail->prev;

    if ( vlist1->current == vlist1->tail )
	/* current pointer of vlist1 points beyond end,
	   set it to first node of vlist2 */
	vlist1->current = vlist2->head->next;

    vlist1->tail = vlist2->tail;
    vlist1->count += vlist2->count;

    VFree (free_me);
    VFree (free_me_too);
    VFree (vlist2);
}


/*
 * VListDestroy
 *
 * Delete vlist. item_free is a pointer to a routine
 * that frees an item.
 */

void VListDestroy (VList vlist, void (*item_free) ())
{
   VPointer free_me;

   vlist->current = vlist->head->next;
   while ( vlist->current != vlist->tail )
   {
       free_me = VListRemove (vlist);
       (*item_free)(free_me);
   }

   VFree (vlist->head);
   VFree (vlist->tail);
   VFree (vlist);
}


/*
 * VListTrim
 *
 * Return last item and take it out of vlist. Make
 * the new last item the current one.
 */

VPointer VListTrim (VList vlist)
{
    VPointer return_me;

    return_me = VListLast(vlist);
    VListRemove (vlist);
    VListLast (vlist);

    return return_me;
}


/*
 * VListSearch
 *
 * Searche vlist starting at the current item until
 * the end is reached or a match is found.
 */

VPointer VListSearch (VList vlist, int (*comp) (), VPointer comp_arg)
{
    if ( vlist->current == vlist->head )
	/* before beginning, go to next node */
	vlist->current = vlist->current->next;

    while ( vlist->current != vlist->tail ) {
	if ( (*comp)(vlist->current->item, comp_arg) )
	    /* a match is found */
	    return (vlist->current->item);
	else
	    vlist->current = vlist->current->next;
    }

    /* no match */
    return vlist->current->item;
}
