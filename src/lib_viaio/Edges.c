/*
 *  $Id: Edges.c 3177 2008-04-01 14:47:24Z karstenm $
 *
 *  This file provides basic support for edges (the VEdges class).
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
 *  Author: David Lowe, UBC Laboratory for Computational Intelligence
 */

/* From the Vista library: */
#include "viaio/Vlib.h"
#include "viaio/file.h"
#include "viaio/os.h"
#include "viaio/VEdges.h"

/* File identification string: */
VRcsId ("$Id: Edges.c 3177 2008-04-01 14:47:24Z karstenm $");


/*
 *  VCreateEdges
 *
 *  Allocates memory for a VEdges structure and initializes its fields.
 *    Initially, this contains zero edges, so each edge must still
 *    be created and added to this record.
 *  Returns a pointer to the edges if successful, NULL otherwise.
 */

VEdges VCreateEdges (int nrows, int ncolumns,
		     int nedge_fields, int npoint_fields)
{
    VEdges edges;

    /* Check parameters: */
    if (nrows < 1  ||  ncolumns < 1)
	VWarning ("VCreateEdges: Invalid number of rows or columns.");

    /* Allocate memory for the VEdges, its indices, and pixel values: */
    edges = VMalloc (sizeof (VEdgesRec));

    /* Initialize the VEdges: */
    edges->nrows = nrows;
    edges->ncolumns = ncolumns;
    edges->attributes = VCreateAttrList ();
    edges->nedge_fields = nedge_fields;
    edges->npoint_fields = npoint_fields;
    edges->nedges = edges->npoints = 0;
    edges->first = edges->last = NULL;
    edges->free = NULL;

    return edges;
}


/*
 *  VAddEdge
 *
 *  Add an edge to the given Edges record.  If the "copy" argument is
 *    TRUE, then new space is allocated to copy the points and the fields
 *    of this edge.  Otherwise, a pointer is created to their current 
 *    location.  
 *  "npoints" is the number of points in this edge, and "closed" indicates
 *    if this is a closed edge.
 */

VEdge VAddEdge (VEdges edges, VFloat *edge_fields, int npoints, VFloat *points,
		VBooleanPromoted closed, VBooleanPromoted copy)
{
    VEdge edge = VMalloc (sizeof (VEdgeRec));
    size_t fsize, psize, isize;
    int i;
    VPointer p;
    VFloat *pdata;
    
    /* Add the edge to the end of the current list of edges in order to
       maintain a consistent ordering of edges during IO. */
    if (edges->last == NULL)
	edges->first = edge;
    else
	edges->last->next = edge;
    edges->last = edge;
    edge->next = NULL;
    edges->nedges += 1;
    edges->npoints += npoints;
    edge->npoints = npoints;
    edge->closed = closed;
    isize = sizeof (VFloat *) * npoints;  /* Size of points index array. */
    
    /* If copying data, enough space is allocated to hold everything. */
    if (copy) {
#ifndef __alpha
	fsize = sizeof (VFloat) * edges->nedge_fields;
	psize = sizeof (VFloat) * npoints * edges->npoint_fields;
#else
	/* pointers must be quadword-aligned on a DEC alpha */
#define quadalign(a)	((((a)-1)/8+1)*8)
	fsize = quadalign(sizeof (VFloat) * edges->nedge_fields);
	psize = quadalign(sizeof (VFloat) * npoints * edges->npoint_fields);
#endif
	p = VMalloc (fsize + psize + isize);
	edge->free = p;
	edge->edge_fields = (VFloat *) p;
	if (fsize > 0)
	    memcpy (p, edge_fields, fsize);
	pdata = (VFloat *) ((char *) p + fsize);
	memcpy (pdata, points, psize);
	edge->point_index = (VFloat **) ((char *) p + fsize + psize);
    } else {
	p = VMalloc (isize);
	edge->free = p;
	edge->edge_fields = edge_fields;
	pdata = points;
	edge->point_index = (VFloat **) p;
    }

    /* Initialize index array into set of points. */
    for (i = 0; i < npoints; i++)
	edge->point_index[i] = pdata + i * edges->npoint_fields;

    return edge;
}


/*
 *  VCopyEdges
 *
 *  Copy a VEdges object.
 */

VEdges VCopyEdges (VEdges src)
{
    VEdges result;
    VEdge e;

    result = VCreateEdges (src->nrows, src->ncolumns, 
			   src->nedge_fields, src->npoint_fields);
    for (e = src->first; e != NULL; e = e->next)
	VAddEdge (result, e->edge_fields, e->npoints, e->point_index[0],
		  e->closed, TRUE);
    if (VEdgesAttrList (result))
	VDestroyAttrList (VEdgesAttrList (result));
    if (VEdgesAttrList (src))
	VEdgesAttrList (result) = VCopyAttrList (VEdgesAttrList (src));
    return result;
}


/*
 *  VDestroyEdges
 *
 *  Frees memory occupied by set of edges.
 */

void VDestroyEdges (VEdges edges)
{
    VEdge edge, next_edge;

    for (edge = edges->first; edge; edge = next_edge) {
	next_edge = edge->next;
	if (edge->free)
	    VFree (edge->free);
    }
    if (edges->free)
	VFree (edges->free);
    VDestroyAttrList (edges->attributes);
    VFree (edges);
}


/*
 *  VReadEdges
 *
 *  Read a Vista data file, extract the edge sets from it, and return a list
 *  of them.
 */

int VReadEdges (FILE *file, VAttrList *attributes, VEdges **edge_sets)
{
    return VReadObjects (file, VEdgesRepn, attributes,
			 (VPointer **) edge_sets);
}


/*
 *  VWriteEdges
 *
 *  Write a list of edge sets to a Vista data file.
 */

VBoolean VWriteEdges (FILE *file, VAttrList attributes,
		      int nedge_sets, VEdges edge_sets[])
{
    return VWriteObjects (file, VEdgesRepn, attributes, nedge_sets,
			  (VPointer *) edge_sets);
}
