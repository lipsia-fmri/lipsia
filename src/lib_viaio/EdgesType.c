/*
 *  $Id: EdgesType.c 3177 2008-04-01 14:47:24Z karstenm $
 *
 *  This file contains methods for the edge set (VEdges) type.
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
#include "viaio/os.h"
#include "viaio/VEdges.h"

/* File identification string: */
VRcsId ("$Id: EdgesType.c 3177 2008-04-01 14:47:24Z karstenm $");


/*
 *  Table of methods.
 */

/* Later in this file: */
static VDecodeMethod VEdgesDecodeMethod;
static VEncodeAttrMethod VEdgesEncodeAttrMethod;
static VEncodeDataMethod VEdgesEncodeDataMethod;

/* Used in Type.c to register this type: */
VTypeMethods VEdgesMethods = {
  (VCopyMethod *) VCopyEdges,		/* copy a VEdges */
  (VDestroyMethod *) VDestroyEdges,	/* destroy a VEdges */
  VEdgesDecodeMethod,			/* decode a VEdges's value */
  VEdgesEncodeAttrMethod,		/* encode a VEdges's attr list */
  VEdgesEncodeDataMethod		/* encode a VEdges's binary data */
};


/*
 *  VEdgesDecodeMethod
 *
 *  The "decode" method registered for the "edges" type.
 *  Convert an attribute list plus binary data to a VEdges object.
 */

static VPointer VEdgesDecodeMethod (VStringConst name, VBundle b)
{
  VEdges edges;
  VLong nrows, ncolumns, nedge_fields, npoint_fields, nedges, npoints;
  VAttrList list;
  int npnts, closed, i, nels;
  size_t length;
  VFloat *p;

#define Extract(name, dict, locn, required)	\
	VExtractAttr (b->list, name, dict, VLongRepn, & locn, required)
    
  /* Extract the required attribute values for edges. */
  if (! Extract (VNRowsAttr, NULL, nrows, TRUE) ||
      ! Extract (VNColumnsAttr, NULL, ncolumns, TRUE) ||
      ! Extract (VNEdgeFieldsAttr, NULL, nedge_fields, TRUE) ||
      ! Extract (VNPointFieldsAttr, NULL, npoint_fields, TRUE) ||
      ! Extract (VNEdgesAttr, NULL, nedges, TRUE) ||
      ! Extract (VNPointsAttr, NULL, npoints, TRUE))
    return NULL;
  if (npoint_fields <= 0) {
    VWarning ("VEdgesReadDataMethod: Bad edges file attributes");
    return NULL;
  }

  /* Create the edges data structure. */
  edges = VCreateEdges ((int) nrows, (int) ncolumns,
			(int) nedge_fields, (int) npoint_fields);
  if (! edges)
    return NULL;

  /* Give it whatever attributes remain: */
  list = VEdgesAttrList (edges);
  VEdgesAttrList (edges) = b->list;
  b->list = list;

  /* Check that the expected amount of binary data was read: */
  nels = (nedges * (nedge_fields + 1) + npoints * npoint_fields);
  length = nels * VRepnPrecision (VFloatRepn) / 8;
  if (length != b->length) {
    VWarning ("VEdgesDecodeMethod: %s image has wrong data length", name);
  Fail:   VDestroyEdges (edges);
    return NULL;
  }

  /* Allocate storage for the edges binary data: */
  edges->free = VMalloc (length = nels * sizeof (VFloat));

  /* Unpack the binary data: */
  if (! VUnpackData (VFloatRepn, nels, b->data, VMsbFirst, & length,
		     & edges->free, NULL))
    goto Fail;

  /* Create an edge record for each edge in this list.  The first number
     in each edge's data encodes number of points and closed flag. */
  p = (VFloat *) edges->free;
  for (i = 0; i < nedges; i++) {
    closed = (int) (*p < 0.0);
    npnts = (int) (closed ? - *p++ : *p++);
    VAddEdge (edges, p, npnts, p + nedge_fields, closed, FALSE);
    p += nedge_fields + npnts * npoint_fields;
  }
  return edges;

#undef Extract
}


/*
 *  VEdgesEncodeAttrMethod
 *
 *  The "encode_attrs" method registered for the "edges" type.
 *  Encode an attribute list value for a VEdges object.
 */

static VAttrList VEdgesEncodeAttrMethod (VPointer value, size_t *lengthp)
{
  VEdges edges = value;
  VAttrList list;

  /* Temporarily prepend several attributes to the edge set's list: */
  if ((list = VEdgesAttrList (edges)) == NULL)
    list = VEdgesAttrList (edges) = VCreateAttrList ();
  VPrependAttr (list, VNColumnsAttr, NULL, VLongRepn,
		(VLong) edges->ncolumns);
  VPrependAttr (list, VNRowsAttr, NULL, VLongRepn,
		(VLong) edges->nrows);
  VPrependAttr (list, VNPointFieldsAttr, NULL, VLongRepn,
		(VLong) edges->npoint_fields);
  VPrependAttr (list, VNEdgeFieldsAttr, NULL, VLongRepn,
		(VLong) edges->nedge_fields);
  VPrependAttr (list, VNPointsAttr, NULL, VLongRepn,
		(VLong) edges->npoints);
  VPrependAttr (list, VNEdgesAttr, NULL, VLongRepn,
		(VLong) edges->nedges);

  /* Compute the file space needed for the edges's binary data: */
  *lengthp = VRepnPrecision (VFloatRepn) / 8 *
    (edges->nedges * (1 + edges->nedge_fields) +
     edges->npoints * edges->npoint_fields);

  return list;
}


/*
 *  VEdgesEncodeDataMethod
 *
 *  The "encode_data" method registered for the "edges" type.
 *  Encode the edge and point fields for a VEdges object.
 */

static VPointer VEdgesEncodeDataMethod (VPointer value, VAttrList list,
					size_t length, VBoolean *free_itp)
{
  VEdges edges = value;
  VAttrListPosn posn;
  VEdge e;
  float pcount;
  size_t len;
  VPointer ptr, p;

  /* Remove the attributes prepended by the VEdgesEncodeAttrsMethod: */
  for (VFirstAttr (list, & posn);
       strcmp (VGetAttrName (& posn), VNColumnsAttr) != 0;
       VDeleteAttr (& posn)) ;
  VDeleteAttr (& posn);

  /* Allocate a buffer for the encoded data: */
  p = ptr = VMalloc (length);

  /* Pack each edge: */
  for (e = edges->first; e != NULL; e = e->next) {

    /* Pack the number of points (negated if closed): */
    pcount = (float) e->npoints;
    if (e->closed)
      pcount = - pcount;
    len = length;
    if (! VPackData (VFloatRepn, 1, & pcount, VMsbFirst, & len, & p, NULL))
      return NULL;
    p = (char *) p + len;
    length -= len;

    /* Pack any edge fields: */
    if (edges->nedge_fields > 0) {
      len = length;
      if (! VPackData (VFloatRepn, edges->nedge_fields, e->edge_fields,
		       VMsbFirst, & len, & p, NULL))
	return NULL;
      p = (char *) p + len;
      length -= len;
    }

    /* Pack point fields: */
    len = length;
    if (! VPackData (VFloatRepn, edges->npoint_fields * e->npoints,
		     e->point_index[0], VMsbFirst, & len, & p, NULL))
      return NULL;
    p = (char *) p + len;
    length -= len;
  }

  *free_itp = TRUE;
  return ptr;
}
