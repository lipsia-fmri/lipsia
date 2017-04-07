/*
 *  $Id: VEdges.h 726 2004-03-08 13:12:45Z lohmann $
 *
 *  Definitions associated with edges: their representation in files and
 *  in memory, and operations that can be performed with them.
 */

#ifndef V_VEdges_h
#define V_VEdges_h 1

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


#ifdef __cplusplus
extern "C" {
#endif


/*  Structures for representing edges in memory.
 *  -------------------------------------------
 *
 *  The edges are stored as a linked list of edge records, as this makes
 *     it easy to allocate and add individual edges.
 *  The point data for each edge are stored in a 2-D array of floating point
 *     values.
 */

typedef struct V_EdgesRec {
    int nrows;			/* number of rows */
    int ncolumns;		/* number of columns */
    VAttrList attributes;	/* list of other attributes */
    int nedge_fields;		/* number of fields in each edge record */
    int npoint_fields;		/* number of fields in each point record */
    int nedges;			/* number of edges */
    int npoints;		/* total number of points */
    struct VEdgeStruct *first;	/* first edge in linked list of edges */
    struct VEdgeStruct *last;	/* last edge in linked list of edges */
    VPointer free;		/* free this storage when destroying edges */
} VEdgesRec;


typedef struct VEdgeStruct {
    struct VEdgeStruct *next;	/* next edge in linked list of edges */
    VFloat *edge_fields;	/* vector of field entries for this edge */
    VBoolean closed;		/* indicates closed edge (a loop) */
    int npoints;		/* number of points in this edge */
    VFloat **point_index;	/* pointers to start of each point */
    VPointer free;		/* free this storage when destroying edges */
} VEdgeRec, *VEdge;


/*
 *  Attributes used to represent edges
 *  ----------------------------------
 */

/* Attribute type names: */
#define VEdgesAttr		"edges"
#define VNEdgeFieldsAttr	"nedge_fields"
#define VNPointFieldsAttr	"npoint_fields"
#define VNEdgesAttr		"nedges"
#define VNPointsAttr		"npoints"


/*
 *  Macros for accessing edges attributes in memory.
 *  -----------------------------------------------
 */

#define VEdgesNRows(edges)	((edges)->nrows)

#define VEdgesNColumns(edges)	((edges)->ncolumns)

#define VEdgesAttrList(edges)	((edges)->attributes)

#define VNEdgeFields(edges)	((edges)->nedge_fields)

#define VNPointFields(edges)	((edges)->npoint_fields)

#define VNEdges(edges)		((edges)->nedges)

#define VFirstEdge(edges)	((edges)->first)

#define VNextEdge(edge)		((edge)->next)

#define VEdgeExists(edge)	((edge) != NULL)

#define VEdgeFields(edge)	((edge)->edge_fields)

#define VEdgeNPoints(edge)	((edge)->npoints)

#define VEdgeClosed(edge)	((edge)->closed)

#define VEdgePointArray(edge)	((edge)->point_index)

/* Following are old macro names which should no longer be used.
   They can be removed in a future version once all of the documentation
   is in place and has been announced. */
#define VEdgesCount(edges)	((edges)->nedges)
#define VEdgePoints(edge)	((edge)->point_index)
#define VEdgesEdgeFields(edges) ((edges)->nedge_fields)
#define VEdgesPointFields(edges) ((edges)->npoint_fields)
#define VEdgesRows(edges)	((edges)->nrows)
#define VEdgesColumns(edges)	((edges)->ncolumns)
#define VEdgePointCount(edge)	((edge)->npoints)


/*
 *  Declarations of library routines.
 */

/* From Edges.c: */

extern VEdges VCreateEdges (
#if NeedFunctionPrototypes
    int			/* nrows */,
    int			/* ncols */,
    int			/* nedge_fields */,
    int			/* npoint_fields */
#endif
);

extern VEdge VAddEdge (
#if NeedFunctionPrototypes
    VEdges		/* edges */,
    VFloat *		/* edge_fields */,
    int			/* npoints */,
    VFloat *		/* points */,
    VBooleanPromoted	/* closed */,
    VBooleanPromoted	/* copy */
#endif
);


extern VEdges VCopyEdges (
#if NeedFunctionPrototypes
    VEdges		/* edges */
#endif
);

extern void VDestroyEdges (
#if NeedFunctionPrototypes
    VEdges		/* edges */
#endif
);

/* From EdgesIO.c: */

extern int VReadEdges (
#if NeedFunctionPrototypes
    FILE *		/* file */,
    VAttrList *		/* attributes */,
    VEdges **		/* edge_sets */
#endif
);

extern VBoolean VWriteEdges (
#if NeedFunctionPrototypes
    FILE *		/* file */,
    VAttrList		/* attributes */,
    int			/* nedge_sets */,
    VEdges *		/* edge_sets */
#endif
);

/* From EdgesToPS.c: */

extern VBoolean VEdgesToPS (
#if NeedFunctionPrototypes
    FILE *		/* f */,
    VEdges		/* edge_set */,
    VBooleanPromoted	/* endpoints */
#endif
);

/* From Link.c: */

extern VEdges VLinkImage (
#if NeedFunctionPrototypes
    VImage		/* image */,
    int			/* low */,
    int			/* high */,
    int			/* minlength */
#endif
);

/* From SegEdges.c: */

extern VEdges VSegEdgesIntoLines (
#if NeedFunctionPrototypes
    VEdges		/* edge_set */,
    double 		/* accuracy */,
    int			/* granularity */,
    double 		/* magnitude */,
    int			/* product */
#endif
);

#ifdef __cplusplus
}
#endif

#endif /* V_VEdges_h */
