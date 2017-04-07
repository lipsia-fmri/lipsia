/*! \file
  2D connected component labelling

Each foreground voxel receives a label indicating
its membership in a connected component.
The algorithm is based on depth first search.

\par Reference:
G. Lohmann (1998). "Volumetric Image Analysis",
John Wiley & Sons, Chichester, England.

\par Author:
Gabriele Lohmann, MPI-CBS
*/


/* From the Vista library: */
#include <viaio/Vlib.h>
#include <viaio/mu.h>

/* From the standard C libaray: */
#include <stdio.h>
#include <stdlib.h>
#include <via/via.h>


typedef struct {
  Pixel *A;
  int front,rear;
} Queue2D;

static Queue2D queue2D;
static int QueueSize2D;
static int msize2D;
typedef int BOOLEAN;

static void QueueClear2D(Queue2D *);
static BOOLEAN QueueEmpty2D(Queue2D *);
static BOOLEAN deQueue2D(Queue2D *, Pixel *);
static BOOLEAN enQueue2D(Queue2D *, Pixel);


#define VPutPixel(image,b,r,c,repn,value) \
  if (repn == VShortRepn) VPixel((image),(b),(r),(c),VShort) = (value); \
  else VPixel((image),(b),(r),(c),VUByte) = (value);


/*!
\fn VImage VLabelImage2d(VImage src, VImage dest,int neighb,VRepnKind repn,int *numlabels)
\param src  input image (bit repn)
\param dest output image (ubyte repn)
\param neighb adjacency type (4 or 8).
\param repn output pixel representation (VUByteRepn or VShortRepn)
\param numlabels ptr to the number of labels found.
*/

VImage
VLabelImage2d(VImage src, VImage dest,long neighb,VRepnKind repn,long *numlabels)
{
  long i,nbands,nrows,ncols;
  Pixel v,vv;
  long label=0,n,nblack=0;
  long r0,r1,c0,c1,b,r,c,rr,cc;
  long ra[6],ca[6],m;
  
  msize2D = 0;
  if (VPixelRepn(src) != VBitRepn) 
    VError("Input image must be of type VBit");

  nbands  = VImageNBands(src);
  nrows   = VImageNRows(src);
  ncols   = VImageNColumns(src);


  switch(repn) {
  case VUByteRepn:
    dest = VSelectDestImage("VLabel2d",dest,nbands,nrows,ncols,VUByteRepn);
    if (! dest) return NULL;
    VFillImage(dest,VAllBands,0);
    break;

  case VShortRepn:
    dest = VSelectDestImage("VLabel2d",dest,nbands,nrows,ncols,VShortRepn);
    if (! dest) return NULL;
    VFillImage(dest,VAllBands,0);
    break;

  default:
    VError("Output image representation must be either ubyte or short.");
  }


  /*
  ** set up Queue
  */
  nblack = nrows * ncols;  
  QueueSize2D = nblack;
  queue2D.A = (Pixel *) VMalloc(sizeof(Pixel) * QueueSize2D);
  if (numlabels != NULL) *numlabels = 0;

  /*
  ** breadth first search
  */
  for (b=0; b<nbands; b++) {

    label = 0;
    switch(neighb) {

    /*
    ** 8 - adjacency
    */
    case 8:

      for (r=0; r<nrows; r++) {
	for (c=0; c<ncols; c++) {
	  if (VPixel(src,b,r,c,VBit) == 0) continue;

	  if (repn == VUByteRepn 
	      && VPixel(dest,b,r,c,VUByte) > 0) continue;
	  if (repn == VShortRepn 
	      && VPixel(dest,b,r,c,VShort) > 0) continue;

	  QueueClear2D(&queue2D);

	  label++;
	  if (label >= VPixelMaxValue (dest)) {
	    VWarning("Number of labels exceeds maximum (%d)",label);
	    return dest;
	  }
	  n = 1;

	  v.r = r;
	  v.c = c;
	  if (enQueue2D(&queue2D,v) == FALSE) VError(" error in enQueue2D");
	  VPutPixel(dest,b,r,c,repn,label);

	  while (! QueueEmpty2D(&queue2D)) {

	    deQueue2D(&queue2D,&v);

	    r0 = (v.r > 0) ? v.r - 1 : 0;
	    r1 = (v.r + 1 < nrows) ? v.r + 1 : nrows - 1;
	    for (rr=r0; rr<=r1; rr++) {

	      c0 = (v.c > 0) ? v.c - 1 : 0;
	      c1 = (v.c + 1 < ncols) ? v.c + 1 : ncols - 1;
	      for (cc=c0; cc<=c1; cc++) {

		if (rr == v.r && cc == v.c) continue;
		if (VPixel(src,b,rr,cc,VBit) == 0) continue;
		if (repn == VUByteRepn
		    && VPixel(dest,b,rr,cc,VUByte) > 0) continue;
		if (repn == VShortRepn 
		    && VPixel(dest,b,rr,cc,VShort) > 0) continue;
		vv.r = rr;
		vv.c = cc;
		if (enQueue2D(&queue2D,vv) == FALSE) VError(" error in enQueue2D");
		VPutPixel(dest,b,rr,cc,repn,label);
		n++;
	      }
	    }
	  }
	}
      }
      break;

    /*
    ** 4 -adjacency
    */
    case 4:

      for (r=0; r<nrows; r++) {
	for (c=0; c<ncols; c++) {
	  if (VPixel(src,b,r,c,VBit) == 0) continue;

	  if (repn == VUByteRepn 
	      && VPixel(dest,b,r,c,VUByte) > 0) continue;
	  if (repn == VShortRepn 
	      && VPixel(dest,b,r,c,VShort) > 0) continue;
	  
	  QueueClear2D(&queue2D);
	
	  label++;
	  if (label >= VPixelMaxValue (dest)) {
	    VWarning("Number of labels exceeds maximum (%d)",label);
	    VFree(queue2D.A);
	    return dest;
	  }
	  n = 1;

	  v.r = r;
	  v.c = c;
	  if (enQueue2D(&queue2D,v) == FALSE) VError(" error in enQueue2D");
	  VPutPixel(dest,b,r,c,repn,label);

	  while (! QueueEmpty2D(&queue2D)) {

	    deQueue2D(&queue2D,&v);

	    for (i=0; i<4; i++) {
	      ra[i] = v.r;
	      ca[i] = v.c;
	    }
	    m = 0;
	    if (v.r > 0) ra[m++] = v.r - 1;
	    if (v.c > 0) ca[m++] = v.c - 1;

	    if (v.r < nrows - 1)  ra[m++] = v.r + 1;
	    if (v.c < ncols - 1)  ca[m++] = v.c + 1;

	    for (i=0; i<m; i++) {
	      if (VPixel(src,b,ra[i],ca[i],VBit) == 0) continue;
	      if (repn == VUByteRepn 
		  && VPixel(dest,b,ra[i],ca[i],VUByte) > 0) continue;
	      if (repn == VShortRepn 
		  && VPixel(dest,b,ra[i],ca[i],VShort) > 0) continue;

	      VPutPixel(dest,b,ra[i],ca[i],repn,label);
	      vv.r = ra[i];
	      vv.c = ca[i];
	      if (enQueue2D(&queue2D,vv) == FALSE) VError(" error in enQueue2D");
	      n++;
	    }
	  }
	}
      }
    
      break;
    
    default:
      ;
    }
  }

  if (queue2D.A) VFree(queue2D.A);

  if (numlabels != NULL) *numlabels = label;
  VCopyImageAttrs (src, dest);
  return dest;
}


/*
** array implementation of queues, following
** Aho, Ullman, "Foundations of Computer Science in C"
** Computer Science Press, 1995, pp. 319 ff
**
*/

#define TRUE  1
#define FALSE 0

static void QueueClear2D(Queue2D *pQ)
{
  pQ->front = 0;
  pQ->rear  = 0;
}

static BOOLEAN
QueueEmpty2D(Queue2D *pQ) 
{
  return (pQ->front == pQ->rear);
}


static BOOLEAN
deQueue2D(Queue2D *pQ, Pixel *pe)
{
  if (pQ->front == pQ->rear) {
    fprintf(stderr," deQueue2D: empty\n");
    return FALSE;
  }
  else {
    (*pe) = pQ->A[(pQ->rear)++];
    return TRUE;
  }
}

static BOOLEAN 
enQueue2D(Queue2D *pQ, Pixel e)
{
  if (pQ->front > msize2D) msize2D = pQ->front;
  
  if (pQ->front < QueueSize2D - 1) {
    pQ->A[(pQ->front)++] = e;
    return TRUE;
  }
  else if (pQ->rear > 2) {
    pQ->A[--(pQ->rear)] = e;
    return TRUE;
  }
  else {
    QueueSize2D += QueueSize2D * 0.2;
    fprintf(stderr," realloc: %d\n",QueueSize2D);
    queue2D.A = (Pixel *) VRealloc(queue2D.A,sizeof(Pixel) * QueueSize2D);
    pQ->A[(pQ->front)++] = e;
    return TRUE;
  }
}
