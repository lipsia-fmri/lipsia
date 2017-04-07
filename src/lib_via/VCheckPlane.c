/*
**  test "checking plane"-condition (Tsao, and Fu's algorithm)
**  used in topological thinning
**
** G.Lohmann, Januar 1997
*/
#include <stdio.h>
#include <stdlib.h>

/* From the Vista library: */
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#define BAND 0
#define ROW  1
#define COL  2

enum MARKTYPE {VISITED, UNVISITED, EMPTY};
enum SKELTYPE {SURFACE, AXIS};

/*
** get number of 8-connected components in 3x3 neighbourhood
*/

static int nad8[9] = {3,5,3,5,8,5,3,5,3};

static int ad8[9][9] = {
  {1,3,4},     /* 0 */
  {0,2,3,4,5}, /* 1 */
  {1,4,5},     /* 2 */
  {0,1,4,6,7}, /* 3 */
  {0,1,2,3,5,6,7,8}, /* 4 */
  {1,2,4,7,8}, /* 5 */
  {3,4,7},     /* 6 */
  {3,4,5,6,8}, /* 7 */
  {4,5,7}      /* 8 */
};

int
VNumComp8(VBit src_pp[9])
{
  int npixels;
  int i,j,k,l;
  int n;
  int label;
  int dest_pp[9];

  npixels = 9;
  for (i=0; i<npixels; i++) dest_pp[i] = 0;

  label = 0;
  for (i=0; i<npixels; i++) {

    if ((src_pp[i] == 0) || (dest_pp[i] > 0)) continue;

    label++;
    dest_pp[i] = label;

    n = 1;
    while (n > 0) {
      n = 0;

      for (j=i; j<npixels; j++) {

	if ((src_pp[j] == 0) || (dest_pp[j] > 0)) continue;

	for (k=0; k<nad8[j]; k++) {
	  l = ad8[j][k];
	  if (dest_pp[l] == label) {
	    dest_pp[j] = label;
	    n = 1;
	    goto next;
	  }
	}
	next: ;
      }
    }
  }

  return (int) label;
}



int
VCheckPlane(VImage src,int b, int r, int c,
	    int plane_type,int skeltype)

{
  int i=0,n=0,u,n1,n2;
  int bb,rr,cc;
  int nbands,nrows,ncols;
  VBit plane[9];

  nbands = VImageNBands(src);
  nrows  = VImageNRows(src);
  ncols  = VImageNColumns(src);


  /* set up checking planes */

  switch(plane_type) {
    
  case BAND:
	    
    n = 0;
    i = 0;
    for (rr=r-1; rr<=r+1; rr++) {
      for (cc=c-1; cc<=c+1; cc++) {
	if (rr >= 0 && rr < nrows
	    && cc >= 0 && cc < ncols) {
	  u = VPixel(src,b,rr,cc,VBit);
	}
	else 
	  u = 0;
	if (u > 0) n++;
	plane[i++] = u;
      }
    }
    break;

  case ROW:

    n = 0;
    i = 0;
    for (bb=b-1; bb<=b+1; bb++) {
      for (cc=c-1; cc<=c+1; cc++) {
	if (bb >= 0 && bb < nbands
	    && cc >= 0 && cc < ncols) {
	  u = VPixel(src,bb,r,cc,VBit);
	}
	else 
	  u = 0;
	if (u > 0) n++;
	plane[i++] = u;
      }
    }
    break;

  case COL:

    n = 0;
    i = 0;
    for (bb=b-1; bb<=b+1; bb++) {
      for (rr=r-1; rr<=r+1; rr++) {
	if (bb >= 0 && bb < nbands 
	    && rr >= 0 && rr < nrows) {
	  u = VPixel(src,bb,rr,c,VBit);
	}
	else 
	  u = 0;
	if (u > 0) n++;
	plane[i++] = u;
      }
    }
    break;

  default:
    VError("error in Checkplane");
  }
  n--;
  if (n < 1) return 0;
  if (skeltype == SURFACE && n < 3) return 0;
/*
  if (skeltype == AXIS && n <= 1 && testvoxel == 0) return 0;
  */

  /*
  ** count number of components before and after deletion 
  */

  n1 = 1;

  /* after deletion */
  plane[4] = 0;
  n2 = VNumComp8(plane);

  if (n1 == n2) return (int) 1;
  else return (int) 0;
}


/*
** test checking plane conditions for a point (Tsao)
*/
int
VCheckPoint(VImage src,int b,int r,int c,int dir)
{
  int check;
  int skeltype = SURFACE;
 
  check = 0;
  if (dir == 0 || dir == 1) {
    check = (VCheckPlane(src,b,r,c,BAND,skeltype) 
	     && VCheckPlane(src,b,r,c,COL,skeltype));
  }
  else if (dir == 2 || dir == 3) {
    check = (VCheckPlane(src,b,r,c,BAND,skeltype) 
	     && VCheckPlane(src,b,r,c,ROW,skeltype));
  }
  else if (dir == 4 || dir == 5) {
    check = (VCheckPlane(src,b,r,c,COL,skeltype) 
	     && VCheckPlane(src,b,r,c,ROW,skeltype));
  }
  return check;
}
