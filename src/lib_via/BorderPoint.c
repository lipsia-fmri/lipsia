/*! \file
  check if border point.

  Used in topological thinning (see Gong,Bertrand, 1997)

\par Author:
Gabriele Lohmann, MPI-CBS
*/


#include <stdio.h>
#include <stdlib.h>

/* From the Vista library: */
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>


/*!
\fn int VBorderPoint(VImage src,int b,int r,int c)
\param src  input image (bit repn)
\param b    slice address
\param r    row address
\param c    column address
*/
int
VBorderPoint(VImage src,int b,int r,int c)
{
  int adj[27];
  int i,n=0,u=0,bb,rr,cc;
  int nblack,nbands,nrows,ncols;

  nbands = VImageNBands(src);
  nrows  = VImageNRows(src);
  ncols  = VImageNColumns(src);

  i = 0;
  nblack = 0;
  for (bb=b-1; bb<=b+1; bb++) {
    for (rr=r-1; rr<=r+1; rr++) {
      for (cc=c-1; cc<=c+1; cc++) {
	if (bb >= 0 && bb < nbands
	    && rr >= 0 && rr < nrows 
	    && cc >= 0 && cc < ncols)
	  u = VPixel(src,bb,rr,cc,VBit);
	else
	  u = 0;
	adj[i++] = u;
	if (u > 0) nblack++;
      }
    }
  }

  if (nblack >= 8) return 1;
  if (nblack < 4) return 0;

  n = adj[4] + adj[10] + adj[12];
  if (n == 3) return 1;

  n = adj[4] + adj[10] + adj[14];
  if (n == 3) return 1;

  n = adj[4] + adj[12] + adj[16];
  if (n == 3) return 1;

  n = adj[4] + adj[14] + adj[16];
  if (n == 3) return 1;


  n = adj[22] + adj[10] + adj[12];
  if (n == 3) return 1;

  n = adj[22] + adj[10] + adj[14];
  if (n == 3) return 1;

  n = adj[22] + adj[12] + adj[16];
  if (n == 3) return 1;

  n = adj[22] + adj[14] + adj[16];
  if (n == 3) return 1;

  return 0;
}


