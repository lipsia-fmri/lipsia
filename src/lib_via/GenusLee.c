/*! \file
  3D Euler characteristic using Lee's Algorithm

\par Reference:
 C.N. Lee, A. Rosenfeld (1987).
 ICCV 1987, London, pp.567-571.

\par Author:
 Gabriele Lohmann, MPI-CBS
*/


/* From the Vista library: */
#include <viaio/Vlib.h>
#include <viaio/file.h>
#include <viaio/mu.h>
#include <viaio/option.h>
#include <viaio/VImage.h>

/* From the standard C library: */
#include <stdio.h>
#include <stdlib.h>

#define NOT(x) ((x) > 0 ? 0 : 1)

/*!
  \fn int VGenusLee (VImage src,VShort neighb)
  \brief compute the 3D Euler characteristic
  \param src   input image (bit repn)
  \param neighb adjacency type (6 or 26).
*/
int
VGenusLee (VImage src,VShort neighb)
{
  int b,r,c,nbands,ncols,nrows,npixels;
  VBit *tmpbuf;
  int x,x1,x2,x3,x4,x5,x6,x7;
  int nx1,nx2,nx3,nx4,nx5,nx6;
  int psi,euler;
  int i;
  VImage tmp=NULL;

  nrows  = VImageNRows (src) + 2;
  ncols  = VImageNColumns (src) + 2;
  nbands = VImageNBands (src) + 2;
  npixels = nbands * ncols * nrows;

  tmp = VCreateImage(nbands,nrows,ncols,VBitRepn);
  VFillImage(tmp,VAllBands,0);

  for (b=1; b<nbands-1; b++) {
    for (r=1; r<nrows-1; r++) {
      for (c=1; c<ncols-1; c++) {
	VPixel(tmp,b,r,c,VBit) = VPixel(src,b-1,r-1,c-1,VBit);
      }
    }
  }

  if (neighb == 26) {
    tmpbuf  = (VBit *) VPixelPtr(tmp,0,0,0);
    for (i=0; i<npixels; i++) {
      x = *tmpbuf;
      *tmpbuf++ = NOT(x);
    }
  }


  euler = psi = 0;

  for (b=1; b<nbands-1; b++) {
    for (r=1; r<nrows-1; r++) {
      for (c=1; c<ncols-1; c++) {

	x = VPixel(tmp, b, r, c, VBit);
	if (x == 0) continue;
		
	x1 = VPixel(tmp, b, r+1, c, VBit);
	x2 = VPixel(tmp, b, r, c+1, VBit);
	x3 = VPixel(tmp, b+1, r, c,VBit);
	x4 = VPixel(tmp, b, r+1, c+1, VBit);
	x5 = VPixel(tmp, b+1, r+1, c, VBit);
	x6 = VPixel(tmp, b+1, r, c+1, VBit);
	x7 = VPixel(tmp, b+1, r+1, c+1, VBit);

	nx1 = NOT(x1);
	nx2 = NOT(x2);
	nx3 = NOT(x3);
	nx4 = NOT(x4);
	nx5 = NOT(x5);
	nx6 = NOT(x6);
	
	if (nx1 * nx2 * nx3 == 1) {
	  psi = 1;
	}

	else if ( x7 == 0 || nx4 + nx5 + nx6 > 1) {
	  psi = x1 * x2 * x3  -  x1 * x2 * nx4  -  x1 * x3 * nx5
	    - x2 * x3 * nx6;
	}

	else {
	  psi = - x1 * x2 * nx3 * nx4
	    - x1 * nx2 * x3 * nx5  - nx1 * x2 * x3 * nx6;
	}
	euler += psi;
	/*	VPixel(tmp,b,r,c,VBit) = 0; */
      }
    }
  }

  VDestroyImage(tmp);
  return euler;
}


