/*! \file
  2D scaling using bilinear interpolation.

Scale all slices of a 3D image using slicewise bilinear interpolation.
The image center remains fixed. 

\par Author:
Gabriele Lohmann, MPI-CBS
*/


/* From the Vista library: */
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <viaio/os.h>
#include <viaio/VImage.h>

#include <stdio.h>
#include <math.h>


/*!
\fn VImage VBiLinearScale2d (VImage src,VImage dest,int dest_nrows,int dest_ncols,
                             VFloat yscale,VFloat xscale)
\brief Scale all slices of a 3D image using a 2D slicewise bilinear interpolation.
\param src         input image (any repn)
\param dest        output image (any repn)
\param dest_nrows  number of rows in output image
\param dest_ncols  number of columns in output image
\param shift[2]    translation vector (row,column)
\param scale[2]    scaling vector (row,column)

*/
VImage 
VBiLinearScale2d (VImage src,VImage dest,int dest_nrows,int dest_ncols,
		  VFloat shift[2], VFloat scale[2])
{
  int   b,r,c;
  float rp=0,cp=0,x,y,yscale,xscale;
  int   right,left,top,bottom;
  float p1,p2,p3,p4;
  float val;
  int   nrows,ncols,nbands; 
  VRepnKind repn;


  /* Extract data from source image */
  nrows  = VImageNRows(src);
  ncols  = VImageNColumns(src);
  nbands = VImageNBands(src);
  repn   = VPixelRepn(src);


  dest = VSelectDestImage("VBiLinearScale2d",dest,nbands,dest_nrows,dest_ncols,repn);
  if (! dest) return NULL;
  VFillImage(dest,VAllBands,0);


#define Top(r)    ( (int) (r) )
#define Bottom(r) ( (int) (r) + 1 )
#define Left(c)   ( (int) (c) )
#define Right(c)  ( (int) (c) + 1 )
    
  yscale = 1.0 / scale[0];
  xscale = 1.0 / scale[1];

  for (b = 0; b < nbands; b++) {

    for (r = 0; r < dest_nrows; r++) {
      rp = yscale * ((float) r - shift[0]);
      if (rp < 0 || rp >= nrows) continue;

      for (c = 0; c < dest_ncols; c++) {
	cp = xscale * ((float) c - shift[1]);
	if (cp < 0 || cp >= ncols) continue;
	
	right = Right(cp);
	left  = Left(cp);
	if (left < 0 || right >= ncols) continue;

	top     = Top(rp);
	bottom  = Bottom(rp);
	if (top < 0 || bottom >= nrows) continue;

	x = right - cp;
	y = bottom - rp;

	p1 = VGetPixel(src,b,top,left);
	p2 = VGetPixel(src,b,top,right);
	p3 = VGetPixel(src,b,bottom,left);
	p4 = VGetPixel(src,b,bottom,right);

	val = x * y * p1
	  + (1-x) * y * p2
	  + x * (1-y) * p3
	  + (1-x) * (1-y) * p4;

	if (repn == VUByteRepn) {
	  val = (int)(val + 0.5);
	  if (val <   0) val = 0;
	  if (val > 255) val = 255;
	}
	VSetPixel(dest,b,r,c,(VDouble) val);
      }
    }
  }

  VCopyImageAttrs (src, dest);
  return dest;
}
