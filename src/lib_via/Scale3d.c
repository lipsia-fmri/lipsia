#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/file.h>
#include <viaio/mu.h>

VImage VTriLinearScale3d(VImage src,VImage dest,
			 int dst_nbands,int dst_nrows,int dst_ncolumns,
			 float *shift,float *scale)
{
  int   b,r,c,bb,rr,cc;
  double c000,c100,c001,c101,c010,c110,c011,c111;
  double c00,c01,c10,c11,c0,c1;
  double x,y,z,xd,yd,zd,x0,y0,z0;
  double val;

  VRepnKind repn = VPixelRepn(src);
  int src_nrows  = VImageNRows(src);
  int src_ncols  = VImageNColumns(src);
  int src_nbands = VImageNBands(src);

  double zscale = 1.0 / scale[2];
  double yscale = 1.0 / scale[1];
  double xscale = 1.0 / scale[0];


  /* create output image */
  if (dest == NULL) {
    dest = VCreateImage(dst_nbands,dst_nrows,dst_ncolumns,repn);
  }
  if (! dest) return NULL;
  VFillImage(dest,VAllBands,0);


  /* Determines the value of each pixel in the destination image: */
  for (b=0; b<dst_nbands; b++) {
    z = zscale * ((double)b - shift[2]);
    if (z < 0 || z > src_nbands) continue;

    for (r=0; r<dst_nrows; r++) {
      y = yscale * ((double) r - shift[1]);
      if (y < 0 || y > src_nrows) continue;

      for (c=0; c<dst_ncolumns; c++) {
	x = xscale * ((double) c - shift[0]);
	if (x < 0 || x > src_ncols) continue;

	
	/* check range */
	cc = (int)x;
	rr = (int)y;
	bb = (int)z;
	if (cc < 1 || cc >= src_ncols-1) continue;
	if (rr < 1 || rr >= src_nrows-1) continue;
	if (bb < 1 || bb >= src_nbands-1) continue;

	/* function values at corners of cube */
	c000 = VGetPixel(src,bb,rr,cc);
	c100 = VGetPixel(src,bb,rr,cc+1);
	c001 = VGetPixel(src,bb+1,rr,cc);
	c101 = VGetPixel(src,bb+1,rr,cc+1);
	c010 = VGetPixel(src,bb,rr+1,cc);
	c110 = VGetPixel(src,bb,rr+1,cc+1);
	c011 = VGetPixel(src,bb+1,rr+1,cc);
	c111 = VGetPixel(src,bb+1,rr+1,cc+1);


	/* lower corner of cube */
	x0 = (int) (x);
	y0 = (int) (y);
	z0 = (int) (z);


	/* distance to corners of cube */
	xd = (x-x0);
	yd = (y-y0);
	zd = (z-z0);

	/* interpolate */
	c00 = c000*(1.0-xd) + c100*xd;
	c01 = c001*(1.0-xd) + c101*xd;
	c10 = c010*(1.0-xd) + c110*xd;
	c11 = c011*(1.0-xd) + c111*xd;

	c0 = c00*(1.0-yd) + c10*yd;
	c1 = c01*(1.0-yd) + c11*yd;
	val = c0*(1.0-zd) + c1*zd;

	VSetPixel(dest,b,r,c,(double)val);
      }
    }
  }

  return dest;
}
