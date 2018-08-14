/*! \file
  3D grey level morphology.

\par Reference:
  P. Maragos, R.W. Schafer (1990):
  "Morphological Systems for multidimensional signal processing",
  Proc. of the IEEE, Vol. 78, No. 4, pp. 690--709.

\par Author:
 Gabriele Lohmann, MPI-CBS
*/

#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <stdio.h>
#include <math.h>


/*!
  \fn VImage VGreyDilation3d(VImage src,VImage se,VImage dest)
  \brief 3D greylevel morphological dilation
  \param src    input image (any repn)
  \param se     raster image containing the structural element (bit repn)
  \param dest   output image (any repn)
  \param
*/
VImage VGreyDilation3d(VImage src,VImage se,VImage dest)

{
  int nbands=VImageNBands(src);
  int nrows=VImageNRows(src);
  int ncols=VImageNColumns(src);
  int b,r,c,bb,rr,cc,b0,b1,r0,r1,c0,c1,wnc,wnr,wnb;
  int xsize,ysize,zsize,x,y,z;
  double v,umax,tiny=1.0e-6;
  VRepnKind repn;


  repn = VPixelRepn(src);
  zsize = VImageNBands(se);
  ysize = VImageNRows(se);
  xsize = VImageNColumns(se);

  wnc = xsize / 2;
  wnr = ysize / 2;
  wnb = zsize / 2;

  dest = VSelectDestImage("VGreyDilation3d",dest,nbands,nrows,ncols,repn);
  if (! dest) VError(" err creating dest image");
  VFillImage(dest,VAllBands,0);


  for (b=0; b < nbands; b++) {
    for (r=0; r < nrows; r++) {
      for (c=0; c < ncols; c++) {

	v = VGetPixel(src,b,r,c);
	if (fabs(v) < tiny) continue;
	umax = VPixelMinValue(src);

	b0 = b-wnb;
	if (b0 < 0) b0 = 0;
	b1 = b + wnb;
	if (b1 >= nbands) b1 = nbands-1;
	for (bb=b0; bb<=b1; bb++) {
	  
	  z = bb-b+wnb;

	  r0 = r-wnr;
	  if (r0 < 0) r0 = 0;
	  r1 = r + wnr;
	  if (r1 >= nrows) r1 = nrows-1;
	  for (rr=r0; rr <= r1; rr++) {

	    y = rr-r+wnr;

	    c0 = c-wnc;
	    if (c0 < 0) c0 = 0;
	    c1 = c + wnc;
	    if (c1 >= ncols) c1 = ncols-1;
	    for (cc=c0; cc <= c1; cc++) {

	      x = cc-c+wnc;

	      if (VPixel(se,z,y,x,VBit) > 0) {
		v = VGetPixel(src,bb,rr,cc);
		if (fabs(v) < tiny) continue;
		if (v > umax) umax=v;
	      }
	    }
	  }
	}
	if (repn == VUByteRepn) {
	  umax = (int) (umax + 0.5);
	  if (umax <   0) umax=0;
	  if (umax > 255) umax=255;
	}
	VSetPixel(dest,b,r,c,(VDouble)umax);
      }
    }
  }
  VCopyImageAttrs(src, dest);
  return dest;
}


/*!
  \fn VImage VGreyErosion3d(VImage src,VImage se,VImage dest)
  \brief 3D greylevel morphological erosion
  \param src    input image (any repn)
  \param se     image containing the structural element (bit repn)
  \param dest   output image (any repn)
  \param
*/
VImage VGreyErosion3d(VImage src,VImage se,VImage dest)
{
  int nbands=VImageNBands(src);
  int nrows=VImageNRows(src);
  int ncols=VImageNColumns(src);
  int b,r,c,bb,rr,cc,b0,b1,r0,r1,c0,c1,wnc,wnr,wnb;
  int xsize,ysize,zsize,x,y,z;
  double v,umin,tiny=1.0e-6;
  VRepnKind repn;

  repn  = VPixelRepn(src);
  zsize = VImageNBands(se);
  ysize = VImageNRows(se);
  xsize = VImageNColumns(se);

  wnc = xsize / 2;
  wnr = ysize / 2;
  wnb = zsize / 2;

  dest = VSelectDestImage("VGreyErosion3d",dest,nbands,nrows,ncols,repn);
  if (! dest) VError("err creating dest image");
  VFillImage(dest,VAllBands,0);

  for (b=0; b < nbands; b++) {
    for (r=0; r < nrows; r++) {
      for (c=0; c < ncols; c++) {

	v = VGetPixel(src,b,r,c);	
	if (fabs(v) < tiny) continue;
	umin = VPixelMaxValue(src);

      	b0 = b-wnb;
	if (b0 < 0) b0 = 0;
	b1 = b + wnb;
	if (b1 >= nbands) b1 = nbands-1;
	for (bb=b0; bb <= b1; bb++) {

	  z = bb-b+wnb;

	  r0 = r-wnr;
	  if (r0 < 0) r0 = 0;
	  r1 = r + wnr;
	  if (r1 >= nrows) r1 = nrows-1;
	  for (rr=r0; rr <= r1; rr++) {

	    y = rr-r+wnr;

	    c0 = c-wnc;
	    if (c0 < 0) c0 = 0;
	    c1 = c + wnc;
	    if (c1 >= ncols) c1 = ncols-1;
	    for (cc=c0; cc <= c1; cc++) {

	      x = cc-c+wnc;

	      if (VPixel(se,z,y,x,VBit) > 0) {
		v = VGetPixel(src,bb,rr,cc);	
		if (fabs(v) < tiny) continue;
		if (v < umin) umin=v;
	      }
	    }
	  }
	}
	VSetPixel(dest,b,r,c,(VDouble)umin);
      }
    }
  }

  VCopyImageAttrs(src, dest);
  return dest;
}
