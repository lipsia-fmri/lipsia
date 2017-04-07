/*! \file
  3D grey level morphology.

\par Reference:
  P. Maragos, R.W. Schafer (1990):
  "Morphological Systems for multidimensional signal processing",
  Proc. of the IEEE, Vol. 78, No. 4, pp. 690--709.

\par Author:
 Gabriele Lohmann, MPI-CBS
*/

#include <viaio/Vlib.h>
#include <via.h>
#include <stdio.h>


/*!
  \fn VImage VGreyDilation3d(VImage src,VImage se,VImage dest)
  \brief 3D greylevel morphological dilation
  \param src    input image (any repn)
  \param se     raster image containing the structural element (bit repn)
  \param dest   output image (any repn)
  \param
*/
VImage
VGreyDilation3d(VImage src,VImage se,VImage dest)

{
  int nbands=VImageNBands(src), 
    nrows=VImageNRows(src), 
    ncols=VImageNColumns(src);
  int b,r,c,bb,rr,cc,b0,b1,r0,r1,c0,c1,wnc,wnr,wnb;
  int xsize,ysize,zsize,x,y,z;
  double v,umax;
  VRepnKind repn;
  int background=0;

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

	v = VReadPixel(src,b,r,c);
	if ((int)v == background && repn == VUByteRepn) continue;

	umax=VPixelMinValue(src);

	b0 = (b > wnb) ? b0 = b-wnb : 0;
	b1 = (b < nbands - wnb) ? b1 = b + wnb : nbands - 1;
	for (bb=b0; bb <= b1; bb++) {
	  
	  z = bb-b+wnb;

	  r0 = (r > wnr) ? r0 = r-wnr : 0;
	  r1 = (r < nrows - wnr) ? r1 = r + wnr : nrows - 1;
	  for (rr=r0; rr <= r1; rr++) {

	    y = rr-r+wnr;

	    c0 = (c > wnc) ? c0 = c-wnc : 0;
	    c1 = (c < ncols - wnc) ? c1 = c + wnc : ncols - 1;
	    for (cc=c0; cc <= c1; cc++) {

	      x = cc-c+wnc;

	      if (VPixel(se,z,y,x,VBit) > 0) {
		v = VReadPixel(src,bb,rr,cc);
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


  /* Let the destination inherit any attributes of the source image: */
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
VImage
VGreyErosion3d(VImage src,VImage se,VImage dest)
{
  int nbands=VImageNBands(src), 
    nrows=VImageNRows(src), 
    ncols=VImageNColumns(src);
  int b,r,c,bb,rr,cc,b0,b1,r0,r1,c0,c1,wnc,wnr,wnb;
  int xsize,ysize,zsize,x,y,z;
  double v,umin;
  VRepnKind repn;
  int background=0;


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

	v = VReadPixel(src,b,r,c);
	if ((int)v == background && repn == VUByteRepn) continue;

	umin=VPixelMaxValue(src);
      
	b0 = (b > wnb) ? b0 = b-wnb : 0;
	b1 = (b < nbands - wnb) ? b1 = b + wnb : nbands - 1;
	for (bb=b0; bb <= b1; bb++) {

	  z = bb-b+wnb;

	  r0 = (r > wnr) ? r0 = r-wnr : 0;
	  r1 = (r < nrows - wnr) ? r1 = r + wnr : nrows - 1;
	  for (rr=r0; rr <= r1; rr++) {

	    y = rr-r+wnr;

	    c0 = (c > wnc) ? c0 = c-wnc : 0;
	    c1 = (c < ncols - wnc) ? c1 = c + wnc : ncols - 1;
	    for (cc=c0; cc <= c1; cc++) {

	      x = cc-c+wnc;

	      if (VPixel(se,z,y,x,VBit) > 0) {
		v = VReadPixel(src,bb,rr,cc);
		if ((int)v != background || repn != VUByteRepn)
		  if (v < umin) umin=v;
	      }
	    }
	  }
	}
	if (repn == VUByteRepn) {
	  umin = (int) (umin+0.5);
	  if (umin < 0) umin=0;
	  if (umin > 255) umin=255;
	}
	VSetPixel(dest,b,r,c,(VDouble)umin);
      }
    }
  }

  /* Let the destination inherit any attributes of the source image: */
  VCopyImageAttrs(src, dest);
  return dest;
}
