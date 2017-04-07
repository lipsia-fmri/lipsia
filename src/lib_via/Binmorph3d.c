/*! \file 
  3D binary morphology.

This file contains functions for 3D binary morphology:
ersion, dilation, and generation of 3D structuring elements.
  
\par Reference:
  P. Maragos, R.W. Schafer (1990):
  "Morphological Systems for multidimensional signal processing",
  Proc. of the IEEE, Vol. 78, No. 4, pp. 690--709.
  
\par Author:
 Gabriele Lohmann, MPI-CBS
*/

/* From the Vista library: */
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <stdio.h>
#include <stdlib.h>


typedef struct {
  short r;
  short c;
} Pixel, *PixelList;


typedef struct {
  short b;
  short r;
  short c;
} Voxel, *VoxelList;




/*!
  \fn VImage VErodeImage3d(VImage src, VImage dest, VoxelList se, int nse)
  \brief 3D morphological erosion
  \param src   input image (bit repn)
  \param dest  output image (bit repn)
  \param se    structuring element
  \param nse   number of elements in the structuring element
  
 Performs a 3D binary morphological erosion using Minkowski addition.
 The structuring element must be given as an array of voxel-
 addresses, i.e. as SEstruct *list, together with the length
 of that list. A binary raster image can be converted into this structure 
 by calling the function "ConvertSE".
*/
VImage
VErodeImage3d(VImage src, VImage dest, VoxelList se, int nse)
{
  int b,r,c,nbands,nrows,ncols;
  int i,bb,rr,cc;

  nbands = VImageNBands(src);
  nrows  = VImageNRows(src);
  ncols  = VImageNColumns(src);
  
  if (VPixelRepn(src) != VBitRepn) 
    VError("Input image must be of type VBit");

  dest = VSelectDestImage("VErodeImage3d",dest,
                          VImageNBands(src),VImageNRows(src),VImageNColumns(src),
                          VBitRepn);
  VFillImage(dest,VAllBands,0);
  
  for (b=0; b<nbands; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {

	if (VPixel(src,b,r,c,VBit) == 0) {
	  VPixel(dest,b,r,c,VBit) = 0;
	  goto nexte;
	}

	for (i=0; i<nse; i++) {

	  bb = b + se[i].b;
	  if (bb < 0 || bb >= nbands) continue;
	  rr = r + se[i].r;
	  if (rr < 0 || rr >= nrows) continue;
	  cc = c + se[i].c;
	  if (cc < 0 || cc >= ncols) continue;

	  if (VPixel(src,bb,rr,cc,VBit) == 0)
	    goto nexte;
	}

	VPixel(dest,b,r,c,VBit) = 1;
	nexte:
	;
      }
    }
  }
  VCopyImageAttrs (src, dest);
  return dest;
}
		


/*!
  \fn VImage VDilateImage3d(VImage src, VImage dest, VoxelList se, int nse)
  \brief 3D     morphological dilation
  \param src    input image (bit repn)
  \param dest   output image (bit repn)
  \param se     structuring element
  \param nse    number of elements in the structuring element
  
 Performs a 3D binary morphological dilation using Minkowski addition.
 The structuring element must be given as an array of voxel-
 addresses, i.e. as SEstruct *list, together with the length
 of that list. A binary raster image can be converted into this structure 
 by calling the function "ConvertSE".
*/
VImage
VDilateImage3d(VImage src, VImage dest, VoxelList se, int nse)
{
  int b,r,c,nbands,nrows,ncols;
  int bb=0,rr=0,cc=0,i;

  nbands  = VImageNBands(src);
  nrows   = VImageNRows(src);
  ncols   = VImageNColumns(src);


  if (VPixelRepn(src) != VBitRepn) 
    VError("Input image must be of type VBit");

  dest = VSelectDestImage("VDilateImage3d",dest,
                          VImageNBands(src),VImageNRows(src),VImageNColumns(src),
                          VBitRepn);
  VFillImage(dest,VAllBands,0);

  for (b=0; b<nbands; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {

	if (VPixel(src,b,r,c,VBit) == 1) {
	  VPixel(dest,b,r,c,VBit) = 1;
	  goto nextd;
	}

	for (i=0; i<nse; i++) {

	  bb = b + se[i].b;
	  if (bb < 0 || bb >= nbands) continue;
	  rr = r + se[i].r;
	  if (rr < 0 || rr >= nrows) continue;
	  cc = c + se[i].c;
	  if (cc < 0 || cc >= ncols) continue;

	  if (VPixel(src,bb,rr,cc,VBit) == 1) {
	    VPixel(dest,b,r,c,VBit) = 1;
	    goto nextd;
	  }
	}
	nextd:
	;
      }
    }
  }

  VCopyImageAttrs (src, dest);
  return dest;
}



/*!
  \fn VoxelList VConvertSE3d(VImage src,int *nse)
  \brief convert a structuring element given as a raster image into a voxel list.
   The voxel list representation is more efficient computationally.
  \param src   input image (bit repn)
  \param *nse  ptr to the number of elements in the voxel list.
  
  This function generates a structuring element that can be used as input
  into VDilateImage3d or VErodeImage3d.
*/
VoxelList
VConvertSE3d(VImage src,int *nse)
{
  int i,n,b,r,c,nbands,nrows,ncols,npixels;
  int b2,r2,c2;
  VoxelList se;
  VBit *se_pp;

  nbands  = VImageNBands(src);
  nrows   = VImageNRows(src);
  ncols   = VImageNColumns(src);
  npixels = VImageNPixels(src);

  b2 = nbands / 2;
  r2 = nrows  / 2;
  c2 = ncols  / 2;

  se_pp = (VBit *) VPixelPtr(src,0,0,0);
  n = 0;
  for (i=0; i<npixels; i++) 
    if (*se_pp++ > 0) n++;

  se = (VoxelList) VMalloc(sizeof(Voxel) * n);
  if (! se) return NULL;

  i = 0;
  for (b=0; b<nbands; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	if (VPixel(src,b,r,c,VBit) > 0) {
	  se[i].b = b - b2;
	  se[i].r = r - r2;
	  se[i].c = c - c2;
	  i++;
	}
      }
    }
  }
  *nse = n;
  return se;
}
	  


/*!
  \fn VImage VGenSphere3d(VShort radius)
  \brief generate a sphere shaped object. The output can be used as
   a structuring element (after passing it through VConvertSE3d).
  \param radius gives the radius of the SE
*/
VImage
VGenSphere3d(VShort radius)
{
  int b,r,c,bb,rr,cc,dim;
  int d,dmax;
  VImage se;

  if (radius < 0) VError("VGenSphere3d: radius must be larger than 0");

  dmax = radius*radius;
  dim = (int) (2*radius + 1);

  se = VCreateImage(dim,dim,dim,VBitRepn);
  VFillImage(se,VAllBands,0);

  for (b=0; b<dim; b++) {
    for (r=0; r<dim; r++) {
      for (c=0; c<dim; c++) {
	bb = b - radius;
	rr = r - radius;
	cc = c - radius;
	d = bb*bb + rr*rr + cc*cc;
	if (d <= dmax)
	  VPixel(se,b,r,c,VBit) = 1;
      }
    }
  }
  return se;
}



/*!
  \fn VImage VGenSphere2d(VShort radius)
  \brief generate a sphere shaped object, only the center slice is filled with foreground voxels.
   The output can be used as a structuring element (after passing it through VConvertSE3d).
  \param radius gives the radius of the SE
*/
VImage
VGenSphere2d(VShort radius)
{
  int b,r,c,bb,rr,cc,dim;
  int d,dmax;
  VImage se;

  if (radius < 0) VError("VGenSphere2d: radius must be larger than 0");

  dmax = radius*radius;
  dim = (int) (2*radius + 1);

  se = VCreateImage(dim,dim,dim,VBitRepn);
  VFillImage(se,VAllBands,0);

  for (b=0; b<dim; b++) {
    if (b != radius) continue;
    for (r=0; r<dim; r++) {
      for (c=0; c<dim; c++) {
	bb = b - radius;
	rr = r - radius;
	cc = c - radius;
	d = bb*bb + rr*rr + cc*cc;
	if (d <= dmax)
	  VPixel(se,b,r,c,VBit) = 1;
      }
    }
  }
  return se;
}
