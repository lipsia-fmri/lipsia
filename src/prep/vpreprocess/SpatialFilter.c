/****************************************************************
 *
 * vpreprocess: SpatialFilter.c
 *
 * Copyright (C) Max Planck Institute 
 * for Human Cognitive and Brain Sciences, Leipzig
 *
 * <lipsia@cbs.mpg.de>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * 
 * $Id: SpatialFilter.c 3190 2008-04-01 16:06:57Z karstenm $
 *
 *****************************************************************/

#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define NSLICES 10000
#define ABS(x) ((x) < 0 ? -(x) : (x))


/* Gaussian function */
double xsgauss(double x,double sigma)
{
  double y,z,a=2.506628273;
  z = x / sigma;
  y = exp(-z*z*0.5)/(sigma * a);
  return y;
}


VImage VSGaussKernel(double sigma)
{
  int    i,dim,n;
  double x,u,step;
  VImage kernel=NULL;
  double sum;

  dim  = 3.0 * sigma + 1;
  n    = 2*dim+1;
  step = 1;

  kernel = VCreateImage(1,1,n,VFloatRepn);

  sum = 0;
  x = -(float)dim;
  for (i=0; i<n; i++) {
    u = xsgauss(x,sigma);
    sum += u;
    VPixel(kernel,0,0,i,VFloat) = u;
    x += step;
  }

  /* normalize */
  for (i=0; i<n; i++) {
    u = VPixel(kernel,0,0,i,VFloat);
    u /= sum;
    VPixel(kernel,0,0,i,VFloat) = u;
  }

  return kernel;
}



VImage VSConvolveCol (VImage src,VImage dest,VImage kernel)
{
  int b,r,c,nbands,nrows,ncols;
  int c0,c1,cc;
  float sum,x;
  VFloat *float_pp;
  int dim,d;

  nrows  = VImageNRows (src);
  ncols  = VImageNColumns (src);
  nbands = VImageNBands (src);
  dest = VCopyImage(src,dest,VAllBands);

  dim  = VImageNColumns(kernel);
  d    = dim/2;

  for (b=0; b<nbands; b++) {
    for (r=0; r<nrows; r++) {
      for (c=d; c<ncols-d; c++) {

	float_pp = (VFloat *) VImageData(kernel);
	sum = 0;
	c0 = c-d;
	c1 = c+d;
	if (c0 < 0) c0 = 0;
	if (c1 >= ncols) c1 = ncols-1;
	for (cc=c0; cc<=c1; cc++) {
	  x = VPixel(src,b,r,cc,VFloat);
	  sum += x * (*float_pp++);
	}
	VPixel(dest,b,r,c,VFloat) = sum;
      }
    }
  }
  return dest;
}


VImage VSConvolveRow (VImage src,VImage dest,VImage kernel)
{
  int b,r,c,nbands,nrows,ncols;
  int r0,r1,rr;
  float sum,x;
  VFloat *float_pp;
  int d,dim;

  dim = VImageNColumns(kernel);
  d = dim/2;

  nrows  = VImageNRows (src);
  ncols  = VImageNColumns (src);
  nbands = VImageNBands (src);
  dest = VCopyImage(src,dest,VAllBands);

  for (b=0; b<nbands; b++) {
    for (r=d; r<nrows-d; r++) {
      for (c=0; c<ncols; c++) {

	float_pp = (VFloat *) VImageData(kernel);
	sum = 0;
	r0 = r-d;
	r1 = r+d;
	if (r0 < 0) r0 = 0;
	if (r1 >= nrows) r1 = nrows-1;
	for (rr=r0; rr<=r1; rr++) {
	  x = VPixel(src,b,rr,c,VFloat);
	  sum += x * (*float_pp++);
	}
	VPixel(dest,b,r,c,VFloat) = sum;
      }
    }
  }
  return dest;
}


VImage VSConvolveBand (VImage src,VImage dest,VImage kernel)
{
  int b,r,c,nbands,nrows,ncols;
  int b0,b1,bb;
  float sum,x;
  VFloat *float_pp;
  int d,dim;

  dim = VImageNColumns(kernel);
  d = dim/2;

  nrows  = VImageNRows (src);
  ncols  = VImageNColumns (src);
  nbands = VImageNBands (src);
  dest = VCopyImage(src,dest,VAllBands);

  for (b=d; b<nbands-d; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {

	float_pp = (VFloat *) VImageData(kernel);
	sum = 0;
	b0 = b-d;
	b1 = b+d;
	if (b0 < 0) b0 = 0;
	if (b1 >= nbands) b1 = nbands-1;
	for (bb=b0; bb<=b1; bb++) {
	  x = VPixel(src,bb,r,c,VFloat);
	  sum += x * (*float_pp++);
	}
	VPixel(dest,b,r,c,VFloat) = sum;
      }
    }
  }
  return dest;
}



VImage VGauss3d(VImage src,VImage dest,VImage kernel)
{
  static VImage tmp=NULL;

  dest = VSConvolveCol(src,dest,kernel);
  tmp  = VSConvolveRow(dest,tmp,kernel);
  dest = VSConvolveBand(tmp,dest,kernel);
  return dest;
}


void VSpatialFilter(VAttrList list,VDouble fwhm)
{
  VAttrListPosn posn;
  VImage src[NSLICES],xsrc=NULL,tmp=NULL,dest=NULL,kernel=NULL,tmp2d=NULL;
  VString str=NULL;
  float v0,v1,v2,v3;
  int b,r,c,i,size;
  int n,nslices,nrows,ncols,dim;
  double u, sigma=0;
  extern VImage VGaussianConv (VImage,VImage,VBand,double,int);


  /* get image dimensions */
  dim = 3;
  v0 = v1 = v2 = v3 = 1;
  n = i = nrows = ncols = 0;
  str = VMalloc(100);
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (i >= NSLICES) VError(" too many slices");
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    VGetAttrValue (& posn, NULL,VImageRepn, & xsrc);
    if (VPixelRepn(xsrc) != VShortRepn) continue;
    if (VImageNBands(xsrc) > n) n = VImageNBands(xsrc);
    if (VImageNRows(xsrc) > nrows) nrows = VImageNRows(xsrc);
    if (VImageNColumns(xsrc) > ncols) ncols = VImageNColumns(xsrc);
    if (VGetAttr (VImageAttrList (xsrc), "voxel", NULL,VStringRepn, (VPointer) & str) == VAttrFound) {
      sscanf(str,"%f %f %f",&v1,&v2,&v3);
    }
    src[i] = xsrc;
    i++;    
  }
  nslices = i;


  /* in general, apply 3D spatial filtering */
  dim = 3;

  /* only if clearly non-isotropic apply 2D filtering */
  v0 = 0.5*(v1+v2);
  if (ABS(v0-v3) > 0.5) {
    VWarning(" non-isotropic voxel grid, apply 2D filtering");
    dim = 2;
  }


  /*
  ** Sigma
  */
  sigma  = fwhm/sqrt(8.0*log(2.0));
  sigma /= (double)v1;


  /*
  ** 2D gauss filtering
  */
  if (dim==2) {
    fprintf(stderr," 2D spatial filter: fwhm=  %.3f mm sigma= %.3f vox\n",fwhm,sigma); 

    size = (int)(6.0 * sigma + 1.5);
    if ((size & 1) == 0) size++;
    fprintf(stderr," size= %d\n",size);

    for (b=0; b<nslices; b++) {
      if (VImageNRows(src[b]) < 2) continue;
      tmp2d  = VGaussianConv(src[b],tmp2d, VAllBands, sigma, size);
      src[b] = VCopyImagePixels(tmp2d,src[b],VAllBands);
    }
    VDestroyImage(tmp2d);
  }


  /*
  ** 3D gauss filtering
  */
  if (dim==3) {
    fprintf(stderr," 3D spatial filter: fwhm=  %.3f mm\n",fwhm); 

    kernel = VSGaussKernel(sigma);

    tmp = VCreateImage(nslices,nrows,ncols,VFloatRepn);
    VFillImage(tmp,VAllBands,0);
    
    for (i=0; i<n; i++) {
      if (i%20 == 0) fprintf(stderr," t= %5d of %d\r",i,n);
      
      VFillImage(tmp,VAllBands,0);
      for (b=0; b<nslices; b++) {
	if (VImageNRows(src[b]) < 2) continue;
	for (r=0; r<nrows; r++) {
	  for (c=0; c<ncols; c++) {
	    u = VPixel(src[b],i,r,c,VShort);
	    VPixel(tmp,b,r,c,VFloat) = u;
	  }
	}
      }
      
      dest = VGauss3d (tmp,dest,kernel);
      
      for (b=0; b<nslices; b++) {
	if (VImageNRows(src[b]) < 2) continue;
	for (r=0; r<nrows; r++) {
	  for (c=0; c<ncols; c++) {
	    u = VPixel(dest,b,r,c,VFloat);
	    VPixel(src[b],i,r,c,VShort) = u;
	  }
	}
      }
    }
  }
}
