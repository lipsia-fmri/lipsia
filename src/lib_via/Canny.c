/*! \file
  3D Canny edge detection

\par Reference:
J. Canny (1986). "A computational approach to edge detection",
IEEE-PAMI, Vol.8, No.6. pp.679--698.

\par Author:
Gabriele Lohmann, MPI-CBS
*/


/*
** 3D Canny edge detection
**
** G.Lohmann, 1996
*/

/* From the Vista library: */
#include <viaio/Vlib.h>
#include <viaio/VImage.h>


/* From the standard C libaray: */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <via/via.h>


/*
** derivative of a gaussian
*/
void
VDerivGaussian(int dim,VImage *gkernel,VImage *dkernel)
{  
  int i,d,n;
  double x,y;
  double pi = 3.14159265,sigma;
  double norm,u;
  double *gauss,*deriv;

  n = dim-1;
  d = dim / 2;
  sigma = 0.25 * n;


  gauss = (double *) VMalloc(sizeof(double) * dim);
  deriv = (double *) VMalloc(sizeof(double) * dim);

  norm = sqrt((double)2.0 * pi) * sigma;
  for (i=0; i<=n; i++) {
    x = (double) (i - d);
    u = - (x * x) / (2.0 * sigma * sigma);
    y = exp((double) u) / norm;
    gauss[i] = y;
  }

  for (i=0; i<=n; i++) {
    x = (double) (i - d);
    y = - x * gauss[i] / (sigma * sigma);
    deriv[i] = y;
  }

  for (i=0; i<dim; i++) {
    VPixel((*gkernel),0,0,i,VFloat) = gauss[i];
    VPixel((*dkernel),0,0,i,VFloat) = deriv[i];
  }
}


/*
** 3D Canny filter in column direction
*/
VImage
VCanny3d_col (VImage src,VImage dest,VImage tmp,VImage gkernel,VImage dkernel)
{
  dest = VConvolveCol(src,NULL,dkernel);
  tmp  = VConvolveRow(dest,tmp,gkernel);
  dest = VConvolveBand(tmp,dest,gkernel);
  return dest;
}

/*
** 3D canny filter (row)
*/
VImage
VCanny3d_row (VImage src,VImage dest,VImage tmp,VImage gkernel,VImage dkernel)
{
  dest = VConvolveCol(src,NULL,gkernel);
  tmp  = VConvolveRow(dest,tmp,dkernel);
  dest = VConvolveBand(tmp,dest,gkernel);
  return dest;
}

/*
** 3D canny filter (band)
*/
VImage
VCanny3d_band (VImage src,VImage dest,VImage tmp,VImage gkernel,VImage dkernel)
{
  dest = VConvolveCol(src,NULL,gkernel);
  tmp  = VConvolveRow(dest,tmp,gkernel);
  dest = VConvolveBand(tmp,dest,dkernel);
  return dest;
}


/*
** 2D Canny filter in column direction
*/
VImage
VCanny2d_col (VImage src,VImage dest,VImage tmp,VImage gkernel,VImage dkernel)
{
  tmp   = VConvolveCol(src,tmp,dkernel);
  dest  = VConvolveRow(tmp,dest,gkernel);
  return dest;
}

/*
** 2D canny filter (row)
*/
VImage
VCanny2d_row (VImage src,VImage dest,VImage tmp,VImage gkernel,VImage dkernel)
{
  tmp   = VConvolveCol(src,tmp,gkernel);
  dest  = VConvolveRow(tmp,dest,dkernel);
  return dest;
}


/*!
\fn void VCanny3d (VImage src,int wsize,VImage *gradb,VImage *gradr,VImage *gradc)
\param src  input image
\param *gradb  output gradient in slice direction (float repn)
\param *gradr  output gradient in row direction (float repn)
\param *gradc  output gradient in column direction (float repn)
\param wsize   dimension of convolution kernel
*/
void
VCanny3d(VImage src,int dim,VImage *gradb,VImage *gradr,VImage *gradc)
{
  VImage xsrc=NULL,tmp=NULL;
  VImage gkernel=NULL,dkernel=NULL;

  /* get convolution kernels */
  gkernel = VCreateImage(1,1,dim,VFloatRepn);
  dkernel = VCreateImage(1,1,dim,VFloatRepn);
  VDerivGaussian(dim,&gkernel,&dkernel);


  /* convert image to float repn */
  xsrc = VConvertImageCopy(src,NULL,VAllBands,VFloatRepn);
  tmp = VCopyImage(xsrc,NULL,VAllBands);

  /* 1D convolutions */
  *gradc = VCanny3d_col(xsrc,NULL,tmp,gkernel,dkernel);
  *gradr = VCanny3d_row(xsrc,NULL,tmp,gkernel,dkernel);
  *gradb = VCanny3d_band(xsrc,NULL,tmp,gkernel,dkernel);

  VCopyImageAttrs (src, (*gradc));
  VCopyImageAttrs (src, (*gradr));
  VCopyImageAttrs (src, (*gradb));

  VDestroyImage(tmp);
  VDestroyImage(xsrc);
  VDestroyImage(dkernel);
  VDestroyImage(gkernel);
}



/*!
\fn void VCanny2d (VImage src,int wsize,VImage *gradr,VImage *gradc)
\param src  input image
\param *gradr  output gradient in row direction (float repn)
\param *gradc  output gradient in column direction (float repn)
\param wsize   dimension of convolution kernel
*/
void
VCanny2d(VImage src,int dim,VImage *gradr,VImage *gradc)
{
  VImage xsrc=NULL,tmp=NULL;
  VImage gkernel=NULL,dkernel=NULL;

  /* get convolution kernels */
  gkernel = VCreateImage(1,1,dim,VFloatRepn);
  dkernel = VCreateImage(1,1,dim,VFloatRepn);
  VDerivGaussian(dim,&gkernel,&dkernel);


  /* convert image to float repn */
  xsrc = VConvertImageCopy(src,NULL,VAllBands,VFloatRepn);
  tmp = VCopyImage(xsrc,NULL,VAllBands);

  /* 1D convolutions */
  *gradc = VCanny2d_col(xsrc,NULL,tmp,gkernel,dkernel);
  *gradr = VCanny2d_row(xsrc,NULL,tmp,gkernel,dkernel);

  VCopyImageAttrs (src, (*gradc));
  VCopyImageAttrs (src, (*gradr));
  
  VDestroyImage(tmp);
  VDestroyImage(xsrc);
  VDestroyImage(dkernel);
  VDestroyImage(gkernel);
}








