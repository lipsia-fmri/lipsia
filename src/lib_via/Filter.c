/*! \file
  Smoothing filters using convolutions

\par Author:
Gabriele Lohmann, MPI-CBS
*/

/* From the Vista library: */
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

/* From the standard C libaray: */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <via/via.h>

#define ABS(x) ((x) < 0 ? -(x) : (x))



/*!
\fn VImage VFilterBox3d (VImage src,VImage dest,int dim)
\param src  input image (any pixel repn)
\param dest output image
\param dim the 1D dimension of the convolution kernel
*/
VImage
VFilterBox3d(VImage src,VImage dest,int dim)
{
  VImage xsrc=NULL,xdest=NULL;
  VImage kernel=NULL;
  VDouble x;

  kernel = VCreateImage(dim,dim,dim,VFloatRepn);
  x = 1.0 / (float) (dim*dim*dim);
  VFillImage(kernel,VAllBands,x);

  /* convert image to float repn */
  xsrc  = VConvertImageCopy(src,NULL,VAllBands,VFloatRepn);
  xdest = VConvolve3d (xsrc,NULL,kernel);
  xdest = VConvertImageCopy(xdest,NULL,VAllBands,VPixelRepn(src));

  VDestroyImage(xsrc);
  VDestroyImage(xdest);

  return dest;
}



/* Gaussian function */
double 
xxgauss(double x,double sigma)
{
  double y,z,a=2.506628273;
  z = x / sigma;
  y = exp((double)-z*z*0.5)/(sigma * a);
  return y;
}


VImage
VGaussKernel(double sigma)
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
    u = xxgauss(x,sigma);
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



/*!
\fn VImage VFilterGauss2d (VImage src,VImage dest,double sigma)
\brief 2d Gauss filter
\param src    input image (any pixel repn)
\param dest   output image
\param double sigma
*/
VImage
VFilterGauss2d (VImage src,VImage dest,double sigma)
{
  VImage xsrc=NULL,xdest=NULL,kernel=NULL;

  if (sigma <= 0) VError("VFilterGauss2d: sigma must be positive");
  kernel = VGaussKernel(sigma);

  /* convert image to float repn */
  xsrc  = VConvertImageCopy(src,xsrc,VAllBands,VFloatRepn);

  xdest = VConvolveCol(xsrc,NULL,kernel);
  xdest = VConvolveRow(xdest,NULL,kernel);

  VDestroyImage(xsrc);

  dest = VConvertImageCopy(xdest,dest,VAllBands,VPixelRepn(src));
  VDestroyImage(xdest);

  return dest;
}


/*!
\fn VImage VFilterGauss3d (VImage src,VImage dest,double sigma)
\brief 3d Gauss filter
\param src    input image (any pixel repn)
\param dest   output image
\param double sigma
*/
VImage
VFilterGauss3d (VImage src,VImage dest,double sigma)
{
  VImage xsrc=NULL,xdest=NULL,tmp=NULL,kernel=NULL;

  if (sigma <= 0) VError("VFilterGauss3d: sigma must be positive");
  kernel = VGaussKernel(sigma);

  /* convert image to float repn */
  xsrc  = VConvertImageCopy(src,NULL,VAllBands,VFloatRepn);

  xdest = VConvolveCol(xsrc,xdest,kernel);
  tmp   = VConvolveRow(xdest,tmp,kernel);
  xdest = VConvolveBand(tmp,xdest,kernel);

  VDestroyImage(xsrc);

  dest = VConvertImageCopy(xdest,NULL,VAllBands,VPixelRepn(src));
  VDestroyImage(xdest);

  return dest;
}



  

