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

#define ABS(x) ((x) < 0 ? -(x) : (x))

extern VImage VConvolveCol(VImage,VImage,float *,int);
extern VImage VConvolveRow(VImage,VImage,float *,int);
extern VImage VConvolveBand(VImage,VImage,float *,int);


/* Gaussian function */
double  xxgauss(double x,double sigma)
{
  double y,z,a=2.506628273;
  z = x / sigma;
  y = exp((double)-z*z*0.5)/(sigma * a);
  return y;
}


float *VGaussKernel(double sigma,int *xdim)
{
  int    i,dim,n;
  double x,u,step;
  double sum;

  dim  = 3.0 * sigma + 1;
  n    = 2*dim+1;
  step = 1;

  float *kernel = (float *) VCalloc(n,sizeof(float));
  
  sum = 0;
  x = -(float)dim;
  for (i=0; i<n; i++) {
    u = xxgauss(x,sigma);
    sum += u;
    kernel[i] = u;
    x += step;
  }

  /* normalize */
  for (i=0; i<n; i++) {
    kernel[i] /= sum;
  }
  *xdim = n;
  return kernel;
}

  

VImage VFilterGauss3d (VImage src,VImage dest,double *sigma)
{
  VImage xsrc=NULL,xdest=NULL,tmp=NULL;
  int dim_col,dim_row,dim_slice;
  
  float *kernel_col = VGaussKernel(sigma[0],&dim_col);
  float *kernel_row = VGaussKernel(sigma[1],&dim_row);
  float *kernel_slice = VGaussKernel(sigma[2],&dim_slice);

  /* convert to float if needed */
  if (VPixelRepn(src) == VFloatRepn) {
    xsrc = src;
  }
  else {
    xsrc = VConvertImageCopy(src,NULL,VAllBands,VFloatRepn);
  }
  

  /* separable filter */
  xdest = VConvolveCol(xsrc,xdest,kernel_col,dim_col);
  tmp   = VConvolveRow(xdest,tmp,kernel_row,dim_row);
  xdest = VConvolveBand(tmp,xdest,kernel_slice,dim_slice);

  VFree(kernel_col);
  VFree(kernel_row);
  VFree(kernel_slice);


  /* convert back to original pixel repn */
  if (VPixelRepn(src) == VFloatRepn) {
    return xdest;
  }
  else {
    VDestroyImage(xsrc);
    dest = VConvertImageCopy(xdest,NULL,VAllBands,VPixelRepn(src));
    VDestroyImage(xdest);
    return dest;
  }
  return NULL;
}
