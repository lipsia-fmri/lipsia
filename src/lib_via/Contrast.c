/*! \file
  
Contrast enhancement.


A piecewise-linear contrast enhancement is performed.
The contrast stretching function is derived from the mean
value and the standard deviation of the input image.
The output pixel repn need not be identical to the input 
image. Thus, it is possible to convert to ubyte repn
for easier visualization.


\par Author:
Gabriele Lohmann, MPI-CBS
*/


/* From the Vista library: */
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <via/via.h>

#define ABS(x) ((x) > 0 ? (x) : -(x))

/*!
\fn VImage VContrast(VImage src,VImage dest,VRepnKind repn,VFloat alpha,VFloat background)
\param src   input image  (any repn)
\param dest  output image (any repn)
\param repn  the output pixel repn (e.g. VUByteRepn)
\param alpha contrast stretching factor. 
The function stretches grey values between mean-alpha*sigma and mean+alpha*sigma.
\param background input grey values with absolute values less than <background> are
assumed to be image background and are not used for computing the image mean and sigma.
If set to zero, it has no effect.
*/
VImage
VContrast(VImage src,VImage dest,VRepnKind repn,VFloat alpha,VFloat background)
{
  int i,nbands,nrows,ncols,npixels;
  float u,xmin,xmax,ymin,ymax,slope,smin,smax;
  float sum1,sum2,nx;
  float mean,sigma;

  if (repn == VDoubleRepn || repn == VLongRepn || repn == VSByteRepn)
    VError(" double, long and sbyte are not supported by VContrast");

  if (alpha <= 0) VError("alpha must be positive");


  nbands = VImageNBands(src);
  nrows  = VImageNRows(src);
  ncols  = VImageNColumns(src);
  npixels = VImageNPixels(src);

  dest = VSelectDestImage("VContrast",dest,nbands,nrows,ncols,repn);
  if (! dest) VError(" err creating dest image");
  VFillImage(dest,VAllBands,0);

  smin = VPixelMaxValue(src);
  smax = VPixelMinValue(src);

  sum1 = sum2 = nx = 0;
  for (i=0; i<npixels; i++) {
    u = VGetPixelValue (src,i);
    if (ABS(u) < background) continue;
    if (u < smin) smin = u;
    if (u > smax) smax = u;
    sum1 += u;
    sum2 += u*u;
    nx++;
  }
  if (nx < 2) VError(" no foreground pixels found");

  mean  = sum1/nx;
  sigma = sqrt((double)((sum2 - nx * mean * mean) / (nx - 1.0)));

  ymax = VPixelMaxValue(dest);
  ymin = VPixelMinValue(dest);

  xmax = mean + alpha * sigma;
  if (xmax > smax) xmax = smax;

  xmin = mean - alpha * sigma;
  if (xmin < smin) xmin = smin;

  slope = (ymax - ymin) / (xmax - xmin);

  for (i=0; i<npixels; i++) {
    u = VGetPixelValue (src,i);
    if (u != 0) u = slope * (u - xmin);
    if (u > ymax) u = ymax;
    if (u < ymin) u = ymin;
    VSetPixelValue(dest,i,(VFloat) u);
  }

  VCopyImageAttrs (src, dest);
  return dest;
}



/*!
\fn VImage VContrastUByte(VImage src,VImage dest,VFloat percent,VFloat background)
\param src   input image  (ubyte repn)
\param dest  output image  (ubyte repn)
\param percent percentage of pixels to ignore at either end of the histogram.
\param background input grey values with absolute values less than <background> are
assumed to be image background.
*/
VImage
VContrastUByte(VImage src,VImage dest,VFloat low,VFloat high)
{
  int nbands,nrows,ncols,npixels;
  float u,v,xmin,xmax,slope,sum,background;
  float histo[256];
  int i,j;
  VUByte *src_pp,*dest_pp;

  if (VPixelRepn(src) != VUByteRepn) VError(" input pixel repn must be ubyte");

  nbands = VImageNBands(src);
  nrows  = VImageNRows(src);
  ncols  = VImageNColumns(src);
  npixels = nbands * nrows * ncols;

  dest = VSelectDestImage("VContrastUByte",dest,nbands,nrows,ncols,VUByteRepn);
  if (! dest) VError(" err creating dest image");
  VFillImage(dest,VAllBands,0);

  for (j=0; j<256; j++) histo[j] = 0;

  xmin = VPixelMaxValue(src);
  xmax = VPixelMinValue(src);
  background = xmin;
  
  src_pp = (VUByte *) VImageData(src);
  for (i=0; i<npixels; i++) {
    j = *src_pp++;
    if (j <= background) continue;
    histo[j]++;
  }
 
  sum = 0;
  for (j=0; j<256; j++) sum += histo[j];
  for (j=0; j<256; j++) histo[j] /= sum;

  xmin = 0;
  sum  = 0;
  for (j=0; j<256; j++) {
    if (j > background) sum += histo[j];
    if (sum > low) break;
  }
  xmin = j;

  xmax = 255.0;
  sum = 0;
  for (j=255; j>0; j--) {
    if (j > background) sum += histo[j];
    if (sum > high) break;
  }
  xmax = j;


  slope = (float) (255.0) / ((float) (xmax - xmin));
  
  src_pp  = (VUByte *) VImageData(src);
  dest_pp = (VUByte *) VImageData(dest);

  for (i=0; i<npixels; i++) {
    u = *src_pp;
 
    if (u <= background) {
      v = 0;
    }
    else {
      v = (int) (slope * (u - xmin) + 0.5);
      if (v <   0) v = 0;
      if (v > 255) v = 255;
    }

    *dest_pp = (VUByte) v;
    src_pp++;
    dest_pp++;
  }

  VFree(histo);
  VCopyImageAttrs (src, dest);
  return dest;
}



/*!
\fn VImage VContrastShort(VImage src,VImage dest,VFloat percent,VFloat background)
\param src   input image  (short repn)
\param dest  output image  (ubyte repn)
*/
VImage
VContrastShort(VImage src,VImage dest,VFloat low,VFloat high)
{
  int nbands,nrows,ncols,npixels;
  float u,v,xmin,xmax,slope,sum;
  float *histo;
  int i,j,dim;
  VShort *src_pp;
  VUByte *dest_pp;
  double smin,smax;
  double percent1,percent2;

  if (VPixelRepn(src) != VShortRepn) VError(" input pixel repn must be short");

  nbands = VImageNBands(src);
  nrows  = VImageNRows(src);
  ncols  = VImageNColumns(src);
  npixels = nbands * nrows * ncols;

  dest = VSelectDestImage("VContrastShort",dest,nbands,nrows,ncols,VUByteRepn);
  if (! dest) VError(" err creating dest image");
  VFillImage(dest,VAllBands,0);


  smin = VRepnMinValue(VShortRepn);
  smax = VRepnMaxValue(VShortRepn);

  percent1 = low;    /* unten  */
  percent2 = high;   /* oben   */


  dim = 2.0 * smax + 1.0;
  histo = (float *) VCalloc(dim,sizeof(float));
  for (j=0; j<dim; j++) histo[j] = 0;


  src_pp = (VShort *) VImageData(src);
  for (i=0; i<npixels; i++) {
    j = *src_pp++;
    j -= smin;
    histo[j]++;
  }
 
  sum = 0;
  for (j=0; j<dim; j++) sum += histo[j];
  for (j=0; j<dim; j++) histo[j] /= sum;

  xmin = 0;
  sum  = 0;
  for (j=0; j<dim; j++) {
    sum += histo[j];
    if (sum > percent1) break;
  }
  xmin = j+smin;

  xmax = dim;
  sum = 0;
  for (j=dim; j>0; j--) {
    sum += histo[j];
    if (sum > percent2) break;
  }
  xmax = j+smin;


  slope = 255.0f / (xmax - xmin);
  
  src_pp  = (VShort *) VImageData(src);
  dest_pp = (VUByte *) VImageData(dest);
  for (i=0; i<npixels; i++) {
    u = *src_pp++;
    v = (int) (slope * (u - xmin) + 0.5);
    if (v < 0) v = 0;
    if (v > 255) v = 255;
    *dest_pp++ = (VUByte) v;
  }

  VFree(histo);
  VCopyImageAttrs (src, dest);
  return dest;
}



VImage
VHistoEqualize(VImage src,VImage dest,VFloat exponent)
{
  int nbands,nrows,ncols,npixels;
  float u,v,sum;
  float *histo,*p;
  int i,j,dim;
  VShort *src_pp;
  VUByte *dest_pp;
  double smin,smax,x,y;


  if (VPixelRepn(src) != VShortRepn) VError(" input pixel repn must be short");
  if (exponent < 0.5) VError("parameter '-exponent' should be >= 0.5"); 
  if (exponent > 10) VWarning("parameter '-exponent' should be < 10"); 

  nbands = VImageNBands(src);
  nrows  = VImageNRows(src);
  ncols  = VImageNColumns(src);
  npixels = nbands * nrows * ncols;

  dest = VSelectDestImage("VContrastShort",dest,nbands,nrows,ncols,VUByteRepn);
  if (! dest) VError(" err creating dest image");
  VFillImage(dest,VAllBands,0);


  y = (double) exponent;
  smin = VRepnMinValue(VShortRepn);
  smax = VRepnMaxValue(VShortRepn);

  dim = 2.0 * smax + 1.0;
  histo = (float *) VCalloc(dim,sizeof(float));
  for (j=0; j<dim; j++) histo[j] = 0;

  p = (float *) VCalloc(dim,sizeof(float));


  src_pp = (VShort *) VImageData(src);
  for (i=0; i<npixels; i++) {
    j = *src_pp++;
    j -= smin;
    if (j == 0) continue;
    histo[j]++;
  }
 
  sum = 0;
  for (j=0; j<dim; j++) sum += histo[j];
  for (j=0; j<dim; j++) histo[j] /= sum;


  /* cumulative hist */
  for (i=0; i<dim; i++) {
    sum = 0;
    for (j=0; j<=i; j++) sum += histo[j];
    p[i] = sum;
  }


  /* make lut */
  for (i=0; i<dim; i++) {
    x = (double)p[i];
    if (x > 0)
      p[i] = (float)(pow(x,y) * 255.0);
  }

  /* apply lut */
  src_pp  = (VShort *) VImageData(src);
  dest_pp = (VUByte *) VImageData(dest);
  for (i=0; i<npixels; i++) {
    u = *src_pp++;
    j = (int) (u-smin);
    v = (double)p[j];
    if (v < 0) v = 0;
    if (v > 255) v = 255;
    *dest_pp++ = (VUByte) v;
  }

  VFree(p);
  VFree(histo);
  VCopyImageAttrs (src, dest);
  return dest;
}




/*!
\fn VImage VMapImageRange(VImage src,VImage dest,VRepnKind repn)
\brief
The minimum and maximum grey values of the input images are computed,
and a linear mapping is performed mapping the input minimum(maximum) 
of the input image to the min(max) value of the output pixel repn.
E.g. if the input image min is -17 and its max is +2376, and the output repn
is VUByteRepn, then the linear mapping function maps -17 to 0 and
+2376 to 255.
\param src   input image  (any repn)
\param dest  output image
\param repn  the output pixel repn (e.g. VUByteRepn)
*/
VImage
VMapImageRange(VImage src,VImage dest,VRepnKind repn)
{
  int i,nbands,nrows,ncols,npixels;
  float u,ymin,ymax,dmin,dmax,slope;

  if (repn == VDoubleRepn || repn == VLongRepn || repn == VSByteRepn)
    VError(" VMapImageRange: double, long and sbyte are not supported");

  nbands = VImageNBands(src);
  nrows  = VImageNRows(src);
  ncols  = VImageNColumns(src);
  npixels = VImageNPixels(src);

  dest = VSelectDestImage("VMapImageRange",dest,nbands,nrows,ncols,repn);
  if (! dest) VError(" err creating dest image");
  VFillImage(dest,VAllBands,0);


  ymin = VPixelMaxValue(src);
  ymax = VPixelMinValue(src);

  for (i=0; i<npixels; i++) {
    u = VGetPixelValue (src,i);
    if (u < ymin) ymin = u;
    if (u > ymax) ymax = u;
  }

  dmax = VPixelMaxValue(dest);
  dmin = VPixelMinValue(dest);

  slope = (dmax - dmin) / (ymax - ymin);

  for (i=0; i<npixels; i++) {
    u = VGetPixelValue (src,i);
    u = slope * (u - ymin);
    if (u < dmin) u = dmin;
    if (u > dmax) u = dmax;
    VSetPixelValue(dest,i,(VFloat) u);
  }

  VCopyImageAttrs (src, dest);
  return dest;
}




VImage
VContrastAny(VImage src,VImage dest,VFloat low,VFloat high)
{
  int nbands,nrows,ncols;
  double xmin,xmax,slope,sum,a;
  float *histo;
  int b,r,c,j,dim;
  VUByte *dest_pp;
  double smin,smax,u,v,tiny;

  nbands = VImageNBands(src);
  nrows  = VImageNRows(src);
  ncols  = VImageNColumns(src);

  dest = VSelectDestImage("VContrastAny",dest,nbands,nrows,ncols,VUByteRepn);
  if (! dest) VError(" err creating dest image");
  VFillImage(dest,VAllBands,0);

  
  smin = VPixelMaxValue(src);
  smax = VPixelMinValue(src);

  for (b=0; b<nbands; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	u = VGetPixel(src,b,r,c);
	if (u < smin) smin = u;
	if (u > smax) smax = u;
      }
    }
  }
  
  dim = 10000;
  if (VPixelRepn(src) == VUByteRepn) dim = 256;
  histo = (float *) VCalloc(dim,sizeof(float));
  for (j=0; j<dim; j++) histo[j] = 0;
  tiny = 2.0/(double)dim;

  a = ((double) dim) / (smax - smin);

  for (b=0; b<nbands; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	u = VGetPixel(src,b,r,c);
	if (ABS(u) < tiny) continue;
	j = (int) (a * (u - smin) + 0.5);
	if (j < 0) j = 0;
	if (j >= dim) j = dim-1;
	histo[j]++;
      }
    }
  }

 
  sum = 0;
  for (j=0; j<dim; j++) sum += histo[j];
  for (j=0; j<dim; j++) histo[j] /= sum;

  xmin = 0;
  sum  = 0;
  for (j=0; j<dim; j++) {
    sum += histo[j];
    if (sum > low) break;
  }
  xmin = ((double)j)/a + smin;
  
  xmax = dim;
  sum = 0;
  for (j=dim; j>0; j--) {
    sum += histo[j];
    if (sum > high) break;
  }
  xmax = ((double)j)/a + smin;

  slope = 255.0 / (xmax - xmin);

  dest_pp = (VUByte *) VImageData(dest);
  for (b=0; b<nbands; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	u = VGetPixel(src,b,r,c);
	v = (int) (slope * (u - xmin) + 0.5);
	if (ABS(u) < tiny) v = 0;
	if (v < 0) v = 0;
	if (v > 255) v = 255;
	*dest_pp++ = (VUByte) v;
      }
    }
  }

  VCopyImageAttrs (src, dest);
  return dest;
}
