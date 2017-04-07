/*! \file
  2D bicubic spline interpolation for image zooming.


\par Reference:
Michael J. Aramini.
"Efficient image magnification by bicubic spline interpolation".
http://members.bellatlantic.net/~vze2vrva/

\par Author:
Gabriele Lohmann, MPI-CBS
*/


/* From the Vista library: */
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <viaio/os.h>
#include <viaio/VImage.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


float a;


float C0(float t)
{
  return -a * t * t * t + a * t * t;
}

float C1(float t)
{
  return -(a + 2.0f) * t * t * t + (2.0f * a + 3.0f) * t * t - a * t;
}

float C2(float t)
{
  return (a + 2.0f) * t * t * t - (a + 3.0f) * t * t + 1.0f;
}

float C3(float t)
{
  return a * t * t * t - 2.0f * a * t * t + a * t;
}

/*!
\fn VImage VBicubicScale2d(VImage src,VImage dest,VFloat zoom)
\param src   input image (any repn)
\param dest  output image (any repn)
\param zoom  zoom factor
*/


VImage
VBicubicScale2d(VImage src,VImage dest,VFloat zoom)
{
  int  band=0,nbands;
  int  in_width, in_height, out_width, out_height, larger_out_dimension;
  int  n, d, j, k, l, m, index;
  float x;
  static int *L=NULL;
  static float **f=NULL, **g=NULL;
  static float *c[4];
  static float *h=NULL;
  VDouble xmin,xmax;


  nbands    = VImageNBands(src);
  in_height = VImageNRows(src);
  in_width  = VImageNColumns(src);

  out_height = (int)((float)in_height * zoom + 0.499);
  out_width  = (int)((float)in_width * zoom + 0.499);

  larger_out_dimension = (out_width > out_height) ? out_width : out_height;


  /*
  **  get n,d  so that zoom = n / d;
  */
  d = 10000;
  n = (int)(10000.0f * zoom);
  while (d%2 == 0 && n%2 == 0) {
    d /= 2;
    n /= 2;
  }
  while (d%3 == 0 && n%3 == 0) {
    d /= 3;
    n /= 3;
  }
  while (d%5 == 0 && n%5 == 0) {
    d /= 5;
    n /= 5;
  }
  while (d%7 == 0 && n%7 == 0) {
    d /= 7;
    n /= 7;
  }
  while (d%11 == 0 && n%11 == 0) {
    d /= 11;
    n /= 11;
  }
  while (d%13 == 0 && n%13 == 0) {
    d /= 13;
    n /= 13;
  }


  if (L == NULL) {
    L = (int *) VMalloc(sizeof(int) * larger_out_dimension);

    k = n;
    if (larger_out_dimension > k) k = larger_out_dimension;
    for (j=0; j<4; j++)
      c[j] = (float *) VMalloc(sizeof(int) * k);
  
    h = (float *) VMalloc(sizeof(int) * in_width);

    f = (float **) VMalloc(sizeof(float *) * in_height);
    for (j=0; j<in_height; j++)
      f[j] = (float *) VMalloc(sizeof(int) * in_width);
  
    g = (float **) VMalloc(sizeof(float *) * out_height);
    for (j=0; j<out_height; j++)
      g[j] = (float *) VMalloc(sizeof(int) * out_width);
  }


  for (k = 0; k < larger_out_dimension; k++)
    L[k] = (k * d) / n;

  for (k = 0; k < n; k++) {
    x = (float)((k * d) % n) / (float)n;
    c[0][k] = C0(x);
    c[1][k] = C1(x);
    c[2][k] = C2(x);
    c[3][k] = C3(x);
  }
  for (k = n; k < larger_out_dimension; k++)
    for (l = 0; l < 4; l++)
      c[l][k] = c[l][k % n];


  dest = VSelectDestImage("VBicubicScale2d",dest,nbands,out_height,out_width,VUByteRepn);
  VFillImage(dest,VAllBands,0);
  VCopyImageAttrs (src,dest);

  xmax = VPixelMaxValue(dest);
  xmin = VPixelMinValue(dest);

  for (band=0; band<nbands; band++) {

    for (j=0; j<in_height; j++) {
      for (k=0; k<in_width; k++) {
	f[j][k] = VGetPixel(src,band,j,k);
      }
    }

    for (k = 0; k < out_height; k++) {
      for (j = 0; j < in_width; j++) {
	h[j] = 0.0f;
	for (l = 0; l < 4; l++) {
	  index = L[k] + l - 1;
	  if ((index >= 0) && (index < in_height))
	    h[j] += f[index][j] * c[3 - l][k];
	}
      }
      for (m = 0; m < out_width; m++) {
	x = 0.5f;
	for (l = 0; l < 4; l++) {
	  index = L[m] + l - 1;
	  if ((index >= 0) && (index < in_width))
	    x += h[index] * c[3 - l][m];
	}

	if (x < xmin) 
	  x = xmin;
	else if (x > xmax)
	  x = xmax;

	g[k][m] = x;
      }
    }


    /*
    ** write to dest image
    */
    for (j=0; j<out_height; j++) {
      for (k=0; k<out_width; k++) {
	VSetPixel(dest,band,j,k,(VDouble) g[j][k]);
      }
    }
  }

  return dest;
}
