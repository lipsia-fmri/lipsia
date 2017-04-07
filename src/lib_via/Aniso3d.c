/*! \file
  3D Anisotropic diffusion

\par Reference:
  G.Gerig, O. Kuebler, R. Kikinis, F.Jolesz:
  Nonlinear anisotropic filtering of MRI data.
  IEEE Trans on Medical Imaging, Vol.11, No.2, June 1992.

\par Author:
  G.Lohmann, MPI-CBS
*/



#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define SQR(x) ((x) * (x))
#define ABS(x) ((x) > 0 ? (x) : -(x))


float
diffusion3d(float dx,float dy,float dz,VShort type,VFloat kappa,VFloat alpha)
{
  float u,v;
  double a;

  if (type == 0) {
    u = (dx*dx + dy*dy + dz*dz)/(kappa*kappa);
    v = (float) exp((double)-u);
  }
  
  else {
    a = 1.0 + alpha;
    u = sqrt((double)(dx*dx + dy*dy + dz*dz));
    u = u/kappa;
    v = (float) pow((double)u,a);
    v = 1.0 / (1.0 + v);
  }
	 
  return v;
}


/*!
\fn VImage VAniso3d(VImage src,VImage dest,VShort numiter,
           VShort type,VFloat kappa,VFloat alpha);
\param src     input image
\param dest    output image
\param numiter number of iterations
\param type    type of diffusion function (0 or 1)
\param kappa   parameter for diffusion function
\param alpha   parameter for diffusion function
*/

VImage 
VAniso3d(VImage src,VImage dest,VShort numiter,
	 VShort type,VFloat kappa,VFloat alpha)
{
  VImage tmp1=NULL,tmp2=NULL;
  int nbands,nrows,ncols;
  int b,r,c,iter;
  float delta;
  float dx,dy,dz,d,u,v;
  float ux1,ux2,uy1,uy2,uz1,uz2;
  float b1,b2,r1,r2,c1,c2;
  VDouble xmax,xmin;
  VBoolean ignore = TRUE;


  nbands = VImageNBands(src);
  nrows  = VImageNRows(src);
  ncols  = VImageNColumns(src);

  if (nbands < 3) VError(" min number of slices is 3");

  tmp1 = VConvertImageCopy(src,NULL,VAllBands,VFloatRepn);
  tmp2 = VCreateImage(nbands,nrows,ncols,VFloatRepn);
  VFillImage(tmp2,VAllBands,0);

  xmax = VPixelMaxValue (tmp1);
  xmin = VPixelMinValue (tmp1);

  delta = 1.0 / 7.0;

  dx = dy = dz = 0;

  for (iter=0; iter < numiter; iter++) {

    for (b=1; b<nbands-1; b++) {
      for (r=1; r<nrows-1; r++) {
	for (c=1; c<ncols-1; c++) {

	  u  = VPixel(tmp1,b,r,c,VFloat);
	  if (ignore && ABS(u) < 1.0e-10) continue;

	  c1 = VPixel(tmp1,b,r,c+1,VFloat);
	  c2 = VPixel(tmp1,b,r,c-1,VFloat);

	  r1 = VPixel(tmp1,b,r+1,c,VFloat);
	  r2 = VPixel(tmp1,b,r-1,c,VFloat);

	  b1 = VPixel(tmp1,b+1,r,c,VFloat);
	  b2 = VPixel(tmp1,b-1,r,c,VFloat);

	  /* col-dir */
	  dx = c1-u;
	  dy = r1-r2;
	  dz = b1-b2;
	  d  = diffusion3d(dx,dy,dz,type,kappa,alpha);
	  ux1 = d*(c1 - u);

	  dx = u-c2;
	  d  = diffusion3d(dx,dy,dz,type,kappa,alpha);
	  ux2 = d*(u - c2);


	  /* row-dir */
	  dx = c1-c2;
	  dy = r1-u;
	  dz = b1-b2;
	  d  = diffusion3d(dx,dy,dz,type,kappa,alpha);
	  uy1 = d*(r1 - u);

	  dy = u-r2;
	  d  = diffusion3d(dx,dy,dz,type,kappa,alpha);
	  uy2 = d*(u - r2);


	  /* slice-dir */
	  dx = c1-c2;
	  dy = r1-r2;
	  dz = b1-u;
	  d  = diffusion3d(dx,dy,dz,type,kappa,alpha);
	  uz1 = d*(b1 - u);

	  dz = u-b2;
	  d  = diffusion3d(dx,dy,dz,type,kappa,alpha);
	  uz2 = d*(u - b2);

	  /* sum */
	  v = u + delta*(ux1 - ux2 + uy1 - uy2 + uz1 - uz2);

	  if (v > xmax) v = xmax;
	  if (v < xmin) v = xmin;
	  VPixel(tmp2,b,r,c,VFloat) = v;
	}
      }
    }
    tmp1 = VCopyImagePixels(tmp2,tmp1,VAllBands);
  }


  /*
  ** output
  */
  dest = VCopyImage(src,dest,VAllBands);

  xmax = VPixelMaxValue (dest);
  xmin = VPixelMinValue (dest);

  for (b=1; b<nbands-1; b++) {
    for (r=1; r<nrows-1; r++) {
      for (c=1; c<ncols-1; c++) {
	v = VPixel(tmp2,b,r,c,VFloat);
	if (v > xmax) v = xmax;
	if (v < xmin) v = xmin;
	VSetPixel(dest,b,r,c,(VDouble) v);
      }
    }
  }

  VDestroyImage(tmp1);
  VDestroyImage(tmp2);

  return dest;
}
