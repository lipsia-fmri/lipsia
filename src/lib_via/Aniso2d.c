/*! \file
  2D Anisotropic diffusion

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
diffusion2d(float dx,float dy, VShort type,VFloat kappa,VFloat alpha)
{
  float u,v;
  double a;

  if (type == 0) {
    u = (dx*dx + dy*dy)/(kappa*kappa);
    v = (float) exp((double)-u);
  }
  
  else {
    a = 1.0 + alpha;
    u = sqrt((double)(dx*dx + dy*dy));
    u = u/kappa;
    v = pow((double)u,a);
    v = 1.0 / (1.0 + v);
  }
	 
  return v;
}


/*!
\fn VImage VAniso2d(VImage src,VImage dest,VShort numiter,
           VShort type,VFloat kappa,VFloat alpha);
\param src  input image 
\param dest  output image
\param numiter number of iterations
\param type type of diffusion function (0 or 1)
\param kappa parameter for diffusion function
\param alpha parameter for diffusion function
*/
VImage 
VAniso2d(VImage src,VImage dest,VShort numiter,
	 VShort type,VFloat kappa,VFloat alpha)
{
  VImage tmp1=NULL,tmp2=NULL;
  int nbands,nrows,ncols;
  int b,r,c,iter;
  float delta;
  float dx,dy,d,u,v;
  float ux1,ux2,uy1,uy2;
  float r1,r2,c1,c2;
  VDouble xmax,xmin;


  nbands = VImageNBands(src);
  nrows  = VImageNRows(src);
  ncols  = VImageNColumns(src);


  tmp1 = VCreateImage(1,nrows,ncols,VFloatRepn);
  tmp2 = VCreateImage(1,nrows,ncols,VFloatRepn);

  xmax = VPixelMaxValue (tmp1);
  xmin = VPixelMinValue (tmp1);

  dest = VCopyImage(src,dest,VAllBands);

  delta = 1.0 / 6.0;
  for (b=0; b<nbands; b++) {
    
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	VPixel(tmp1,0,r,c,VFloat) = VGetPixel(src,b,r,c);
      }
    }


    for (iter=0; iter < numiter; iter++) {

      for (r=1; r<nrows-1; r++) {
	for (c=1; c<ncols-1; c++) {

	  u  = VPixel(tmp1,0,r,c,VFloat);

	  c1 = VPixel(tmp1,0,r,c+1,VFloat);
	  c2 = VPixel(tmp1,0,r,c-1,VFloat);

	  r1 = VPixel(tmp1,0,r+1,c,VFloat);
	  r2 = VPixel(tmp1,0,r-1,c,VFloat);

	  /* col-dir */
	  dx = c1-u;
	  dy = r1-r2;
	  d  = diffusion2d(dx,dy,type,kappa,alpha);
	  ux1 = d*(c1 - u);

	  dx = u-c2;
	  d  = diffusion2d(dx,dy,type,kappa,alpha);
	  ux2 = d*(u - c2);


	  /* row-dir */
	  dx = c1-c2;
	  dy = r1-u;
	  d  = diffusion2d(dx,dy,type,kappa,alpha);
	  uy1 = d*(r1 - u);

	  dy = u-r2;
	  d  = diffusion2d(dx,dy,type,kappa,alpha);
	  uy2 = d*(u - r2);


	  /* sum */
	  v = u + delta*(ux1 - ux2 + uy1 - uy2);

	  if (v > xmax) v = xmax;
	  if (v < xmin) v = xmin;
	  VPixel(tmp2,0,r,c,VFloat) = v;
	}
      }
      tmp1 = VCopyImagePixels(tmp2,tmp1,VAllBands);
    }
    
    /*
    ** output
    */
    for (r=1; r<nrows-1; r++) {
      for (c=1; c<ncols-1; c++) {
	v = VPixel(tmp2,0,r,c,VFloat);
	if (v > VPixelMaxValue (dest)) v = VPixelMaxValue (dest);
	if (v < VPixelMinValue (dest)) v = VPixelMinValue (dest);
	VSetPixel(dest,b,r,c,(VDouble) v);
      }
    }
  }

  VDestroyImage(tmp1);
  VDestroyImage(tmp2);

  return dest;
}
