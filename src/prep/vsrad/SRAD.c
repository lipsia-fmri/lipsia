/*
**
** Lit: 
**   Y.Yu, S.T. Acton:
**   Speckle reducing anisotropic diffusion.
**   IEEE Trans Image Proc., Vol. 11, No. 11, Nov 2002
**
**  G.Lohmann, MPI-CBS, 2008
*/


#include <viaio/Vlib.h>
#include <viaio/VImage.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>


#define SQR(x) ((x) * (x))
#define ABS(x) ((x) > 0 ? (x) : -(x))


/*
** repair border
*/
void
Border(VImage src)
{
  int b=0,r,c,nrows,ncols;
  double u;

  nrows = VImageNRows(src);
  ncols = VImageNColumns(src);

  for (r=0; r<VImageNRows(src); r++) {
    u = VGetPixel(src,b,r,1);
    VSetPixel(src,b,r,0,u);
    u = VGetPixel(src,b,r,ncols-2);
    VSetPixel(src,b,r,ncols-1,u);
  }

  for (c=0; c<VImageNColumns(src); c++) {
    u = VGetPixel(src,b,1,c);
    VSetPixel(src,b,0,c,u);
    u = VGetPixel(src,b,nrows-2,c);
    VSetPixel(src,b,nrows-1,c,u);
  }
}


VImage 
VSRAD(VImage src,VImage dest,VShort numiter,VShort type,VFloat rho)
{
  static VImage tmp1=NULL,tmp2=NULL,ctmp=NULL;
  int nbands,nrows,ncols;
  int b,r,c,iter;
  double u,v,w,x,y;
  double d,da,dd,dt,t,q,q0;
  double cx,c0,c1,c2,c3;
  double tiny=1.0e-6;


  nbands = VImageNBands(src);
  nrows  = VImageNRows(src);
  ncols  = VImageNColumns(src);


  if (tmp1 == NULL) {
    tmp1 = VCreateImage(1,nrows,ncols,VFloatRepn);
    tmp2 = VCreateImage(1,nrows,ncols,VFloatRepn);
    ctmp = VCreateImage(1,nrows,ncols,VFloatRepn);
  }

  dt = 0.05;

  dest = VCopyImage(src,dest,VAllBands);

  for (b=0; b<nbands; b++) {
    
    t = 0;
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	VPixel(tmp1,0,r,c,VFloat) = VGetPixel(src,b,r,c);
      }
    }

    VFillImage(ctmp,VAllBands,0);
    VFillImage(tmp2,VAllBands,0);


    q0 = 1;
    for (iter=0; iter < numiter; iter++) {

      for (r=1; r<nrows-1; r++) {
	for (c=1; c<ncols-1; c++) {


	  u  = VPixel(tmp1,0,r,c,VFloat);
	  if (ABS(u) < tiny) continue;

	  v  = VPixel(tmp1,0,r+1,c,VFloat);
	  w  = VPixel(tmp1,0,r,c+1,VFloat);
	  x  = VPixel(tmp1,0,r-1,c,VFloat);
	  y  = VPixel(tmp1,0,r,c-1,VFloat);

	  da = (SQR(u-v) + SQR(u-w) + SQR(u-x) + SQR(u-y))/(u*u);
	  dd = (v+w+x+y - 4.0*u)/(u*u);
	  q = (0.5*da - 0.0625*dd) / (SQR(1.0 + 0.25*dd));
	  if (q > 0) {
	    q = sqrt(q);
	    u = (q*q - q0*q0) / ((q0*q0)*(1.0 + q0*q0));
	  }
	  else
	    u = 0;

	  cx = 0;
	  if (type == 0) 
	    cx = 1.0 / (1.0 + u);
	  else 
	    cx = exp(-u);

	  if (gsl_isnan(cx) || gsl_isinf(cx)) {
	    cx = 0;
	    VWarning(" cx, insnan, isinf");
	  }

	  VPixel(ctmp,0,r,c,VFloat) = cx;
	}
      }
      Border(ctmp);
      
      for (r=1; r<nrows-1; r++) {
	for (c=1; c<ncols-1; c++) {
	  
	  u  = VPixel(tmp1,0,r,c,VFloat);
	  if (ABS(u) < tiny) continue;

	  v  = VPixel(tmp1,0,r+1,c,VFloat);
	  w  = VPixel(tmp1,0,r,c+1,VFloat);
	  x  = VPixel(tmp1,0,r-1,c,VFloat);
	  y  = VPixel(tmp1,0,r,c-1,VFloat);

	  c0  = VPixel(ctmp,0,r+1,c,VFloat);
	  c1  = VPixel(ctmp,0,r,c+1,VFloat);
	  c2  = VPixel(ctmp,0,r-1,c,VFloat);
	  c3  = VPixel(ctmp,0,r,c-1,VFloat);

          /*  ???
          cx =  VPixel(ctmp,0,r,c,VFloat);
	  d = c0*(v-u) + c1*(w-u) + cx*(x-u) + cx*(y-u); 
          */
	  
          d = c0*(v-u) + c1*(w-u) + c2*(x-u) + c3*(y-u);


	  u = u + 0.25*dt * d;
	  if (gsl_isnan(u) || gsl_isinf(u)) {
	    u = 0;
	    VWarning(" insnan, isinf");
	  }

	  if (u > VPixelMaxValue (tmp2)) u = VPixelMaxValue (tmp2);
	  if (u < VPixelMinValue (tmp2)) u = VPixelMinValue (tmp2);
	  VPixel(tmp2,0,r,c,VFloat) = u;

	}
      }
      q0 *= exp(-rho*t);

      Border(tmp2);
      tmp1 = VCopyImagePixels(tmp2,tmp1,VAllBands);
      t += dt;
    }
    
    /*
    ** output
    */
    for (r=1; r<nrows-1; r++) {
      for (c=1; c<ncols-1; c++) {
	v = VPixel(tmp2,0,r,c,VFloat);
	if (gsl_isnan(v) || gsl_isinf(v)) {
	  v = 0;
	  VWarning(" insnan, isinf");
	}

	if (v > VPixelMaxValue (dest)) v = VPixelMaxValue (dest);
	if (v < VPixelMinValue (dest)) v = VPixelMinValue (dest);
	VSetPixel(dest,b,r,c,(VDouble) v);
      }
    }
  }

  /*
  VDestroyImage(tmp1);
  VDestroyImage(tmp2);
  */
  return dest;
}
