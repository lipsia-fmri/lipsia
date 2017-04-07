/*! \file
 compute gradient magnitude

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

/*!
\fn VImage VMagnitude3d(VImage gradz,VImage grady,VImage gradx,VImage dest)
\param gradz  input gradient in slice direction (float repn)
\param grady  input gradient in row direction (float repn)
\param gradx  input gradient in column direction (float repn)
\param dest   output image (float repn)
*/
VImage
VMagnitude3d(VImage gradz,VImage grady,VImage gradx,VImage dest)
{
  int nbands,nrows,ncols;
  int i,npixels;
  VFloat *gradx_pp,*grady_pp,*gradz_pp,*dest_pp;
  float x,y,z,u;

  nbands  = VImageNBands (gradz);
  nrows   = VImageNRows (gradz);
  ncols   = VImageNColumns (gradz);
  npixels = VImageNPixels (gradz);

  gradx_pp = (VFloat *) VImageData(gradx); 
  grady_pp = (VFloat *) VImageData(grady); 
  gradz_pp = (VFloat *) VImageData(gradz); 

  dest = VSelectDestImage("VMagnitude3d",dest,nbands,nrows,ncols,VFloatRepn);
  VCopyImageAttrs (gradx,dest);

  dest_pp = (VFloat *) VImageData(dest);

  for (i=0; i<npixels; i++) {
    x = (float) (*gradx_pp++); 
    y = (float) (*grady_pp++);
    z = (float) (*gradz_pp++);
    u = sqrt((double)(x*x + y*y + z*z));
    *dest_pp++ = u;
  }

  return dest;
}


/*!
\fn VImage VMagnitude2d(VImage grady,VImage gradx,VImage dest)
\param grady  input gradient in row direction (float repn)
\param gradx  input gradient in column direction (float repn)
\param dest   output image (float repn)
*/
VImage
VMagnitude2d(VImage grady,VImage gradx,VImage dest)
{
  int nbands,nrows,ncols;
  int i,npixels;
  VFloat *gradx_pp,*grady_pp,*dest_pp;
  float x,y,u;

  nbands  = VImageNBands (gradx);
  nrows   = VImageNRows (gradx);
  ncols   = VImageNColumns (gradx);
  npixels = VImageNPixels (gradx);

  dest = VSelectDestImage("VMagnitude2d",dest,nbands,nrows,ncols,VFloatRepn);
  VCopyImageAttrs (gradx,dest);


  gradx_pp = (VFloat *) VImageData(gradx); 
  grady_pp = (VFloat *) VImageData(grady); 
  dest_pp  = (VFloat *) VImageData(dest); 

  for (i=0; i<npixels; i++) {
    x = (float) (*gradx_pp++); 
    y = (float) (*grady_pp++);
    u = sqrt((double)(x*x + y*y));
    *dest_pp++ = u;
  }

  return dest;
}
