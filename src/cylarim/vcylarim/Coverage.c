/*
** Cylarim:  Cylinder-based Algorithm for Laminar fMRI Data Analysis
**
** G.Lohmann, MPI-KYB, Oct 2024
*/


/* From the standard C library: */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "../cylutils/cyl.h"


VImage VCylCover(VImage zmap,VImage rim,VImage wimage,Cylinders *cyl)
{
  size_t i,j,k,l;
  int b,r,c;
  double z=0,zx=0,zmax=0,zmin=0;
  VFloat *px=NULL,*pw=NULL;
  
  /* fill dest image */
  VFillImage(wimage,VAllBands,0);
  VImage dest = VCreateImageLike(zmap);
  VFillImage(dest,VAllBands,0);

  for (k=0; k<cyl->numcylinders; k++) {
    zmax=-99999;
    zmin=99999;
    z=zx=0;
    for (l=0; l<cyl->addr[k]->size; l++) {
      j = cyl->addr[k]->data[l];
      b = gsl_matrix_int_get(cyl->xmap,j,0);
      r = gsl_matrix_int_get(cyl->xmap,j,1);
      c = gsl_matrix_int_get(cyl->xmap,j,2);
      z = VPixel(zmap,b,r,c,VFloat);
      if (z > zmax) zmax = z;
      if (z < zmin) zmin = z;
    }
    zx = zmax;
    if (fabs(zmin) > zmax) zx = zmin;
    
    for (l=0; l<cyl->addr[k]->size; l++) {
      j = cyl->addr[k]->data[l];
      b = gsl_matrix_int_get(cyl->xmap,j,0);
      r = gsl_matrix_int_get(cyl->xmap,j,1);
      c = gsl_matrix_int_get(cyl->xmap,j,2);
      VPixel(dest,b,r,c,VFloat) += zx;
      VPixel(wimage,b,r,c,VFloat) += 1;
    }
  }

  pw = VImageData(wimage);
  px = VImageData(dest);
  for (i=0; i<VImageNPixels(dest); i++) {
    if (pw[i] > 0.001) { px[i] /= pw[i]; }
    else px[i] = 0;
  }  
  return dest;
}
