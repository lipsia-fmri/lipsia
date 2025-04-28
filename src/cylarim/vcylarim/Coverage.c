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
  double z=0,zmax=0;
  VFloat *px=NULL,*pw=NULL;
  
  /* fill dest image */
  VFillImage(wimage,VAllBands,0);
  VImage dest = VCreateImageLike(zmap);
  VFillImage(dest,VAllBands,0);

  for (k=0; k<cyl->numcylinders; k++) {
    zmax=0;
    for (l=0; l<cyl->addr[k]->size; l++) {
      j = cyl->addr[k]->data[l];
      b = gsl_matrix_int_get(cyl->xmap,j,0);
      r = gsl_matrix_int_get(cyl->xmap,j,1);
      c = gsl_matrix_int_get(cyl->xmap,j,2);
      if (VPixel(rim,b,r,c,VUByte) != 3) continue;
      z = VPixel(zmap,b,r,c,VFloat);
      if (z > zmax) zmax = z;
    }
    for (l=0; l<cyl->addr[k]->size; l++) {
      j = cyl->addr[k]->data[l];
      b = gsl_matrix_int_get(cyl->xmap,j,0);
      r = gsl_matrix_int_get(cyl->xmap,j,1);
      c = gsl_matrix_int_get(cyl->xmap,j,2);
      if (VPixel(rim,b,r,c,VUByte) != 3) continue;
      VPixel(dest,b,r,c,VFloat) += zmax;
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


void ZmapNormalize(VImage zmap,VImage rim)
{
  size_t i;
  VFloat *pp = VImageData(zmap);
  VUByte *pu = VImageData(rim);

  double s1=0,s2=0,nx=0;
  for (i=0; i<VImageNPixels(zmap); i++) {
    if (pu[i] > 0) {
      s1 += pp[i];
      s2 += pp[i]*pp[i];
      nx++;
    }
  }
  if (nx < 0.1) VError(" zmap is zero");
  double mean = s1/nx;
  double var = (s2 - nx * mean * mean) / (nx - 1.0);
  if (var < TINY) VError(" no variance in zmap");
  double sd = sqrt(var);

  for (i=0; i<VImageNPixels(zmap); i++) {
    if (pu[i] > 0) {
      pp[i] = (pp[i] - mean)/sd;
    }
    else {
      pp[i] = 0;
    }
  }
}
