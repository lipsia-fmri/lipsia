
/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/Vlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../cylutils/cyl.h"

extern void gsl_sort_float(float *,size_t,size_t);

/*
** equivolume layering via histogram equalization
*/
void HistEqualize(Cylinders *cyl,VImage wimage,VImage metric,VImage rim)
{
  int b,r,c;
  size_t i,j,k,n;
  float u=0,cdf=0;

  /* ini output image */
  VImage newmetric = VCreateImageLike(metric);
  VFillImage(newmetric,VAllBands,0);
  VFillImage(wimage,VAllBands,0);
  
  size_t npix=0;
  for (i=0; i<cyl->numcylinders; i++) {
    if (cyl->addr[i]->size > npix) npix = cyl->addr[i]->size;
  }
  npix++;
  float *sorted = VCalloc(npix,sizeof(float));
  float *pp = VCalloc(npix,sizeof(float));

  
  /* hist equalization via sorting */
  for (i=0; i<cyl->numcylinders; i++) {

    n = cyl->addr[i]->size;
    for (j=0; j<n; j++) {
      k = cyl->addr[i]->data[j];
      b = gsl_matrix_int_get(cyl->xmap,k,0);
      r = gsl_matrix_int_get(cyl->xmap,k,1);
      c = gsl_matrix_int_get(cyl->xmap,k,2);
      u = VPixel(metric,b,r,c,VFloat);
      pp[j] = sorted[j] = u;
    }
    gsl_sort_float(sorted,1,j);
    
    /* refill  */   
    for (j=0; j<n; j++) {
      float val = pp[j];        
      size_t left = 0, right = n;
      while (left < right) {
	size_t mid = left + (right - left) / 2;
	if (sorted[mid] <= val) {
	  left = mid + 1;
	}
	else {
	  right = mid;
	}
      }
      cdf = (float)left / (float)n;

      k = cyl->addr[i]->data[j];
      b = gsl_matrix_int_get(cyl->xmap,k,0);
      r = gsl_matrix_int_get(cyl->xmap,k,1);
      c = gsl_matrix_int_get(cyl->xmap,k,2);
      
      VPixel(newmetric,b,r,c,VFloat) += cdf;
      VPixel(wimage,b,r,c,VFloat) += 1.0;
    }
  }

  
  /* rescale */
  VFloat *pd = VImageData(newmetric);
  VFloat *px = VImageData(metric);
  VFloat *pw = VImageData(wimage);
  for (i=0; i<VImageNPixels(wimage); i++) {
    if (pw[i] > 0.1) px[i] = pd[i]/pw[i];
  }

  VFillImage(wimage,VAllBands,0);
  VDestroyImage(newmetric);
  VFree(sorted);
  VFree(pp);
}

