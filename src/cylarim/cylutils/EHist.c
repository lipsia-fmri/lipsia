
/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

  
/* From the standard C library: */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../cylutils/cyl.h"


void HistEqualize(Cylinders *cyl,VImage wimage,VImage metric)
{
  int b,r,c;
  size_t i,j,k,l,nbins = 10000;
  double u=0,w=0,hmin=-0.01,hmax=1.01;

  
  /* ini weight image */
  VFillImage(wimage,VAllBands,0);
  for (k=0; k<cyl->numcylinders; k++) {
    for (l=0; l<cyl->addr[k]->size; l++) {
      j = cyl->addr[k]->data[l];
      b = gsl_matrix_int_get(cyl->xmap,j,0);
      r = gsl_matrix_int_get(cyl->xmap,j,1);
      c = gsl_matrix_int_get(cyl->xmap,j,2);
      VPixel(wimage,b,r,c,VFloat) += 1.0;
    }
  }
  
  double wmax=-1;
  VFloat *pp = VImageData(wimage);
  for (i=0; i<VImageNPixels(wimage); i++) {
    if (pp[i] > wmax) wmax = pp[i];
  }

  
  /* ini output image */
  VImage dest = VCreateImageLike(metric);
  VFillImage(dest,VAllBands,0);

  
  /* ini histogram */
  gsl_histogram *hist= gsl_histogram_alloc (nbins);
  gsl_histogram_set_ranges_uniform (hist,hmin,hmax);
  gsl_histogram_pdf *cdf  = gsl_histogram_pdf_alloc(nbins);

  for (i=0; i<cyl->numcylinders; i++) {

    /* hist equalization */
    gsl_histogram_reset(hist);
    for (j=0; j<cyl->addr[i]->size; j++) {

      k = cyl->addr[i]->data[j];
      b = gsl_matrix_int_get(cyl->xmap,k,0);
      r = gsl_matrix_int_get(cyl->xmap,k,1);
      c = gsl_matrix_int_get(cyl->xmap,k,2);

      u = VPixel(metric,b,r,c,VFloat); 
      w = wmax - VPixel(wimage,b,r,c,VFloat);
      if (w < 0) w = 0;
      if (u <= hmin || u >= hmax) continue;
      gsl_histogram_accumulate(hist,u,w);
    }
    gsl_histogram_pdf_init(cdf,hist);

    
    /* fill output image */
    for (j=0; j<cyl->addr[i]->size; j++) {
      
      k = cyl->addr[i]->data[j];
      b = gsl_matrix_int_get(cyl->xmap,k,0);
      r = gsl_matrix_int_get(cyl->xmap,k,1);
      c = gsl_matrix_int_get(cyl->xmap,k,2);
      u = VPixel(metric,b,r,c,VFloat);
       
      k=0;
      if (gsl_histogram_find(hist,u,&k) == GSL_SUCCESS) {
	VPixel(dest,b,r,c,VFloat) += cdf->sum[k];
      }
      else VWarning(" hist");
    }
  }

  
  VFloat *pd = VImageData(dest);
  VFloat *pw = VImageData(wimage);
  for (i=0; i<VImageNPixels(wimage); i++) {
    if (pw[i] > 0.1) {
      pd[i] /= pw[i];
      if (pd[i] < 0.0001) pd[i] = 0.0001;
    }
  }
  VCopyImagePixels(dest,metric,VAllBands);
  
  VDestroyImage(dest);
  gsl_histogram_pdf_free(cdf);
  gsl_histogram_free(hist);
}

