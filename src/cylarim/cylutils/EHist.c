
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


void HistEqualize(Cylinders *cyl,VImage wimage,VImage metric,VImage rim)
{
  int b,r,c;
  size_t i,j,k,nbins = 10000;
  double u=0,hmin=-0.05,hmax=1.05,eps=0.0001;

  /* ini output image */
  VImage newmetric = VCreateImageLike(metric);
  VFillImage(newmetric,VAllBands,0);
  VFillImage(wimage,VAllBands,0);
  
  /* ini histogram */
  gsl_histogram *hist= gsl_histogram_alloc (nbins);
  gsl_histogram_set_ranges_uniform (hist,hmin,hmax);
  gsl_histogram_pdf *cdf  = gsl_histogram_pdf_alloc(nbins);

  for (i=0; i<cyl->numcylinders; i++) {

    hmin = 9999;
    hmax = 0;
    for (j=0; j<cyl->addr[i]->size; j++) {      
      k = cyl->addr[i]->data[j];
      b = gsl_matrix_int_get(cyl->xmap,k,0);
      r = gsl_matrix_int_get(cyl->xmap,k,1);
      c = gsl_matrix_int_get(cyl->xmap,k,2);
      u = VPixel(metric,b,r,c,VFloat);
      if (u < TINY) continue;
      if (u > hmax) hmax = u;
      if (u < hmin) hmin = u;
    }

    /* hist equalization */
    gsl_histogram_reset(hist);
    gsl_histogram_set_ranges_uniform (hist,hmin-eps,hmax+eps);

    for (j=0; j<cyl->addr[i]->size; j++) {

      k = cyl->addr[i]->data[j];
      b = gsl_matrix_int_get(cyl->xmap,k,0);
      r = gsl_matrix_int_get(cyl->xmap,k,1);
      c = gsl_matrix_int_get(cyl->xmap,k,2);
      u = VPixel(metric,b,r,c,VFloat);
      if (u < TINY) continue;
      gsl_histogram_increment(hist,u);
    }
    gsl_histogram_pdf_init(cdf,hist);
    
    /* fill output image */
    for (j=0; j<cyl->addr[i]->size; j++) {
      
      k = cyl->addr[i]->data[j];
      b = gsl_matrix_int_get(cyl->xmap,k,0);
      r = gsl_matrix_int_get(cyl->xmap,k,1);
      c = gsl_matrix_int_get(cyl->xmap,k,2);
      u = VPixel(metric,b,r,c,VFloat);
      if (u < TINY) continue;
      VPixel(wimage,b,r,c,VFloat) += 1.0;
      
      k=0;
      if (gsl_histogram_find(hist,u,&k) == GSL_SUCCESS) {
	VPixel(newmetric,b,r,c,VFloat) += (hmax*cdf->sum[k] + hmin);
      }
      else VWarning(" hist");
    }
  }

  VFloat *pd = VImageData(newmetric);
  VFloat *pw = VImageData(wimage);
  VFloat *px = VImageData(metric);
  for (i=0; i<VImageNPixels(wimage); i++) {
    if (pw[i] > 0.1) pd[i] /= pw[i];
    if (pd[i] < TINY && px[i] > 0 && px[i] < 0.1) pd[i] = hmin;
    if (pd[i] < TINY && px[i] > 0.9) pd[i] = hmax;
  }
  VCopyImagePixels(newmetric,metric,VAllBands);

  VDestroyImage(newmetric);
  gsl_histogram_pdf_free(cdf);
  gsl_histogram_free(hist);
}

