
/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/Vlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../cylutils/cyl.h"

/*
** equivolume layering via histogram equalization
*/
void HistEqualize(Cylinders *cyl,VImage wimage,VImage metric,VImage rim,VAttrList geolist)
{
  int b,r,c;
  size_t i,j,k,nbins = 10000;
  double u=0;
  double xmin = 0.0001;
  double hmin = -0.0001;
  double hmax = 1.02;

  /* ini output image */
  VImage newmetric = VCreateImageLike(metric);
  VFillImage(newmetric,VAllBands,0);
  VFillImage(wimage,VAllBands,0);
  
  /* ini histogram */
  gsl_histogram *hist= gsl_histogram_alloc (nbins);
  gsl_histogram_set_ranges_uniform (hist,hmin,hmax);
  gsl_histogram_pdf *cdf  = gsl_histogram_pdf_alloc(nbins);

  double lower,upper;
  gsl_histogram_get_range(hist,0,&lower,&upper);
  double add = (upper-lower)*0.5;  /* half bin size */
  
  for (i=0; i<cyl->numcylinders; i++) {
    if (cyl->addr[i]->size < 10) continue;
    
    gsl_histogram_reset(hist);
    for (j=0; j<cyl->addr[i]->size; j++) {

      k = cyl->addr[i]->data[j];
      b = gsl_matrix_int_get(cyl->xmap,k,0);
      r = gsl_matrix_int_get(cyl->xmap,k,1);
      c = gsl_matrix_int_get(cyl->xmap,k,2);
      u = VPixel(metric,b,r,c,VFloat);
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
      
      k=0;
      if (gsl_histogram_find(hist,u,&k) == GSL_SUCCESS) {
	VPixel(newmetric,b,r,c,VFloat) += (cdf->sum[k] + add);
	VPixel(wimage,b,r,c,VFloat) += 1.0;
      }
      else VWarning(" hist equalization: %f",u);
    }
  }

  /* rescale */
  VFloat *pd = VImageData(newmetric);
  VFloat *pw = VImageData(wimage);
  VFloat *px = VImageData(metric);
  for (i=0; i<VImageNPixels(wimage); i++) {
    if (pw[i] > 0.1) {
      pd[i] /= pw[i];
      if (pd[i] < xmin) pd[i] = xmin;
    }

    /* add rim */
    else if (px[i] > TINY) {
      if (px[i] < 0.1) pd[i] = xmin;
      if (px[i] > 0.9) pd[i] = 1.0;
    }
  }
  VCopyImagePixels(newmetric,metric,VAllBands);

  VDestroyImage(newmetric);
  gsl_histogram_pdf_free(cdf);
  gsl_histogram_free(hist);
}

