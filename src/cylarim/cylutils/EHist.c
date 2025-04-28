
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

extern void XWriteOutput(VImage image,VAttrList geolist,char *filename);

void HistEqualize(Cylinders *cyl,VImage wimage,VImage metric,VImage rim,VAttrList geolist)
{
  int b,r,c;
  size_t i,j,k,nbins = 10000;
  double u=0,hmin=-0.02,hmax=1.02,eps=0.001;
  
  double xmin=9999,xmax=-9999;
  VFloat *px = VImageData(metric);
  for (i=0; i<VImageNPixels(metric); i++) {
    if (px[i] > 0 && px[i] < xmin) xmin = px[i];
    if (px[i] > xmax) xmax = px[i];
  }
  if (xmax < xmin) VError(" HistEqualize, xmin: %f, xmax: %f",xmin,xmax);
  
  /* ini output image */
  VImage newmetric = VCreateImageLike(metric);
  VFillImage(newmetric,VAllBands,0);
  VFillImage(wimage,VAllBands,0);
  
  /* ini histogram */
  gsl_histogram *hist= gsl_histogram_alloc (nbins);
  gsl_histogram_set_ranges_uniform (hist,hmin,hmax);
  gsl_histogram_pdf *cdf  = gsl_histogram_pdf_alloc(nbins);

  for (i=0; i<cyl->numcylinders; i++) {
    if (cyl->addr[i]->size < 10) continue;
    hmin = 9999;
    hmax = -9999;
    for (j=0; j<cyl->addr[i]->size; j++) {      
      k = cyl->addr[i]->data[j];
      b = gsl_matrix_int_get(cyl->xmap,k,0);
      r = gsl_matrix_int_get(cyl->xmap,k,1);
      c = gsl_matrix_int_get(cyl->xmap,k,2);
      u = VPixel(metric,b,r,c,VFloat);
      if (u > hmax) hmax = u;
      if (u < hmin) hmin = u;
    }
    if (hmax < hmin) continue;

    /* hist equalization */
    gsl_histogram_reset(hist);
    gsl_histogram_set_ranges_uniform (hist,hmin-eps,hmax+eps);

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
      VPixel(wimage,b,r,c,VFloat) += 1.0;
      
      k=0;
      if (gsl_histogram_find(hist,u,&k) == GSL_SUCCESS) {
	VPixel(newmetric,b,r,c,VFloat) += (cdf->sum[k]-hmin)/(hmax-hmin);
      }
      else VWarning(" hist");
    }
  }

  /* rescale */
  VFloat qmin = 0.0001;
  VFloat qmax = 1.0;
  VFloat *pd = VImageData(newmetric);
  VFloat *pw = VImageData(wimage);
  for (i=0; i<VImageNPixels(wimage); i++) {
    if (pw[i] > 0.1) {
      pd[i] /= pw[i];
      pd[i] = (pd[i] - xmin)/(xmax-xmin);
      if (pd[i] < xmin) pd[i] = xmin;
      if (pd[i] > xmax) pd[i] = xmax;
    }
    /* add rim, if present */
    if (pw[i] < TINY && px[i] > 0.9) pd[i] = qmax;
    if (pd[i] < TINY && px[i] > 0 && px[i] < 0.1) pd[i] = qmin;
    if (pd[i] > qmax) pd[i] = qmax;
  }
  VCopyImagePixels(newmetric,metric,VAllBands);

  VDestroyImage(newmetric);
  gsl_histogram_pdf_free(cdf);
  gsl_histogram_free(hist);
}

