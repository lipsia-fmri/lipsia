/*
** False disovery rate (FDR)
**
** G.Lohmann, Jan 2017
*/

#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_math.h>


/* get array of FDR values */
double *GetFdr(gsl_histogram *nullhist,gsl_histogram *realhist)
{
  size_t i=0;
  size_t nbins = gsl_histogram_bins (nullhist);

  /* cumulative distribution functions (CDF) */
  gsl_histogram_pdf *cdfz = gsl_histogram_pdf_alloc(nbins);
  gsl_histogram_pdf *cdf0 = gsl_histogram_pdf_alloc(nbins);
  gsl_histogram_pdf_init (cdfz,realhist);
  gsl_histogram_pdf_init (cdf0,nullhist);

  /* fdr array */
  double *fdr = (double *)VCalloc(nbins,sizeof(double));
  for (i=0; i<nbins; i++) {
    double Fz = cdfz->sum[i];
    double F0 = cdf0->sum[i];
    double xFdr = 1.0;
    if (Fz  < 1.0) xFdr = (1.0-F0)/(1.0-Fz);
    if (xFdr > 1.0) xFdr = 1.0;
    if (xFdr < 0.0) xFdr = 0.0;
    fdr[i] = xFdr;
  }

  /* enforce monotonicity */
  fdr[0] = 1.0;
  for (i=1; i<nbins; i++) {
    if (fdr[i] > fdr[i-1]) fdr[i] = fdr[i-1];
  }
  return fdr;
}


/* thresholding */
void ApplyFdrThreshold(VImage dest,double alpha)
{
  double u=0,beta=1.0-alpha;
  size_t i,n=0;
  VFloat *pp = VImageData(dest); 
  for (i=0; i<VImageNPixels(dest); i++) {
    u = (*pp);
    if (u < beta) *pp = 0;
    else n++;
    pp++;
  }
  fprintf(stderr," Number of voxels at Fdr < %.4f:  %lu\n",alpha,n);
}


/* estimage false discovery rates */
void FDR(VImage src,VImage dest,gsl_histogram *nullhist,gsl_histogram *realhist,double alpha)
{
  size_t j;
  double u=0,tiny=1.0e-8;
  double hmin = gsl_histogram_min (realhist);
  double hmax = gsl_histogram_max (realhist);


  /* compute Fdr */
  double *Fdr = GetFdr(nullhist,realhist);


  /* apply FDR */
  VFillImage(dest,VAllBands,0);
  int b,r,c;
  for (b=0; b<VImageNBands(dest); b++) {
    for (r=0; r<VImageNRows(dest); r++) {
      for (c=0; c<VImageNColumns(dest); c++) {

	u = VGetPixel(src,b,r,c);
	if (u > hmax) u = hmax-tiny;
	if (u < hmin) u = hmax+tiny;
	if (fabs(u) > 0) {
	  if (gsl_histogram_find (realhist,(double)u,&j) == GSL_SUCCESS) {
	    VPixel(dest,b,r,c,VFloat) = 1.0 - Fdr[j];
	  }
	}
      }
    }
  }

  /* apply threshold if needed */
  if (alpha < 0.99) ApplyFdrThreshold(dest,alpha);
  return;
}
