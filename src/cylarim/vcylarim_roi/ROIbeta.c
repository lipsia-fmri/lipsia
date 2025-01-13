/*
** cylarim: laminar statistics in a ROI
**
** G.Lohmann, MPI-KYB, Nov 2024
*/
#include <viaio/Vlib.h>
#include <viaio/VImage.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_histogram.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>



void ROIbeta(VImage *betaimage,VImage metric,VImage roi,VString filename)
{
  size_t i,j,nbins=9;
  double sum=0,hmin=0.01,hmax=1.01;

  
  gsl_histogram *hist[3];
  for (j=0; j<3; j++) {
    hist[j] = gsl_histogram_alloc (nbins);
    gsl_histogram_set_ranges_uniform (hist[j],hmin,hmax);
  }

  VBit *pr = VImageData(roi);
  VFloat *pb = VImageData(betaimage[0]);
  VFloat *pm = VImageData(metric);
  
  for (i=0; i<VImageNPixels(metric); i++) {
    if (pr[i] == 0) continue;
    for (j=0; j<3; j++) {
      pb = VImageData(betaimage[j]);
      if (pb[i] > 0)
	gsl_histogram_accumulate(hist[j],(double)pm[i],(double)pb[i]);
    }
  }
  for (j=0; j<3; j++) {
    sum = gsl_histogram_sum(hist[j]);
    gsl_histogram_scale(hist[j],1.0/sum);
  }

  
  FILE *fp = fopen(filename,"w");
  if (!fp) VError(" err opening %s",filename);

  double x,upper,lower,step;
  gsl_histogram_get_range (hist[0],0,&lower,&upper);
  step = 0.5*(upper-lower);
  
  for (i=0; i<nbins; i++) {
    gsl_histogram_get_range (hist[0],i,&lower,&upper);
    x = lower + step;
    fprintf(fp," %f ",x);
    for (j=0; j<3; j++) {
      fprintf(fp," %f",gsl_histogram_get(hist[j],i));
    }
    fprintf(fp,"\n");
  }
  fclose(fp);



}
