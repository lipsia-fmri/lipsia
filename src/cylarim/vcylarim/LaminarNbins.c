
/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/Vlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>


/* From the standard C library: */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "../cylutils/cyl.h"


double LaminarNBins(gsl_vector *y,gsl_vector *mvec,gsl_vector *beta,int verbose)
{
  size_t i,n=y->size;
  int j,nbins=beta->size;
  double u=0,w=0,jx=0;
  double nx = (double)nbins;

  gsl_vector_set_zero(beta);
 
  double *sum = (double *)VCalloc(nbins,sizeof(double));
  double *kx = (double *)VCalloc(nbins,sizeof(double));
  
  for (i=0; i<n; i++) {
    u = y->data[i];
    w = mvec->data[i];
    for (j=1; j<=nbins; j++) {
      jx = (double)j;
      if (w > (jx-1.0)/nx && w <= jx/nx) {
	sum[j-1] += u;
	kx[j-1] += 1.0;
      }
    }
  }

  for (j=0; j<nbins; j++) {
    beta->data[j] = 0;
    if (kx[j] > 0) beta->data[j] = sum[j]/kx[j];
  }

  /*
  if (verbose) {
    FILE *fp = fopen("l.txt","w");
    for (i=0; i<n; i++) {

      w = mvec->data[i];
      for (j=1; j<=nbins; j++) {
	jx = (double)j;
	if (w > (jx-1.0)/nx && w <= jx/nx) {
	  fprintf(fp," %f  %d\n",w,j-1);
	}
      }
    }
    fclose(fp);
    exit(0);
  }
  */

  VFree(sum);
  VFree(kx);
  return 1;
}

