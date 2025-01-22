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
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>



void ROIplot(VImage *betaimage,VImage roi,VString filename)
{
  size_t i,j,n;
  double q0,q1,q2;;

  VBit *pr = VImageData(roi);
  VFloat *pb = VImageData(betaimage[0]);


  n=0;
  for (i=0; i<VImageNPixels(roi); i++) {
    if (pr[i] == 0) continue;
    n++;
  }

  gsl_matrix *B = gsl_matrix_calloc(3,n);
  
  for (j=0; j<3; j++) {
    pb = VImageData(betaimage[j]);
    n=0;
    for (i=0; i<VImageNPixels(roi); i++) {
      if (pr[i] == 0) continue;
      gsl_matrix_set(B,j,n,(double)pb[i]);    
      n++;
    }
  }

  FILE *fp = fopen(filename,"w");
  if (!fp) VError(" err opening %s",filename);

  double *pp;
  for (j=0; j<3; j++) {
    pp = gsl_matrix_ptr(B,j,0);
    gsl_sort(pp,1,n);
    q0 = gsl_stats_quantile_from_sorted_data(pp,1,n,0.25);
    q1 = gsl_stats_quantile_from_sorted_data(pp,1,n,0.50);
    q2 = gsl_stats_quantile_from_sorted_data(pp,1,n,0.75);
    fprintf(stderr," %3lu   %8.5f  %8.5f  %8.5f\n",j,q1,q0,q2);
    fprintf(fp," %3lu   %8.5f  %8.5f  %8.5f\n",j,q1,q0,q2);
  }

  fclose(fp);
}
