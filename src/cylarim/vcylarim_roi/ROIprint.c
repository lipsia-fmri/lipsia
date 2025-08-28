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



void ROIprint(VImage zmap,VImage metric,VImage roi,int nbins,VString filename)
{
  size_t i,j;
  double u,v,z,jx;

  double nx = (double)nbins;
  double *s1 = (double *)VCalloc(nbins,sizeof(double));
  double *s2 = (double *)VCalloc(nbins,sizeof(double));
  double *kx = (double *)VCalloc(nbins,sizeof(double));
  double *mean = (double *)VCalloc(nbins,sizeof(double));
  double *sd = (double *)VCalloc(nbins,sizeof(double));

  VFloat *pz = VImageData(zmap);
  VFloat *pm = VImageData(metric);
  VBit   *pr = VImageData(roi);
  for (i=0; i<VImageNPixels(metric); i++) {
    if (pr[i] == 0) continue;
    u = pm[i];
    z = pz[i];
    if (u < TINY || u > 1.0) continue;
    if (fabs(z) < TINY) continue;
    for (j=1; j<=nbins; j++) {
      jx = (double)j;
      if (u > (jx-1.0)/nx && u <= jx/nx) {
	s1[j-1] += z;	
	s2[j-1] += z*z;
	kx[j-1]++;
      }
    }
  }
 
  for (j=0; j<nbins; j++) {
    if (kx[j] < 5) continue;
    mean[j] = s1[j]/kx[j];
    sd[j] = sqrt((s2[j] - kx[j] * mean[j] * mean[j]) / (kx[j] - 1.0));
  }
  
  FILE *fp = fopen(filename,"w");
  if (!fp) VError(" err opening %s",filename);

  for (j=0; j<nbins; j++) {
    u = (double)j/(double)nbins;
    v = (double)(j+1)/(double)nbins;
    fprintf(fp," %8.4f  %8.4f  %8.4f   %.2f\n",0.5*(u+v),mean[j],sd[j],kx[j]);
  }

  fclose(fp);
}
