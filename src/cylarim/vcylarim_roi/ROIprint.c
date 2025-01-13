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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


extern double gaussian(double x,int);


void ROIprint(VImage zmap,VImage metric,VImage *betaimage,int nbeta,VImage roi)
{
  size_t i,j,k;


  VBit *pr = VImageData(roi);
  size_t npix=0;
  for (i=0; i<VImageNPixels(roi); i++) {
    if (pr[i] != 0) npix++;
  }

    
  /* print zmap vs metric */
  VFloat *pz = VImageData(zmap);
  VFloat *pm = VImageData(metric);
  FILE *fp = fopen("pts.txt","w");  
  for (i=0; i<VImageNPixels(zmap); i++) {
    if (pr[i] == 0) continue;
    fprintf(fp," %f %f\n",pm[i],pz[i]);
  }
  fclose(fp);
  if (nbeta < 4) return;  /* no model */

  
  /* read beta data */
  gsl_matrix *D = gsl_matrix_calloc(nbeta,npix);
  VFloat *pb = NULL;
  double *pp;

  for (j=0; j<nbeta; j++) {
    pb = VImageData(betaimage[j]);
    pp = gsl_matrix_ptr(D,j,0);
    k=0;
    for (i=0; i<VImageNPixels(roi); i++) {
      if (pr[i] == 0) continue;
      pp[k] = (double)(pb[i]);
      k++;
    }
  }


  /* mean beta values */
  gsl_vector *beta = gsl_vector_calloc(nbeta);
  for (j=0; j<nbeta; j++) {
    pp = gsl_matrix_ptr(D,j,0);
    beta->data[j] = gsl_stats_mean(pp,1,npix);
  }


  /* print model */
  FILE *fq = fopen("model.txt","w");
  double x=0,z=0;
  
  for (x=0.0; x<1.01; x+=0.05) {
    z=0;
    for (j=0; j<nbeta-1; j++) {
      z += beta->data[j]*gaussian(x,(int)j);
    }
    if (nbeta > 3) z += beta->data[3];
    fprintf(fq," %f %f\n",x,z);
  }
  fclose(fq);

  gsl_matrix_free(D);
}
