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

extern double t2z(double,double);
extern void GetResolution(VImage src,gsl_vector *reso);


/* simplified formula for effective degrees of freedom */ 
double DegreesOfFreedom(size_t npix,gsl_vector *reso_orig,gsl_vector *reso)
{
  double nx = (double)npix;
  double vol_orig = reso_orig->data[0]*reso_orig->data[1]*reso_orig->data[2];
  double vol_resampled = reso->data[0]*reso->data[1]*reso->data[2];
  double npix_orig = nx*vol_resampled/vol_orig;
  double dof_orig = npix_orig - 2.0;  /* two model parameters: mean,sd */
  double fwhm = 2.5*vol_orig;  /* default assumption: 2.5 times the voxel size */
  double dof = dof_orig * vol_orig/fwhm;
  if (dof < 1) dof = 1;   /* prevent divide by zero */

  fprintf(stderr,"\n num voxels: %lu,  original number of voxels: %.2f\n",npix,npix_orig);
  fprintf(stderr," estimated degrees of freedom: %.4f\n\n",dof);

  return dof;
}


/* paired t-test */
void ROIstats(VImage *betaimage,VImage *covimage,VImage edfimage,int nbeta,int ncov,VImage roi)
{
  size_t i,j,k;

  VBit *pr = VImageData(roi);
  size_t npix=0;
  for (i=0; i<VImageNPixels(roi); i++) {
    if (pr[i] != 0) npix++;
  }
 
    
  /* read data */
  gsl_matrix *D = gsl_matrix_calloc(nbeta,npix);
  gsl_vector *edf = gsl_vector_calloc(npix);
  VFloat *pb,*pe;
  double *pp;

  k=0;
  pe = VImageData(edfimage);
  for (i=0; i<VImageNPixels(roi); i++) {
    if (pr[i] == 0) continue;
    edf->data[k++] = (double) pe[i];
  }

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


  /* mean effective degrees of freedom */
  double nx=(double)npix;
  double s0=0,meanedf=0;
  for (i=0; i<npix; i++) {
    s0 += edf->data[i];
  }
  meanedf = s0/nx;
  fprintf(stderr,"mean EDF: %.4f\n\n",meanedf);

  /* onesample t-tests */
  double mean=0,sd,se=0,t=0;
  
  fprintf(stderr," Onesample tests:\n");
  fprintf(stderr,"         mean       stddev     stderr       T\n");
  fprintf(stderr," ------------------------------------------------\n");
  for (j=0; j<3; j++) {
    pp = gsl_matrix_ptr(D,j,0);
    mean = gsl_stats_mean(pp,1,npix);
    sd = gsl_stats_sd_m(pp,1,npix,mean);
    se = sd/meanedf;
    t = mean/se;
    fprintf(stderr," %3lu:  %8.4f   %8.4f   %8.4f   %9.3f\n",j,mean,sd,se,t);
  }
  fprintf(stderr,"\n\n");
  

  /* paired t-tests */
  fprintf(stderr," Paired tests:\n");
  fprintf(stderr,"                        mean       stddev     stderr       T\n");
  fprintf(stderr," ---------------------------------------------------------------\n");
  double *p0 = gsl_matrix_ptr(D,0,0);
  double *p1 = gsl_matrix_ptr(D,1,0);
  double *p2 = gsl_matrix_ptr(D,2,0);

  double s1=0,s2=0,d=0;
  s0=s1=s2=0;
  for (i=0; i<npix; i++) {
    d = p1[i]-p0[i];
    s1 += d;
    s2 += d*d;
  }
  mean = s1/nx;
  sd = sqrt((s2 - nx * mean * mean) / (nx - 1.0));
  se = sd/meanedf;
  t = mean/se;
  fprintf(stderr,"        middle-deep:  %8.4f   %8.4f   %8.4f   %9.3f\n",mean,sd,se,t);

  
  s1=s2=0;
  for (i=0; i<npix; i++) {
    d = p1[i]-p2[i];
    s1 += d;
    s2 += d*d;
  }
  mean = s1/nx;
  sd = sqrt((s2 - nx * mean * mean) / (nx - 1.0));
  se = sd/meanedf;
  t = mean/se;
  fprintf(stderr," middle-superficial:  %8.4f   %8.4f   %8.4f   %9.3f\n",mean,sd,se,t);

  s1=s2=0;
  for (i=0; i<npix; i++) {
    d = p2[i]-p0[i];
    s1 += d;
    s2 += d*d;
  }
  mean = s1/nx;
  sd = sqrt((s2 - nx * mean * mean) / (nx - 1.0));
  se = sd/meanedf;
  t = mean/se;
  fprintf(stderr,"   superficial-deep:  %8.4f   %8.4f   %8.4f   %9.3f\n",mean,sd,se,t);

  fprintf(stderr,"\n");
  gsl_matrix_free(D);
}
