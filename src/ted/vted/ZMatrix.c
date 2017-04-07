/*
** compute correlation matrix
**
** G.Lohmann, Feb 2011
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_histogram.h>

#include "viaio/Vlib.h"
#include "viaio/VImage.h"
#include "viaio/mu.h"

#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/

#define SQR(x) ((x) * (x))
#define ABS(x) ((x) > 0 ? (x) : -(x))

extern void VNormalize(float *data,int nt,VBoolean stddev);
extern long GetAddr(VImage mapimage,int bi,int ri,int ci,int m,int k,int l);

int VNumNeigbours(size_t id,VImage map,VImage mapimage,int adjdef)
{
  int nslices = VImageNBands(mapimage);
  int nrows = VImageNRows(mapimage);
  int ncols = VImageNColumns(mapimage);

  int bi = VPixel(map,0,0,id,VShort);
  int ri = VPixel(map,0,1,id,VShort);
  int ci = VPixel(map,0,2,id,VShort);
  if (bi < 1 || bi >= nslices-2) return 0;
  if (ri < 1 || ri >= nrows-2) return 0;
  if (ci < 1 || ci >= ncols-2) return 0;

  int wn=1,rad2=0;
  if (adjdef == 3) {
    wn = 2;
    rad2 = 2*2;
  }

  int n=0;
  int k,l,m;
  for (m=-wn; m<=wn; m++) {
    for (k=-wn; k<=wn; k++) {
      for (l=-wn; l<=wn; l++) {
      	long jj = GetAddr(mapimage,bi,ri,ci,m,k,l);
      	if (jj < 0) continue;  /* outside of brain mask */

      	if (adjdef == 0) {     /* 6 adjacency */
      	  if (ABS(m)+ABS(k)+ABS(l) > 1) continue;
      	}
      	if (adjdef == 1) {     /* 18-adjacency */
      	  if (ABS(m) > 0 && ABS(k) > 0 && ABS(l) > 0) continue;
      	}	
	if (adjdef == 3) {     /* sphere */
	  if (m*m + k*k + l*l > rad2) continue;
	}
      	n++;
      }
    }
  }
  return n;
}

/* convert to ranks */
void GetRank(float *data,gsl_vector *v,gsl_permutation *perm,gsl_permutation *rank)
{
  size_t i;
  size_t n = v->size;
  for (i=0; i<n; i++) gsl_vector_set(v,i,(double)data[i]);
  gsl_sort_vector_index (perm, v);
  gsl_permutation_inverse (rank, perm);
  for (i=0; i<n; i++) data[i] = (float)rank->data[i];
}

double Spearman(const float *data1,const float *data2,int n)
{
  int i;
  double nx = (double)n;
  double kx = nx*(nx*nx-1.0);
  double sxy=0.0;
  for (i=0; i<n; i++) {
    const double u = (double)data1[i];
    const double v = (double)data2[i];
    const double d = (u-v);
    sxy += d*d;
  }
  double rho = 1.0 - 6.0*sxy/kx;
  return rho;
}


double Correlation(const float *data1,const float *data2,int n)
{
  int i;
  double corr=0;
  for (i=0; i<n; i++) {   
    const double u = (double)data1[i];
    const double v = (double)data2[i];
    corr += u*v;
  }
  corr /= (double)n;
  return corr;
}


float EdgeCorr(const float *data1,const float *data2,int n,int metric)
{
  float corr=0.0;
  if (metric == 0) corr = Correlation(data1,data2,n);
  if (metric == 1) corr = Spearman(data1,data2,n);
  if (corr > 0) return corr;
  else return 0.0;
}



void GetSNR(gsl_matrix_float **X1,gsl_matrix_float **X2,int *table,int n,gsl_matrix_float *SNR,int metric)
{
  long j,nvox=X1[0]->size1;
  long k,s,nt=X1[0]->size2;

  if (n < 3) VError(" n: %d",n);
  double nx = (double)n;
  double *sum1 = (double *)VCalloc(nt,sizeof(double));
  double *sum2 = (double *)VCalloc(nt,sizeof(double));

  gsl_vector *vec = NULL;
  gsl_permutation *perm = NULL;
  gsl_permutation *rank = NULL;
  if (metric == 1) {  /* only needed for spearman correlation */
    vec = gsl_vector_calloc(nt);
    perm = gsl_permutation_alloc(nt);
    rank = gsl_permutation_alloc(nt);
  }

  double ave=0,var=0,snr=0;

  for (j=0; j<nvox; j++) {
    memset(sum1,0,nt*sizeof(double));
    memset(sum2,0,nt*sizeof(double));
    const float *pp=NULL;
    for (s=0; s<n; s++) {
      if (table[s] == 0)
	pp = gsl_matrix_float_const_ptr(X1[s],j,0);
      else
	pp = gsl_matrix_float_const_ptr(X2[s],j,0);
      for (k=0; k<nt; k++) {
	const double x = (double)(*pp++);
	sum1[k] += x;
	sum2[k] += x*x;
      }
    }
    for (k=0; k<nt; k++) {
      ave = (sum1[k]/nx);
      var = (sum2[k] - nx*ave*ave) / (nx - 1.0);
      snr = 0;
      if (var > 0) snr = ave / sqrt(var);
      gsl_matrix_float_set(SNR,j,k,snr);
    }

    if (metric == 0) {  /* pearson correlation */
      float *qq = gsl_matrix_float_ptr(SNR,j,0);
      VNormalize(qq,nt,TRUE);
    }
    if (metric == 1) {  /* spearman ranks      */
      float *qq = gsl_matrix_float_ptr(SNR,j,0);
      GetRank(qq,vec,perm,rank);
    }
  }
  VFree(sum1);
  VFree(sum2);
}




float ZMatrix(gsl_histogram *histogram,gsl_matrix_float *SNR1,gsl_matrix_float *SNR2,int n1,int n2,
	      VImage roi,VImage map,VImage mapimage,int adjdef,float elength,float quantile,int step,int metric)
{
  size_t i;
  size_t nvox=SNR1->size1;
  int nt=SNR1->size2;
  size_t progress=0;
  int rad2 = (int)(elength*elength);
  double tiny=1.0e-8;
  gsl_set_error_handler_off ();


  fprintf(stderr," Computing matrix...\n");

  gsl_histogram_reset(histogram);
  size_t nbins = gsl_histogram_bins (histogram);
  double hmax = gsl_histogram_max (histogram);
  double hmin = gsl_histogram_min (histogram);
  double zmin = 99999.0;
  double zmax = -99999.0;
  int minadj = 1;


#pragma omp parallel for shared(progress,histogram) schedule(dynamic) firstprivate(SNR1,SNR2)
  for (i=0; i<nvox; i+=step) {
    if (i%1000 == 0) fprintf(stderr," %ld000\r",(long)progress++);
    int bi = (int)VPixel(map,0,0,i,VShort);
    int ri = (int)VPixel(map,0,1,i,VShort);
    int ci = (int)VPixel(map,0,2,i,VShort);
    int nadjx = VNumNeigbours(i,map,mapimage,adjdef);
    if (nadjx < minadj) continue;

    int roiflagi = 0;
    if (roi != NULL) {
      if (VGetPixel(roi,bi,ri,ci) > 0.5) roiflagi = 1;
    }

    gsl_histogram *tmphist = gsl_histogram_alloc (nbins);
    gsl_histogram_set_ranges_uniform (tmphist,hmin,hmax);
    gsl_histogram_reset(tmphist);
    
    const float *datax1 = gsl_matrix_float_const_ptr(SNR1,i,0);
    const float *datax2 = gsl_matrix_float_const_ptr(SNR2,i,0);


    size_t j=0;
    for (j=0; j<i; j+=step) {
      int bj = (int)VPixel(map,0,0,j,VShort);
      int rj = (int)VPixel(map,0,1,j,VShort);
      int cj = (int)VPixel(map,0,2,j,VShort);
      int d = SQR(bi-bj) + SQR(ri-rj) + SQR(ci-cj);
      if (d < rad2) continue;
      int nadjy = VNumNeigbours(j,map,mapimage,adjdef);
      if (nadjy < minadj) continue;

      int roiflagj = 0;
      if (roi != NULL) {
	if (VGetPixel(roi,bj,rj,cj) > 0.5) roiflagj = 1;
	if (roiflagi + roiflagj != 1) continue; 
      }
      const float *datay1 = gsl_matrix_float_const_ptr(SNR1,j,0);
      const float *datay2 = gsl_matrix_float_const_ptr(SNR2,j,0);

      /* edge z-value */
      double z1 = EdgeCorr(datax1,datay1,nt,metric);
      double z2 = EdgeCorr(datax2,datay2,nt,metric);
      double z = (z1-z2);
      if (z < zmin) zmin = z;
      if (z > zmax) zmax = z;

      if (z < hmin) z = hmin;
      if (z > hmax-tiny) z = hmax-tiny;
      gsl_histogram_increment (tmphist,z);
    }


#pragma omp critical 
    {
      gsl_histogram_add (histogram,tmphist);
    }
    gsl_histogram_free (tmphist);
  }


  /* get quantile cutoff */
  gsl_histogram_pdf *pdf = gsl_histogram_pdf_alloc(nbins);
  gsl_histogram_pdf_init(pdf,histogram);
  double lower=0,upper=0;

  size_t i0=0;
  for (i=nbins-1; i>=0; i--) {
    gsl_histogram_get_range (histogram,i,&lower,&upper);
    if (pdf->sum[i] < quantile) {
      if (gsl_histogram_get_range (histogram,i,&lower,&upper) == GSL_EDOM) VError(" err hist");
      i0 = i;
      break;
    }
  }
  double sum=0;
  for (i=nbins-1; i>i0; i--) {
    sum += gsl_histogram_get(histogram,i);
  }
  return (float)upper;
}
