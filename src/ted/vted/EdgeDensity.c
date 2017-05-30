/*
** compute edge densities
**
** G.Lohmann, Feb 2015
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

extern void GetRank(float *data,gsl_vector *v,gsl_permutation *perm,gsl_permutation *rank);
extern float EdgeCorr(const float *data1,const float *data2,int,int);


int GetAddr(VImage mapimage,int bi,int ri,int ci,int m,int k,int l)
{
  int b = bi+m;
  int r = ri+k;
  int c = ci+l;
  if (b < 0 || r < 0 || c < 0) return -1;
  if (b >= VImageNBands(mapimage) || r >= VImageNRows(mapimage) || c >= VImageNColumns(mapimage)) return -1;
  int ii = (long)VPixel(mapimage,b,r,c,VInteger);
  return ii;
}

int VEdgeNeigbours(size_t id,VImage map,VImage mapimage,int adjdef,int *x)
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
      	int jj = GetAddr(mapimage,bi,ri,ci,m,k,l);
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
      	x[n] = jj;
      	n++;
      }
    }
  }
  return n;
}

void VCheckImage(VImage src)
{
  VAttrList out_list = VCreateAttrList();
  VAppendAttr(out_list,"image",NULL,VImageRepn,src);
  FILE *out_file = fopen("test.v","w");
  VWriteFile (out_file, out_list);
}


size_t EdgeDensity(float *E,int *I,int *J,size_t nedges_estimated,
		   gsl_histogram *TedHist,gsl_matrix_float *SNR1,gsl_matrix_float *SNR2,int n1,int n2,
		   VImage roi,VImage map,VImage mapimage,int adjdef,float elength,float zthreshold,
		   int nperm,int numperm,int step,float noise_cutoff,int metric)
{
  size_t i;
  size_t nvox=SNR1->size1;
  int nt=SNR1->size2;
  size_t progress=0;
  int rad2 = (int)(elength*elength);
  float tiny=1.0e-6;
  gsl_set_error_handler_off ();

  int iisel=4;

  /* size of neighbourhood */
  int nadj=0;
  if (adjdef == 0) nadj = 7;
  else if (adjdef == 1) nadj = 19;
  else if (adjdef == 2) nadj = 27;
  else if (adjdef == 3) nadj = 33;
  else VError(" illegal adjdef %d",adjdef);
  size_t nedges=0;
  int minadj = 1;

  /* histogram init */
  size_t nbins = gsl_histogram_bins (TedHist);
  double hmax = gsl_histogram_max (TedHist);
  double hmin = gsl_histogram_min (TedHist);


#pragma omp parallel for shared(nedges,progress,TedHist,E,I,J) schedule(dynamic)
  for (i=0; i<nvox; i+=step) {
    if (i%1000 == 0) fprintf(stderr," %ld000\r",(long)progress++);

    int bi = (int)VPixel(map,0,0,i,VShort);
    int ri = (int)VPixel(map,0,1,i,VShort);
    int ci = (int)VPixel(map,0,2,i,VShort);

    int iselect=0;

    int roiflagi = 0;
    if (roi != NULL) {
      if (VGetPixel(roi,bi,ri,ci) > 0.5) roiflagi = 1;
    }

    const float *datax1 = gsl_matrix_float_const_ptr(SNR1,i,0);
    const float *datax2 = gsl_matrix_float_const_ptr(SNR2,i,0);

    int *x = (int *) VCalloc(nadj,sizeof(int));
    int *y = (int *) VCalloc(nadj,sizeof(int));

    gsl_histogram *tmphist = gsl_histogram_alloc (nbins);
    gsl_histogram_set_ranges_uniform (tmphist,hmin,hmax);
    gsl_histogram_reset(tmphist);

    int nadjx = VEdgeNeigbours(i,map,mapimage,adjdef,x);
    if (nadjx < minadj) continue;

    size_t j=0;
    for (j=0; j<i; j+=step) {
      int bj = (int)VPixel(map,0,0,j,VShort);
      int rj = (int)VPixel(map,0,1,j,VShort);
      int cj = (int)VPixel(map,0,2,j,VShort);
      int d = SQR(bi-bj) + SQR(ri-rj) + SQR(ci-cj);
      if (d < rad2) continue;

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
      if (z < zthreshold) continue;
 
      int nadjy = VEdgeNeigbours(j,map,mapimage,adjdef,y);
      if (nadjy < minadj) continue;


      /* use every 'step' value, only when estimating noise cutoff, i.e. step > 1 */
      if (iselect%iisel == 0 && step > 1) {
	iselect = 0;
      }
      else if (iselect%iisel != 0 && step > 1) {
	iselect++;
	continue;
      }


      /* inspect local neighbourhood */
      int r,s;
      double kx=0.0,nx=0.0;
      for (s=0; s<nadjx; s++) {
      	size_t ii = (size_t) x[s];
        int bii = (int)VPixel(map,0,0,ii,VShort);
        int rii = (int)VPixel(map,0,1,ii,VShort);
        int cii = (int)VPixel(map,0,2,ii,VShort);

	const float *dataxx1 = gsl_matrix_float_const_ptr(SNR1,ii,0);
	const float *dataxx2 = gsl_matrix_float_const_ptr(SNR2,ii,0);

      	for (r=0; r<nadjy; r++) {
      	  size_t jj = (size_t) y[r];
          int bjj = (int)VPixel(map,0,0,jj,VShort);
          int rjj = (int)VPixel(map,0,1,jj,VShort);
          int cjj = (int)VPixel(map,0,2,jj,VShort);

          int d2 = SQR(bii-bjj) + SQR(rii-rjj) + SQR(cii-cjj);
          if (d2 < rad2) continue;

	  const float *datayy1 = gsl_matrix_float_const_ptr(SNR1,jj,0);
	  const float *datayy2 = gsl_matrix_float_const_ptr(SNR2,jj,0);

	  /* edge z-value */
	  double z1 = EdgeCorr(dataxx1,datayy1,nt,metric);
	  double z2 = EdgeCorr(dataxx2,datayy2,nt,metric);
	  double z = (z1-z2);
          if (ABS(z) < tiny) continue;
      	  if (z > zthreshold) kx++;
          nx++;
     	}
      }
      if (nx < 1.0) continue;
      double ted = kx/nx;

      if (ted < hmin) ted = hmin;
      if (ted > hmax-tiny) ted = hmax-tiny;
      gsl_histogram_increment (tmphist,ted);

      if (E == NULL) continue;
      if (nperm != numperm) continue;
      if (ted < noise_cutoff) continue;


      #pragma omp critical 
      {
	if (nedges >= nedges_estimated) {
	  VError(" alloc, nedges_estimated:  %lu",nedges_estimated);
	}
	E[nedges] = ted;    /* task-based edge density */
	I[nedges] = i;      /* voxel address i */
	J[nedges] = j;      /* voxel address j */
	nedges++;
      }
    }
#pragma omp critical 
    {
      gsl_histogram_add (TedHist,tmphist);
    }
    gsl_histogram_free (tmphist);
    VFree(x);
    VFree(y);
  }
  return nedges;
}


/* estimate number of truly needed edges */
size_t EstimateEdges(float *E,int *I,int *J,size_t old_estimate,
		     gsl_histogram *TedHist,gsl_matrix_float *SNR1,gsl_matrix_float *SNR2,int n1,int n2,
		     VImage roi,VImage map,VImage mapimage,int adjdef,float elength,float zthr,
		     float noise_cutoff,int metric)
{
  float tmp_cutoff = -1.0;
  int step = 2;
  int nperm=0,numperm=1;

  gsl_histogram_reset(TedHist);
  EdgeDensity(E,I,J,old_estimate,TedHist,SNR1,SNR2,n1,n2,roi,map,mapimage,(int)adjdef,elength,zthr,
	      nperm,numperm,step,tmp_cutoff,metric);

  double sd = gsl_histogram_sigma (TedHist);
  if (sd < 1.0e-6) VError(" stdev of TedHist: %f",sd);

  size_t i0 = 0;
  if (gsl_histogram_find (TedHist,(double)noise_cutoff,&i0) != GSL_SUCCESS)
    VError (" err finding noise cutoff index, %f",noise_cutoff);
  
  size_t hbins = gsl_histogram_bins(TedHist);
  gsl_histogram_pdf *cdf = gsl_histogram_pdf_alloc(hbins);
  gsl_histogram_pdf_init (cdf,TedHist);


  double qx = 1.0 - cdf->sum[i0];
  if (qx > 1) qx = 1.0;
  size_t new_estimate = (size_t) (qx*(double)old_estimate);

  gsl_histogram_reset(TedHist);
  gsl_histogram_pdf_free (cdf);
  return new_estimate;
}
