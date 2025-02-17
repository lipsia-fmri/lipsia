#include <viaio/VImage.h>
#include <viaio/Vlib.h>

/* From the standard C library: */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>

#include "../cylutils/cyl.h"


/* Estimate spatial autocorrelation using Moran's I */
double Moran(Cylinders *cyl,size_t cid,gsl_vector *reso,gsl_vector *residuals)
{
  int k;
  double bi,ri,ci,bj,rj,cj;
  size_t i,j,n=residuals->size;
  double d,w,xi,xj,sum1,sum2,sumw,rho=0;
  
  /* voxel addresses scaled by voxel size */
  double *bx = (double *)VCalloc(n,sizeof(double));
  double *rx = (double *)VCalloc(n,sizeof(double));
  double *cx = (double *)VCalloc(n,sizeof(double));
  for (i=0; i<n; i++) {
    k = cyl->addr[cid]->data[i];
    bx[i] = ((double)gsl_matrix_int_get(cyl->xmap,k,0))*reso->data[0];
    rx[i] = ((double)gsl_matrix_int_get(cyl->xmap,k,1))*reso->data[1];
    cx[i] = ((double)gsl_matrix_int_get(cyl->xmap,k,2))*reso->data[2];
  }

  /* mean residuals */
  double mean = gsl_stats_mean(residuals->data,1,residuals->size);

  /* for all voxels in this cylinder */
  sum1=sum2=sumw=0;
  for (i=0; i<n; i++) {
    bi = bx[i];
    ri = rx[i];
    ci = cx[i];
    xi = residuals->data[i] - mean;
    sum2 += xi*xi;
    
    for (j=0; j<i; j++) {  /* symmetric */
      bj = bx[j];
      rj = rx[j];
      cj = cx[j];
      d = (bi-bj)*(bi-bj) + (ri-rj)*(ri-rj) + (ci-cj)*(ci-cj);
      w = exp(-d/2.0);   /* assume spatial gaussian decay */
      sumw += 2.0*w;
      
      xj = residuals->data[j] - mean;
      sum1 += 2.0*w*xi*xj;
    }
  }
  VFree(bx);
  VFree(rx);
  VFree(cx);
  if (sum2 < TINY) return 0;
  
  sumw /= (double)n;
  rho = sum1/sum2;
  rho /= sumw;
  if (rho > 1) rho = 1;
  return rho;
}
