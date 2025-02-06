
/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/Vlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>


/* From the standard C library: */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../cylutils/cyl.h"

extern double Xt2z(double t,double df);


/* simple layer averaging */
void layermean(gsl_vector *y,gsl_vector *mvec,gsl_vector *mean,gsl_vector *var,gsl_vector *kx,double voxel_scale)
{
  size_t i;
  double s0=0,s1=0,s2=0,sd0=0,sd1=0,sd2=0,nx0=0,nx1=0,nx2=0;

  for (i=0; i<y->size; i++){
    if (mvec->data[i] < 0.333) {
      s0 += y->data[i]; sd0 += y->data[i]*y->data[i];  nx0++;
    }
    if (mvec->data[i] >= 0.333 && mvec->data[i] < 0.666) {
      s1 += y->data[i]; sd1 += y->data[i]*y->data[i]; nx1++;
    }
    if (mvec->data[i] >= 0.666) {
      s2 += y->data[i]; sd2 += y->data[i]*y->data[i];  nx2++;
    }
  }


  /* estimate original number of voxels prior to upsampling */
  gsl_vector_set_zero(mean);
  gsl_vector_set_zero(var);
  gsl_vector_set_zero(kx);

  nx0 /= voxel_scale;
  nx1 /= voxel_scale;
  nx2 /= voxel_scale;
  
  if (nx0 > 2.1) {
    mean->data[0] = s0/nx0;
    var->data[0] = (sd0 - nx0 * mean->data[0] * mean->data[0]) / (nx0 - 1.0);
    kx->data[0] = nx0;
  }
  if (nx1 > 2.1) {
    mean->data[1] = s1/nx1;
    var->data[1] = (sd1 - nx1 * mean->data[1] * mean->data[1]) / (nx1 - 1.0);
    kx->data[1] = nx1;
  }
  if (nx2 > 2.1) {
    mean->data[2] = s2/nx2;
    var->data[2] = (sd2 - nx2 * mean->data[2] * mean->data[2]) / (nx2 - 1.0);
    kx->data[2] = nx2;
  }
}


/* two-sample Welsh test */
double Welsh(double mean1,double mean2,double var1,double var2,double nx1,double nx2)
{
  if (nx1 < 5 || nx2 < 5) return 0;
 
  double s1 = (var1*var1)/(nx1*nx1*(nx1-1.0));
  double s2 = (var2*var2)/(nx2*nx2*(nx2-1.0));

  double var = (var1/nx1 + var2/nx2);
  double df = (var*var)/(s1 + s2);
  
  double t = (mean1 - mean2)/sqrt(var);
  double z = Xt2z(t,df);
  return z;
}

double LayerMean(VImage zmap,VImage metric,Cylinders *cyl,size_t cid,double voxel_scale,
		 gsl_vector *beta,gsl_vector *zval)
{
  int b,r,c;
  size_t i,k;
  size_t dim=beta->size;
  size_t n=cyl->addr[cid]->size;
  double rtcode = -1;

  gsl_vector_set_zero(beta);
  gsl_vector_set_zero(zval);

  
  /* cylinder must be big enough for sufficient stats */
  if (n < dim*10) return -1;
  
  gsl_vector *y = gsl_vector_calloc(n);
  gsl_vector *mvec = gsl_vector_calloc(n);

  
  /* fill y-vector with activation map values */
  for (i=0; i<n; i++) {
    k = cyl->addr[cid]->data[i];
    b = gsl_matrix_int_get(cyl->xmap,k,0);
    r = gsl_matrix_int_get(cyl->xmap,k,1);
    c = gsl_matrix_int_get(cyl->xmap,k,2);
    y->data[i] = VPixel(zmap,b,r,c,VFloat);
    mvec->data[i] = VPixel(metric,b,r,c,VFloat);
  }

  /* return if z-values in this cylinder have no variance */
  double tss = gsl_stats_tss(y->data,1,y->size);
  if (tss < 0.001) goto ende;


  gsl_vector *var = gsl_vector_calloc(3);
  gsl_vector *kx = gsl_vector_calloc(3);
    
  layermean(y,mvec,beta,var,kx,voxel_scale);
 
  zval->data[0] = Welsh(beta->data[0],beta->data[1],var->data[0],var->data[1],
			kx->data[0],kx->data[1]);
 
  zval->data[1] = Welsh(beta->data[0],beta->data[2],var->data[0],var->data[2],
			kx->data[0],kx->data[2]);
 
  zval->data[2] = Welsh(beta->data[1],beta->data[2],var->data[1],var->data[2],
			kx->data[1],kx->data[2]);
 
  gsl_vector_free(kx);
  gsl_vector_free(var);
  rtcode = 1;
 
  
  /* free memory */
 ende: ;
  gsl_vector_free(y);
  gsl_vector_free(mvec);
  return rtcode;
}



