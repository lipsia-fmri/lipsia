
/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/Vlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>


/* From the standard C library: */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../cylutils/cyl.h"


int LaminarSmooth(VImage zmap,VImage metric,Cylinders *cyl,size_t cid,VImage dest,VImage ndest)
{
  size_t i,k,nzero=0;
  size_t n = cyl->addr[cid]->size;
  int b,r,c,j;
  double u,jx;
  
  gsl_vector *y = gsl_vector_calloc(n);
  gsl_vector *mvec = gsl_vector_calloc(n);

  for (i=0; i<n; i++) {
    k = cyl->addr[cid]->data[i];
    b = gsl_matrix_int_get(cyl->xmap,k,0);
    r = gsl_matrix_int_get(cyl->xmap,k,1);
    c = gsl_matrix_int_get(cyl->xmap,k,2);
    mvec->data[i] = VPixel(metric,b,r,c,VFloat);
    y->data[i] = VPixel(zmap,b,r,c,VFloat);
    if (fabs(y->data[i]) < TINY) nzero++;
  }
  if (nzero > n/2) {
    gsl_vector_free(mvec);
    gsl_vector_free(y);
    return -1;
  }
  
  int numlayers = 31;
  double nx = (double)numlayers;
  double sum=0,kx=0;
  for (j=0; j<numlayers; j++) {
    jx = (double)j;

    sum=kx=0;
    for (i=0; i<n; i++) {
      u = mvec->data[i];
      if (u > jx/nx && u <= (jx+1)/nx) {
	sum += y->data[i];
	kx++;
      }
    }

    for (i=0; i<n; i++) {
      u = mvec->data[i];
      if (u > jx/nx && u <= (jx+1)/nx) {
	if (kx > 0.001) y->data[i] = sum/kx;
      }
    }
  }

  

  for (i=0; i<n; i++) {
    k = cyl->addr[cid]->data[i];
    b = gsl_matrix_int_get(cyl->xmap,k,0);
    r = gsl_matrix_int_get(cyl->xmap,k,1);
    c = gsl_matrix_int_get(cyl->xmap,k,2);
    VPixel(dest,b,r,c,VFloat) += y->data[i];
    VPixel(ndest,b,r,c,VFloat) += 1.0;
  }
  
  gsl_vector_free(mvec);
  gsl_vector_free(y);
  return 1;
}

