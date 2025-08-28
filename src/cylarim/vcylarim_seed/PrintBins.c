
/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/Vlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_chebyshev.h>

/* From the standard C library: */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../cylutils/cyl.h"

extern double Xt2z(double t,double df);
extern void YNorm(gsl_vector *x);



void PrintBins(VImage zmap,VImage metric,VImage dest,Cylinders *cyl,size_t cid,char *filename)
{
  int b,r,c;
  size_t n0 = cyl->addr[cid]->size;
  size_t i,j,k;
  double u=0;

  gsl_vector *y0 = gsl_vector_calloc(n0);
  gsl_vector *mvec0 = gsl_vector_calloc(n0);
  for (i=0; i<n0; i++) {
    k = cyl->addr[cid]->data[i];
    b = gsl_matrix_int_get(cyl->xmap,k,0);
    r = gsl_matrix_int_get(cyl->xmap,k,1);
    c = gsl_matrix_int_get(cyl->xmap,k,2);
    mvec0->data[i] = VPixel(metric,b,r,c,VFloat);
    y0->data[i] = VPixel(zmap,b,r,c,VFloat);

    if (mvec0->data[i] > 0.05 && mvec0->data[i] < 0.95)
      VPixel(dest,b,r,c,VFloat) = y0->data[i];
  }
 
  
  size_t n = 0;
  for (i=0; i<mvec0->size; i++) {
    if (mvec0->data[i] > 0.05 && mvec0->data[i] < 0.95) n++;
  }
  gsl_vector *y = gsl_vector_calloc(n);
  gsl_vector *mvec = gsl_vector_calloc(n);

  j=0;
  for (i=0; i<mvec0->size; i++) {
    if (mvec0->data[i] > 0.05 && mvec0->data[i] < 0.95) {
      y->data[j] = y0->data[i];
      mvec->data[j] = mvec0->data[i];
      j++;
    }
  }

  FILE *fp0 = fopen("deep.txt","w");
  for (i=0; i<n; i++) {
    u = mvec->data[i];
    if (u < 0.333) fprintf(fp0," %f %f\n",u,y->data[i]);
  }
  fclose(fp0);

  FILE *fp1 = fopen("middle.txt","w");
  for (i=0; i<n; i++) {
    u = mvec->data[i];
    if (u >= 0.333 && u < 0.666) fprintf(fp1," %f %f\n",u,y->data[i]);
  }
  fclose(fp1);

  FILE *fp2 = fopen("superficial.txt","w");
  for (i=0; i<n; i++) {
    u = mvec->data[i];
    if (u > 0.666) fprintf(fp2," %f %f\n",u,y->data[i]);
  }
  fclose(fp2);
  
  gsl_vector_free(y);
  gsl_vector_free(mvec);
}
