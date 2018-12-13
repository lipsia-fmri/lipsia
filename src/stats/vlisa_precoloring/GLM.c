/*
** 1st level statistical inference using LISA
** GLM (general linear modeling) using precoloring
**
** G.Lohmann, MPI-KYB, 2018
*/
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>


extern gsl_matrix *PseudoInv(gsl_matrix *A,gsl_matrix *B);
extern void printmat(gsl_matrix *R,char *str);
extern void printvec(gsl_vector *x,char *str);


/*
** general linear model (GLM) using pecoloring
*/
void VGLM(gsl_matrix *Data,gsl_matrix *X,gsl_matrix *XInv,gsl_vector *con,VImage map,VImage zmap)
{
  int i;
  int m = Data->size2;
  int n = con->size;
  gsl_set_error_handler_off();
  gsl_vector *y = gsl_vector_calloc (m);
  gsl_vector *beta = gsl_vector_calloc (n);

 
  /* compute pseudoinverse */
  XInv = PseudoInv(X,XInv);


  /* main loop */
  VFillImage(zmap,VAllBands,0);
  size_t nvox=0;
  for (nvox=0; nvox < Data->size1; nvox++) {

    y->data = gsl_matrix_ptr(Data,nvox,0);
    gsl_blas_dgemv(CblasNoTrans,1.0,XInv,y,0.0,beta);

    /* contrast image */
    double sum=0;
    for (i=0; i<beta->size; i++) {
      sum += (beta->data[i]*con->data[i]);
    }
    if (gsl_isnan(sum) || gsl_isinf(sum)) continue;
    int b = VPixel(map,0,0,nvox,VShort);
    int r = VPixel(map,0,1,nvox,VShort);
    int c = VPixel(map,0,2,nvox,VShort);
    VPixel(zmap,b,r,c,VFloat) = sum;
  }

  gsl_vector_free(beta);
  gsl_vector_free(y);
}

