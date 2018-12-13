/*
**  Matrix Projection,
**   Ref:  Lohmann, Loktyushin, Stelzer, Scheffler (2018) bioRXiv
**
** G.Lohmann, MPI-KYB, Nov 2018
*/

#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/

#define SQR(x) ((x) * (x))
#define IMIN(x,y) ((x) < (y) ? (x) : (y))
#define IMAX(x,y) ((x) > (y) ? (x) : (y))

extern double ECMcorrelation(const float *arr1,const float *arr2,size_t nt,int type);
extern void NormVec(float *x,size_t);


void VMatrixProjection(gsl_matrix_float *X,float *ev,int type,int seed)
{
  size_t i,j;
  size_t nvox = X->size1;
  size_t nt = X->size2;
  size_t p = 32;

  fprintf(stderr," VMatrixProjection...\n");
  
  /* generate random matrix */
  const gsl_rng_type *Tx = gsl_rng_default;
  gsl_rng *rx = gsl_rng_alloc (Tx);
  gsl_rng_set(rx,(unsigned long int)seed);
  gsl_matrix *R = gsl_matrix_calloc(nvox,p);
  for (i=0; i<nvox; i++) {
    for (j=0; j<p; j++) {
      gsl_matrix_set(R,i,j,gsl_ran_ugaussian(rx));
    }
  }

  /* matrix product 1 */
  gsl_matrix *Y = gsl_matrix_calloc(nvox,p);
  size_t progress=0;
  size_t step=4;  /* sparse */
  fprintf(stderr," pass 1:\n");


#pragma omp parallel for shared(progress) schedule(guided) firstprivate(X,Y,R)
  for (i=0; i<nvox; i++) {
    if (i%1000 == 0) fprintf(stderr," %d000  of  %lu\r",(int)(++progress),nvox);
    
    size_t j=0,k=0;
    double *tmp1 = (double *)VCalloc(nvox,sizeof(double));
    double *tmp2 = NULL;
    double sum=0;
    const float *arr1 = gsl_matrix_float_const_ptr(X,i,0);

    for (j=0; j<nvox; j++) {
      const float *arr2 = gsl_matrix_float_const_ptr(X,j,0);
      tmp1[j] = ECMcorrelation(arr1,arr2,nt,type);
    }

    for (j=0; j<p; j++) {
      tmp2 = gsl_matrix_ptr(R,j,0);
      sum=0;
      for (k=0; k<nvox; k+=step) {
	sum += tmp1[k]*tmp2[k];
      }
      gsl_matrix_set(Y,i,j,sum);
    }
    VFree(tmp1);
  }
  fprintf(stderr,"\n");
  gsl_matrix_free(R);
  
 
  /* get orthonormal basis for Y  (SVD) */
  gsl_vector *sv = gsl_vector_calloc(p);
  gsl_vector *work = gsl_vector_calloc(p);
  gsl_matrix *V = gsl_matrix_calloc(p,p);
  gsl_linalg_SV_decomp(Y,V,sv,work);


  /* matrix product 2 */
  fprintf(stderr," pass 2:\n");
  gsl_matrix *D = gsl_matrix_calloc(nvox,p);
  progress=0;

  #pragma omp parallel for shared(progress) schedule(guided) firstprivate(X,Y,D)
  for (i=0; i<nvox; i++) {
    if (i%1000 == 0) fprintf(stderr," %d000  of  %lu\r",(int)(++progress),nvox);

    size_t j=0,k=0;
    double *tmp1 = (double *)VCalloc(nvox,sizeof(double));
    double sum=0;
    const float *arr1 = gsl_matrix_float_const_ptr(X,i,0);

    for (j=0; j<nvox; j++) {
      const float *arr2 = gsl_matrix_float_const_ptr(X,j,0);
      tmp1[j] = ECMcorrelation(arr1,arr2,nt,type);
    }

    for (j=0; j<p; j++) {
      sum=0;      
      for (k=0; k<nvox; k+=step) {
	sum += tmp1[k] * gsl_matrix_get(Y,k,j);
      }
      gsl_matrix_set(D,i,j,sum);
    }
    VFree(tmp1);
  }
  fprintf(stderr,"\n");
  gsl_matrix_free(Y);

  
  /* final SVD */
  gsl_linalg_SV_decomp(D,V,sv,work);

  
  /* output */
  size_t k=0;
  for (i=0; i<nvox; i++) {
    ev[i] = gsl_matrix_get(D,i,0);
    if (ev[i] < 0) k++;
  }

  /* sign switch if needed */
  if (k > nvox/2) {
    for (i=0; i<nvox; i++) {
      ev[i] = -ev[i];
    }
  }

  /* normalize */
  NormVec(ev,nvox);
  float kx = sqrt((float)(nvox));
  for (i=0; i<nvox; i++)  ev[i] *= kx;
}
