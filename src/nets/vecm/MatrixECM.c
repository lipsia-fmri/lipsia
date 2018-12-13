/*
** ECM - eigenvector centrality mapping using full matrix
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
#include <gsl/gsl_blas.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/

#define SQR(x) ((x) * (x))


/* check if matrix is irreducible. 
** If it is not, then Perron-Frobenius theorem does not apply
*/
void CheckReducibility(float *A,size_t nvox)
{
  size_t i,j,k,n;
  float umin = 99999;
  
  n=0;
  for (i=0; i<nvox; i++) {
    for (j=0; j<i; j++) {
      k=j+i*(i+1)/2;
      if (A[k] > 0) n++;
      if (A[k] > 0 && A[k] < umin) umin = A[k];
    }
    for (j=i+1; j<nvox; j++) {
      k=i+j*(j+1)/2;
      if (A[k] > 0) n++;
    }
  }
  
  if (n < nvox) {
    VWarning(" matrix is not irreducible, correction term is applied.");
    for (i=0; i<nvox; i++) {
      for (j=0; j<i; j++) {
	k=j+i*(i+1)/2;
	if (A[k] < umin) A[k] = umin*0.5;
      }
    }
  }
}



/* re-implementation of cblas_sspmv,  cblas_sspmv causes problems */
void my_sspmv(float *A,float *x,float *y,size_t n)
{
  size_t i,j,k,kk;
  float tmp1=0,tmp2=0;

  kk = k = 0;
  for (j=0; j<n; j++) {
    y[j] = 0;
    tmp1 = x[j];
    tmp2 = 0;
    k = kk;
    for (i=0; i<j; i++) {
      y[i] += tmp1 * A[k];
      tmp2 += A[k] * x[i];
      k++;
    }
    y[j] += tmp1*A[kk+j] + tmp2;
    kk += (j+1);
  }
}


void MatrixPowerIteration(float *A,float *ev,size_t nvox,int maxiter)
{
  int i,iter;
  float sum,d,nx;

  fprintf(stderr," power iteration...\n");
  float *y = (float *) VCalloc(nvox,sizeof(float));
  nx = (double)nvox;
  for (i=0; i<nvox; i++) ev[i] = y[i] = 1.0/nx;


  for (iter=0; iter < maxiter; iter++) {

    /* y = Ax,  A symmetric and lower-triangular */
    /* cblas_sspmv(CblasRowMajor,CblasLower,(int)n,1.0f,A,ev,1,1.0f,y,1); */
    my_sspmv(A,ev,y,nvox);


    /* normalize */
    sum = 0;
    for (i=0; i<nvox; i++) sum += y[i]*y[i];
    sum = sqrt(sum);

    /* check convergence */
    d = 0;
    for (i=0; i<nvox; i++) {
      y[i] /= sum;
      d += SQR(ev[i] - y[i]);
      ev[i] = y[i];
    }
    fprintf(stderr," %5d   %f\n",(int)iter,d);

    if (iter > 2 && d < 1.0e-6) break;
  }
  VFree(y);
}



float CorrMetric(float z,int type)
{
  switch (type) {
  case 1:  /* add */
    z += 1.0;
    break;
  case 2:  /* pos  */
    if (z < 0) z = 0;
    break;
  case 3:  /* abs */
    z = fabs(z);
    break;
  case 4:  /* neg */
    if (z < 0) z = -z;
    else z = 0;
    break;
  default:
    VError("illegal type");
  }
  return z;
}


double ECMcorrelation(const float *arr1,const float *arr2,size_t nt,int type)
{
  size_t i;
  double sum=0,z=0,kx=(double)nt;

  if ((type > 0 && type < 5) || (type==6)) {
    sum=0;
    for (i=0; i<nt; i++) {
      const double u = (double)arr1[i];
      const double v = (double)arr2[i];
      sum += u*v;
    }
    z = sum/kx;
  }

  switch (type) {

  case 0:  /* RLC positive */
    sum=0;
    for (i=0; i<nt; i++) {
      const double u = (double)arr1[i];
      const double v = (double)arr2[i];
      sum += u*v + fabs(u*v);
    }
    z = sum/(2.0*kx);
    break;

  case 1:  /* add */
     z += 1.0;
     break;
     
  case 2:  /* pos  */
    if (z < 0) z = 0;
    break;

  case 3:  /* abs */
    z = fabs(z);
    break;

  case 4:  /* neg */
    if (z < 0) z = -z;
    else z = 0;
    break;

  case 5:   /* Gaussian of Euclidean distance */
    sum=0;
    for (i=0; i<nt; i++) {
      const double u = (double)arr1[i];
      const double v = (double)arr2[i];
      const double d = u-v;
      sum += d*d;
    }
    z = sum/kx;
    z = exp(-0.5*z*z);
    break;

  case 7:  /* RLC negative */
    sum=0;
    for (i=0; i<nt; i++) {
      const double u = (double)arr1[i];
      const double v = (double)arr2[i];
      sum += (fabs(u*v) - u*v);
    }
    z = sum/(2.0*kx);
    break;


  default:
    ;
  }
  return z;
}



void VMatrixECM(gsl_matrix_float *X,float *ev,int type,int maxiter,VImage map)
{
  size_t i,j;
  size_t nvox = X->size1;
  size_t nt = X->size2;

  
  /* compute similarity matrix */
  size_t m = (nvox*(nvox+1))/2;
  float *A = (float *) calloc(m,sizeof(float));
  if (!A) VError(" err allocating correlation matrix (too big)");
  memset(A,0,m*sizeof(float));
  size_t progress=0;
  

#pragma omp parallel for shared(progress) private(j) schedule(guided) firstprivate(X,A)
  for (i=0; i<nvox; i++) {
    if (i%1000 == 0) fprintf(stderr," %d000  of  %lu\r",(int)(++progress),nvox);

    const float *arr1 = gsl_matrix_float_const_ptr(X,i,0);
    for (j=0; j<=i; j++) {
      if (i == j) continue;
      const float *arr2 = gsl_matrix_float_const_ptr(X,j,0);
      
      const double v = ECMcorrelation(arr1,arr2,nt,type);
      const size_t k=j+i*(i+1)/2;
      A[k] = v;
    }
  }
  fprintf(stderr,"\n");

  
  /* CheckReducibility(A,nvox); */
  /* DMN(map,A,nvox); */
    
  /* power iteration */
  MatrixPowerIteration(A,ev,nvox,maxiter);
  VFree(A);
}
