
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/* Gaussian function */
double tgauss(double x,double sigma)
{
  double y,z,a=2.506628273;
  z = x / sigma;
  y = exp((double)-z*z*0.5)/(sigma * a);
  return y;
}



/* create a Gaussian kernel for convolution */
gsl_vector *CGaussKernel(double sigma)
{
  int i,dim,n;
  double x,u,step;
  gsl_vector *kernel;

  /* no temporal smoothing */
  if (sigma < 0.001) return NULL;

  dim  = 4.0 * sigma + 1;
  n    = 2*dim+1;
  step = 1;

  kernel = gsl_vector_alloc(n);
  
  x = -(float)dim;
  for (i=0; i<n; i++) {
    u = tgauss(x,sigma);
    gsl_vector_set(kernel,i,u);
    x += step;
  }
  return kernel;
}


/* fill gaussian matrix */
void GaussMatrix(gsl_matrix *S,gsl_vector *kernel)
{
  int dim = kernel->size/2;
  int i,j,k;
  gsl_matrix_set_zero(S);
  for (i=0; i<S->size1; i++) {
    double sum = 0;
    for (j=0; j<S->size2; j++) {
      k = i-j+dim;
      if (k < 0 || k >= kernel->size) continue;
      double x = gsl_vector_get(kernel,k);
      sum += x;
      gsl_matrix_set(S,i,j,x);
    }
    /* normalize */
    for (j=0; j<S->size2; j++) {
      double x = gsl_matrix_get(S,i,j);
      gsl_matrix_set(S,i,j,x/sum);
    }
  }
}


void VectorConvolve(const double *src,gsl_vector *dst,gsl_vector *kernel)
{
  int i,j,k,len,n,dim;
  double sum;

  dim = kernel->size/2;
  len = dst->size;
  n   = len - dim;

  for (i=dim; i<n; i++) {
    sum = 0;
    k=0;
    for (j=i-dim; j<=i+dim; j++) {
      sum += src[j] * kernel->data[k];
      k++;
    }
    dst->data[i] = sum;
  }


  /* Randbehandlung */
  for (i=0; i<dim; i++) {
    dst->data[i] = src[i];
  }
  for (i=n; i<len; i++) {
    dst->data[i] = src[i];
  }
}
