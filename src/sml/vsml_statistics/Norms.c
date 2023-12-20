/*
** Various methods for normalizaing data
**
** G.Lohmann, MPI-KYB, Jan 2020
*/
#include "viaio/Vlib.h"
#include "viaio/VImage.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_math.h>




extern void VCheckImage(VImage src);

extern double kth_smallest(double *a, size_t n, size_t k);
#define Median(a,n) kth_smallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))

#define SQR(x) ((x) * (x))
#define ABS(x) ((x) > 0 ? (x) : -(x))


void ZNorm(gsl_vector *x)
{
  double mean = gsl_stats_mean(x->data,1,x->size);
  double sd = gsl_stats_sd_m(x->data,1,x->size,mean);
  gsl_vector_add_constant(x,-mean);
  gsl_vector_scale(x,1.0/sd);
}

void VYNorm(gsl_vector *x)
{
  double mean = gsl_stats_mean(x->data,1,x->size);
  gsl_vector_add_constant(x,-mean);
}

void MeanSubtract(gsl_vector *x0,gsl_vector *x1)
{
  double mean = gsl_stats_mean(x0->data,1,x0->size);
  gsl_vector_add_constant(x0,-mean);
  if (x1 != NULL) gsl_vector_add_constant(x1,-mean);
}

void MeanAdd(gsl_vector *x0,gsl_vector *x1,double mean)
{
  gsl_vector_add_constant(x0,mean);
  if (x1 != NULL) gsl_vector_add_constant(x1,mean);
}



void XNorm(double *x,size_t n)
{
  size_t i;
  double mean = gsl_stats_mean(x,1,n);
  double sd = gsl_stats_sd_m(x,1,n,mean);
  for (i=0; i<n; i++) x[i] = (x[i]-mean)/sd;
}

void YNorm(double *x,size_t n)
{
  size_t i;
  double mean = gsl_stats_mean(x,1,n);
  for (i=0; i<n; i++) x[i] = (x[i]-mean);
}

/* rank norm, range [-1,1] */
void RNorm(double *x,size_t n)
{
  size_t j;
  size_t *perm = (size_t *) VCalloc(n,sizeof(size_t));
  double *wx = (double *)VCalloc(n,sizeof(double));  
  gsl_sort_index (perm,x,1,n);
  for (j=0; j<n; j++) wx[perm[j]] = (double)j/(double)n;
  for (j=0; j<n; j++) x[j] = 2.0*(wx[j]-0.5);
  VFree(wx);
  VFree(perm);
}

/* rank norm, range [0,1] */
void UNorm(double *x,size_t n)
{
  size_t j;
  size_t *perm = (size_t *) VCalloc(n,sizeof(size_t));
  double *wx = (double *)VCalloc(n,sizeof(double));  
  gsl_sort_index (perm,x,1,n);
  for (j=0; j<n; j++) wx[perm[j]] = (double)j/(double)n;
  for (j=0; j<n; j++) x[j] = wx[j];
  VFree(wx);
  VFree(perm);
}

/* rank norm, range [0,1] */
void RankNorm(gsl_vector *x)
{
  size_t i;
  double nx = (double)x->size;
  gsl_permutation *perm = gsl_permutation_alloc(x->size);
  gsl_permutation *rank = gsl_permutation_alloc(x->size);
  gsl_sort_vector_index (perm,x);
  gsl_permutation_inverse (rank, perm);
  for (i=0; i<x->size; i++) x->data[i] = rank->data[i];
  for (i=0; i<x->size; i++) x->data[i] /= nx;
  gsl_permutation_free(perm);
  gsl_permutation_free(rank);
}


/* double centering */
void KNorm(gsl_matrix *X)
{
  size_t i,j;
  double u,v,s,sij,ni=(double)X->size1,nj=(double)X->size2;
  gsl_matrix *A = gsl_matrix_calloc(X->size1,X->size2);
  gsl_matrix_memcpy(A,X);


  sij=0;
  for (i=0; i<X->size1; i++) {
    for (j=0; j<X->size2; j++) {
      sij += gsl_matrix_get(X,i,j);
    }
  }
  sij /= (ni*nj);

  double *tmpj = (double *) VCalloc(X->size2,sizeof(double));
  for (j=0; j<X->size2; j++) {
    s=0;
    for (i=0; i<X->size1; i++) s += gsl_matrix_get(X,i,j);
    tmpj[j] = s/ni;
  }

  double *tmpi = (double *) VCalloc(X->size1,sizeof(double));
  for (i=0; i<X->size1; i++) {
    s=0;
    for (j=0; j<X->size2; j++) s += gsl_matrix_get(X,i,j);
    tmpi[i] = s/nj;
  }
  
  for (i=0; i<X->size1; i++) {
    for (j=0; j<X->size2; j++) {
      u = gsl_matrix_get(X,i,j);
      v = u - tmpi[i] - tmpj[j] + sij;
      gsl_matrix_set(A,i,j,v);
    }
  }

  gsl_matrix_memcpy(X,A);
  gsl_matrix_free(A);
  VFree(tmpi);
  VFree(tmpj);
}


/* row, colum centering */
void ZKXNorm(gsl_matrix *Xtrain,gsl_matrix *Xtest)
{
  size_t i,j,dim=Xtrain->size1*Xtrain->size2;
  double s=0,u=0,v=0;
  double ni = (double)Xtrain->size1;
  double nj = (double)Xtrain->size2;

  /* get stats of training data */
  double *tmpj = (double *) VCalloc(Xtrain->size2,sizeof(double));
  for (j=0; j<Xtrain->size2; j++) {
    s=0;
    for (i=0; i<Xtrain->size1; i++) s += gsl_matrix_get(Xtrain,i,j);
    tmpj[j] = s/ni;
  }
  double sij = gsl_stats_mean(Xtrain->data,1,dim);

  /* apply to training data */
  for (i=0; i<Xtrain->size1; i++) {
    s=0;
    for (j=0; j<Xtrain->size2; j++) s += gsl_matrix_get(Xtrain,i,j);
    s /= nj;

    for (j=0; j<Xtrain->size2; j++) {
      u = gsl_matrix_get(Xtrain,i,j);
      v = u - s - tmpj[j] + sij;
      gsl_matrix_set(Xtrain,i,j,v);
    }
  }

  /* apply to test data */
  for (i=0; i<Xtest->size1; i++) {
    s=0;
    for (j=0; j<Xtest->size2; j++) s += gsl_matrix_get(Xtest,i,j);
    s /= nj;
    for (j=0; j<Xtest->size2; j++) {
      u = gsl_matrix_get(Xtest,i,j);
      v = u - s - tmpj[j] + sij;
      gsl_matrix_set(Xtest,i,j,v);
    }
  }
  VFree(tmpj);
}


void KXNorm(gsl_matrix *Xtrain,gsl_matrix *Xtest)
{
  size_t n = Xtrain->size1+Xtest->size1;
  size_t m = Xtrain->size2;
  size_t i,j;
  double u=0;

  gsl_matrix *X = gsl_matrix_calloc(n,m);

  for (j=0; j<m; j++) {
    for (i=0; i<Xtest->size1; i++) {
      u = gsl_matrix_get(Xtest,i,j);
      gsl_matrix_set(X,i,j,u);
    }
    for (i=0; i<Xtrain->size1; i++) {
      u = gsl_matrix_get(Xtrain,i,j);
      gsl_matrix_set(X,i+Xtest->size1,j,u);
    }
  }

  KNorm(X);
  for (j=0; j<m; j++) {
    for (i=0; i<Xtest->size1; i++) {
      u = gsl_matrix_get(X,i,j);
      gsl_matrix_set(Xtest,i,j,u);
    }
  }
  gsl_matrix_free(X);

  KNorm(Xtrain);
}



/* z-scoring across columns */
void MatNorm1(gsl_matrix *Xtrain,gsl_matrix *Xtest)
{
  size_t i,j;
  size_t n = Xtrain->size1;
  double u=0,z=0,mean=0,sd=1;

  gsl_vector *x = gsl_vector_calloc(n);
  
  /* for each column */
  for (j=0; j<Xtrain->size2; j++) {

    /* training */
    for (i=0; i<n; i++) x->data[i] = gsl_matrix_get(Xtrain,i,j);
    mean = gsl_stats_mean(x->data,1,n);
    sd = gsl_stats_sd_m(x->data,1,n,mean);
    if (sd < TINY) continue;
 
    for (i=0; i<n; i++) {
      u = gsl_matrix_get(Xtrain,i,j);
      z = (u-mean)/sd;
      gsl_matrix_set(Xtrain,i,j,z);
    }

    /* test */
    for (i=0; i<Xtest->size1; i++) {
      u = gsl_matrix_get(Xtest,i,j);
      z = (u-mean)/sd;
      gsl_matrix_set(Xtest,i,j,z);
    }
  }
  gsl_vector_free(x);
}


/* z-scoring */
void MatNorm0(gsl_matrix *X)
{
  size_t i,j;
  double u=0,z=0,mean=0,sd=1;

  gsl_vector *x = gsl_vector_calloc(X->size1);

  for (j=0; j<X->size2; j++) {
    
    for (i=0; i<X->size1; i++) x->data[i] = gsl_matrix_get(X,i,j);
    mean = gsl_stats_mean(x->data,1,x->size);
    sd = gsl_stats_sd_m(x->data,1,x->size,mean);
    if (sd < TINY) VError("MatNorm0, sd=0");

    for (i=0; i<X->size1; i++) {
      u = gsl_matrix_get(X,i,j);
      z = (u-mean)/sd;
      gsl_matrix_set(X,i,j,z);
    }
  }
  gsl_vector_free(x);
}


/* column centering */
void MatColumnCentering(gsl_matrix *X)
{
  size_t i,j;
  double u=0,mean=0;

  gsl_vector *x = gsl_vector_calloc(X->size1);

  for (j=0; j<X->size2; j++) {
    
    for (i=0; i<X->size1; i++) x->data[i] = gsl_matrix_get(X,i,j);
    mean = gsl_stats_mean(x->data,1,x->size);

    for (i=0; i<X->size1; i++) {
      u = gsl_matrix_get(X,i,j);
      gsl_matrix_set(X,i,j,u-mean);
    }
  }
  gsl_vector_free(x);
}


/* z-scoring across columns */
void MatNorm2(gsl_matrix *Xtrain,gsl_matrix *Xtest)
{
  size_t i,j;
  size_t n = Xtrain->size1 + Xtest->size1;
  double u=0,z=0,mean=0,sd=1;

  gsl_vector *x = gsl_vector_calloc(n);
  
  /* for each column */
  for (j=0; j<Xtrain->size2; j++) {

    /* training */
    for (i=0; i<Xtrain->size1; i++) x->data[i] = gsl_matrix_get(Xtrain,i,j);
    for (i=0; i<Xtest->size1; i++) x->data[i+Xtrain->size1] = gsl_matrix_get(Xtest,i,j);
    
    mean = gsl_stats_mean(x->data,1,n);
    sd = gsl_stats_sd_m(x->data,1,n,mean);
    if (sd < TINY) continue;
 
    for (i=0; i<n; i++) {
      u = gsl_matrix_get(Xtrain,i,j);
      z = (u-mean)/sd;
      gsl_matrix_set(Xtrain,i,j,z);
    }

    /* test */
    for (i=0; i<Xtest->size1; i++) {
      u = gsl_matrix_get(Xtest,i,j);
      z = (u-mean)/sd;
      gsl_matrix_set(Xtest,i,j,z);
    }
  }
  gsl_vector_free(x);
}



/* rank transform, gaussianize */
void GNorm(double *x,size_t n)
{
  size_t i;
  double u=0,z=0,nx=(double)n;

  size_t *perm = (size_t *) VCalloc(n,sizeof(size_t));
  gsl_sort_index (perm,x,1,n);

  double eps = 0.5/nx;
  for (i=0; i<n; i++) {
    u = ((double)(i))/nx;
    u += eps;
    if (u > 1.0-eps) u = 1.0-eps;
    if (u < eps) u = eps;
    z = gsl_cdf_ugaussian_Pinv(u);
    x[perm[i]] = z;
  }
  VFree(perm);
}



void MatRowNorm(gsl_matrix *X)
{
  size_t i,j;
  double *pp,mean=0,sd=0;
  for (i=0; i<X->size1; i++) {
    pp = gsl_matrix_ptr(X,i,0);
    mean = gsl_stats_mean(pp,1,X->size2);
    sd = gsl_stats_sd_m(pp,1,X->size2,mean);
    if (sd < TINY) continue;
    for (j=0; j<X->size2; j++) {
      pp[j] = (pp[j]-mean)/sd;
    }
  }
}



void RankMatch(gsl_vector *ytrain,gsl_vector *result)
{
  size_t i;
  size_t ntrain = ytrain->size;
  size_t ntest = result->size;
  double u=0,z=0,nx=(double)ntrain;

  RankNorm(result);
  
  gsl_vector *ztrain = gsl_vector_calloc(ntrain);
  gsl_vector_memcpy(ztrain,ytrain);
  gsl_sort(ztrain->data,1,ntrain);

  double *x = (double *) VCalloc(ntrain,sizeof(double));
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline,ntrain);
  for (i=0; i<ntrain; i++) x[i] = ((double)i)/nx;
  gsl_spline_init (spline,x,ztrain->data,ntrain);


  for (i=0; i<ntest; i++) {
    u = result->data[i];
    if (u < 0.001) u = 0.001;
    if (u > 0.999) u = 0.999;
    int rtcode = gsl_spline_eval_e(spline,u,acc,&z);
    if (rtcode==GSL_EDOM) {
      fprintf(stderr," err spline, u: %f\n",u);
      z = 0;
    }
    /* fprintf(stderr," %8.4f  %8.4f\n",result->data[i],z); */
    result->data[i] = z;
  }

  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  gsl_vector_free(ztrain);
  VFree(x);
}

