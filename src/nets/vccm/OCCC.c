/*
** overall concordance correlation coefficient (OCCC)
** Barnhart et al. (2002), Biometrics, 58, 1020-1027.
** 
** G.Lohmann, July 2010
*/
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

#define SQR(x) ((x) * (x))
#define ABS(x) ((x) > 0 ? (x) : -(x))


double
Cov(double *arr1,double *arr2,int n)
{
  int i;
  double nx,sum,ave1,ave2;

  nx = (double)n;

  sum = 0;
  for (i=0; i<n; i++) sum += arr1[i];
  ave1 = sum/nx;
    
  sum = 0;
  for (i=0; i<n; i++) sum += arr2[i];
  ave2 = sum/nx;
    
  sum = 0;
  for (i=0; i<n; i++)
    sum += (arr1[i] - ave1)*(arr2[i] - ave2);
  
  return sum/(nx-1);
}


void
AveVar(double *arr,int n,double *ave,double *var)
{
  int i;
  double u,sum1,sum2,nx,mean;

  nx = (double)n;
  sum1 = sum2 = 0;
  for (i=0; i<n; i++) {
    u = arr[i];
    sum1 += u;
    sum2 += u*u;
  }
  mean = sum1/nx;

  *ave = mean;
  *var = (sum2 - nx * mean * mean) / (nx - 1.0);
}


double
VOCCC(gsl_matrix *data,int nvoxels,int verbose)
{
  int i,j,k;
  double sum1,sum2,sum3,mean,varx,u,v;
  static gsl_matrix *cov=NULL;
  static gsl_vector *ave=NULL,*var=NULL;

  int nsubjects = data->size1;

  if (ave == NULL) {
    ave = gsl_vector_calloc(nsubjects);
    var = gsl_vector_calloc(nsubjects);
    cov = gsl_matrix_calloc(nsubjects,nsubjects);
  }

  for (i=0; i<nsubjects; i++) {
    double *ptr1 = gsl_matrix_ptr(data,i,0);
    AveVar(ptr1,nvoxels,&mean,&varx);
    gsl_vector_set(ave,i,mean);
    gsl_vector_set(var,i,varx);

    for (j=0; j<nsubjects; j++) {
      double *ptr1 = gsl_matrix_ptr(data,i,0);
      double *ptr2 = gsl_matrix_ptr(data,j,0);
      u = Cov(ptr1,ptr2,nvoxels);
      gsl_matrix_set(cov,i,j,u);
    }
  }
  if (verbose > 0) {
    for (i=0; i<nsubjects; i++) {
      for (j=0; j<nsubjects; j++) {
	fprintf(stderr," %8.5f",gsl_matrix_get(cov,i,j));
      }
      fprintf(stderr,"\n");
    }
    fprintf(stderr,"\n");
    exit(0);
  }
  

  sum1 = 0;
  for (j=1; j<=nsubjects-1; j++) {
    for (k=j+1; k<=nsubjects; k++) {
      u = gsl_matrix_get(cov,j-1,k-1);
      sum1 += u;
    }
  }
  sum1 *= 2.0;


  sum2 = 0;
  for (j=0; j<nsubjects; j++) {
    sum2 += gsl_vector_get(var,j);
  }
  sum2 *= (double)(nsubjects-1);

  
  sum3 = 0;
  for (j=1; j<=nsubjects-1; j++) {
    for (k=j+1; k<=nsubjects; k++) {
      u = gsl_vector_get(ave,j-1);
      v = gsl_vector_get(ave,k-1);
      sum3 += (u-v)*(u-v);
    }
  }

  if (sum2 + sum3 < 1.0e-6) return 0;
  return sum1 / (sum2 + sum3);
}
