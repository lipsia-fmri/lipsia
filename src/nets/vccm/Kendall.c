/*
** Kendall's coefficient of concordance - Kendall's W
** 
** G.Lohmann, July 2010
*/
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

extern void gsl_sort_vector_index (gsl_permutation *,gsl_vector *);

#define SQR(x) ((x) * (x))
#define ABS(x) ((x) > 0 ? (x) : -(x))


double kendall(size_t **xrank,double *r,int n,int m,int otype)
{
  int i,j;
  double S,W,R,chi_square=0;
  double nx=0,mx=0;

  nx = (double)n;
  mx = (double)m;

  for (i=0; i<n; i++) {
    r[i] = 0;
    for (j=0; j<m; j++) r[i] += xrank[j][i];
  }

  R = 0;
  for (i=0; i<n; i++) R += r[i];
  R /= nx;

  S = 0;
  for (i=0; i<n; i++) S += SQR(r[i] - R);
  
  W = 12.0*S/(mx*mx*(nx*nx-1)*nx);
  if (otype == 0) return W;

  /* by central limit theorem, chi_square is approx normal for large n */
  /* see en.wikipedia.org/wiki/Chi-square_distribution */
  chi_square = (nx-1)*mx*W;
  return (chi_square - nx)/(sqrt(2.0*nx));
}

double
VKendall_W(gsl_matrix *data,int nvoxels,int otype)
{
  int i,j,n,m;
  double W0;
  static double *r=NULL;
  static gsl_vector *vec = NULL;
  static size_t **xrank=NULL;
  static gsl_permutation *perm=NULL,*rank=NULL;


  m = data->size1;
  n = nvoxels;

  if (vec == NULL) {
    r    = (double *) calloc(data->size2,sizeof(double));
    vec  = gsl_vector_calloc(data->size2);
    perm = gsl_permutation_alloc(data->size2);
    rank = gsl_permutation_alloc(data->size2);

    xrank = (size_t **) calloc(m,sizeof(size_t *));
    for (i=0; i<m; i++) xrank[i] = (size_t *) calloc(data->size2,sizeof(size_t));
  }

  gsl_vector_set_all (vec,9999);
  for (i=0; i<m; i++) {
    for (j=0; j<n; j++) {
      gsl_vector_set(vec,j,gsl_matrix_get(data,i,j));
    }
    gsl_sort_vector_index (perm, vec);
    gsl_permutation_inverse (rank, perm);
    for (j=0; j<n; j++) xrank[i][j] = rank->data[j];
  }

  /* Kendall's W of measured data */
  W0 = kendall(xrank,r,n,m,otype);
  return W0;
}
