/*
** use median of trials
** 
** G.Lohmann, May 2016
*/
#include <viaio/Vlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics.h>

#include <gsl/gsl_sort.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>

extern void VNormalize(float *data,int nt,VBoolean stddev);
extern void GetRank(float *data,gsl_vector *v,gsl_permutation *perm,gsl_permutation *rank);

void GetMedian(gsl_matrix_float **X1,gsl_matrix_float **X2,int *table,int n,gsl_matrix_float *SNR,int metric)
{
  long j,nvox=X1[0]->size1;
  long k,s,nt=X1[0]->size2;

  double *data = (double *) VCalloc(n,sizeof(double));
  double median = 0;

  gsl_vector *vec = NULL;
  gsl_permutation *perm = NULL;
  gsl_permutation *rank = NULL;
  if (metric == 1) {
    vec = gsl_vector_calloc(nt);
    perm = gsl_permutation_alloc(nt);
    rank = gsl_permutation_alloc(nt);
  }

  /* for each voxel ... */
  for (j=0; j<nvox; j++) {

    for (k=0; k<nt; k++) {
      for (s=0; s<n; s++) {
	if (table[s] == 0)
	  data[s] = (double)gsl_matrix_float_get(X1[s],j,k);
	else
	  data[s] = (double)gsl_matrix_float_get(X2[s],j,k);
      }
      gsl_sort(data,1,n);
      median = gsl_stats_median_from_sorted_data (data,1,(size_t) n);
      gsl_matrix_float_set(SNR,j,k,(float)median);
    }
    if (metric == 0) {  /* pearson correlation */
      float *qq = gsl_matrix_float_ptr(SNR,j,0);
      VNormalize(qq,nt,TRUE);
    }
    if (metric == 1) {  /* spearman ranks      */
      float *qq = gsl_matrix_float_ptr(SNR,j,0);
      GetRank(qq,vec,perm,rank);
    }
  }
  VFree(data);

  if (metric == 1) {
    gsl_permutation_free (perm);
    gsl_permutation_free (rank);
    gsl_vector_free(vec);
  }
}
