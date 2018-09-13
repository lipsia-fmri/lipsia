#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "viaio/Vlib.h"
#include "viaio/VImage.h"
#include "viaio/mu.h"
#include "viaio/option.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>


/* two-pass formula, correct for round-off error */
void VNormalize(gsl_vector *vec)
{
  double ave,var,sd,nx,s,u,tiny=1.0e-6;
  size_t j,nt = vec->size;

  nx = (double)nt;
  ave = 0;
  for (j=0; j<nt; j++) ave += vec->data[j];
  ave /= nx;

  var = u = 0;
  for (j=0; j<nt; j++) {
    s = vec->data[j]-ave;
    u   += s;
    var += s*s;
  }
  var=(var-u*u/nx)/(nx-1);
  sd = sqrt(var);
  if (sd < tiny) {
    for (j=0; j<nt; j++) vec->data[j] = 0;
  }
  else {
    for (j=0; j<nt; j++) {
      u = vec->data[j];
      vec->data[j] = (u-ave)/sd;
    }
  }
}


/* convert to ranks for spearman correlation */
void GetRank(gsl_vector *vec,gsl_permutation *perm,gsl_permutation *rank)
{
  size_t i;
  size_t n = vec->size;
  gsl_sort_vector_index (perm, vec);
  gsl_permutation_inverse (rank, perm);
  for (i=0; i<n; i++) vec->data[i] = (float)rank->data[i];
}


size_t NumVoxels(VImage src)
{
  int b,r,c;
  size_t n=0;
  for (b=0; b<VImageNBands(src); b++) {
    for (r=0; r<VImageNRows(src); r++) {
      for (c=0; c<VImageNColumns(src); c++) {
	if (VGetPixel(src,b,r,c) > 0.01) n++;
      }
    }
  }
  return n;
}

gsl_matrix_float *DataMatrix(VImage *src,int nslices,VImage roi,int metric)
{
  int b,r,c;
  float u=0;
  size_t nt = VImageNBands(src[0]);
  size_t nvox = NumVoxels(roi);
  gsl_matrix_float *X = gsl_matrix_float_calloc(nvox,nt);
  if (X==NULL) VError(" error allocating data matrix, size: %lu x %lu\n",nvox,nt);

  gsl_vector *vec = gsl_vector_calloc(nt);
  gsl_permutation *perm = NULL;
  gsl_permutation *rank = NULL;

  if (metric == 1) {  /* only needed for spearman correlation */
    perm = gsl_permutation_alloc(nt);
    rank = gsl_permutation_alloc(nt);
  }


  size_t i=0,j=0;
  for (b=0; b<VImageNBands(roi); b++) {
    for (r=0; r<VImageNRows(roi); r++) {
      for (c=0; c<VImageNColumns(roi); c++) {

	if (VGetPixel(roi,b,r,c) < 0.01) continue;

	for (j=0; j<nt; j++) {
	  u = (float) VGetPixel(src[b],j,r,c);
	  vec->data[j] = u;
	}

	if (metric == 0) {  /* linear correlation */
	  VNormalize(vec);
	}
	else if (metric == 1) {  /* convert to ranks for Spearman corr */
	  GetRank(vec,perm,rank);
	}
	else {
	  VError(" unknown metric %d",metric);
	}

	for (j=0; j<nt; j++) {
	  gsl_matrix_float_set(X,i,j,(float)vec->data[j]);
	}

	i++;
      }
    }
  }

  gsl_vector_free(vec);
  if (metric==1) {
    gsl_permutation_free (perm);
    gsl_permutation_free (rank);
  }
  return X;
}


/* new voxel map */
VImage VoxelMap(VImage roi)
{
  size_t nvox = NumVoxels(roi);
  VImage map = VCreateImage(1,4,nvox,VShortRepn);
  if (map == NULL) VError(" error allocating addr map");
  VFillImage(map,VAllBands,0); 
  int nslices = VImageNBands(roi);
  int nrows = VImageNRows(roi);
  int ncols = VImageNColumns(roi);
  VSetAttr(VImageAttrList(map),"nvoxels",NULL,VLongRepn,(VLong)nvox);
  VSetAttr(VImageAttrList(map),"nslices",NULL,VLongRepn,(VLong)nslices);
  VSetAttr(VImageAttrList(map),"nrows",NULL,VLongRepn,(VLong)nrows);
  VSetAttr(VImageAttrList(map),"ncols",NULL,VLongRepn,(VLong)ncols);
  VPixel(map,0,3,0,VShort) = nslices;
  VPixel(map,0,3,1,VShort) = nrows;
  VPixel(map,0,3,2,VShort) = ncols;

  int b,r,c;
  size_t i = 0;
  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {

	float u = VGetPixel(roi,b,r,c);
	if (u < 0.1) continue;
	VPixel(map,0,0,i,VShort) = b;
	VPixel(map,0,1,i,VShort) = r;
	VPixel(map,0,2,i,VShort) = c;
	i++;
      }
    }
  }
  return map;
}
