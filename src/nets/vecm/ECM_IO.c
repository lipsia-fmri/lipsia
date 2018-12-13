/*
** ECM - eigenvector centrality mapping, I/O interface
**
** G.Lohmann, MPI-KYB, updated Oct 2018
*/

#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/

#define SQR(x) ((x) * (x))
#define ReLU(x) ((x) > 0 ? (x) : 0)


VImage WriteOutput(VAttrList list,VImage map,float *eigvec,size_t n)
{
  size_t i;
  int b,r,c;
  int nrows=0,ncols=0,ntimesteps=0;
  int nslices = VAttrListNumImages(list);
  VImage *src = VAttrListGetImages(list,nslices);
  VImageDimensions(src,nslices,&ntimesteps,&nrows,&ncols);
  
  VImage dest = VCreateImage(nslices,nrows,ncols,VFloatRepn);
  VFillImage(dest,VAllBands,0);
  VCopyImageAttrs (src[0], dest);
  VSetAttr(VImageAttrList(dest),"modality",NULL,VStringRepn,"ECM");

  for (i=0; i<n; i++) {
    b = VPixel(map,0,0,i,VInteger);
    r = VPixel(map,0,1,i,VInteger);
    c = VPixel(map,0,2,i,VInteger);
    if (b < 0 || r < 0 || c < 0) VError("  WriteOutput: illegal address %d %d %d",b,r,c);
    VPixel(dest,b,r,c,VFloat) = eigvec[i];
  }
  return dest;
}



void NormalizeData(gsl_matrix_float *X,int nt)
{
  int i,j;
  double s1,s2,mean,sd,u,nx,tiny=1.0e-6;
  float *arr=NULL;
  nx = (double)nt;

  for (i=0; i<X->size1; i++) {
    s1=s2=0;
    arr = gsl_matrix_float_ptr(X,i,0);
    for (j=0; j<nt; j++) {
      u = (double)arr[j];
      s1 += u;
      s2 += u*u;
    }
    mean = s1/nx;
    sd = sqrt((s2 - nx * mean * mean) / (nx - 1.0));

    arr = gsl_matrix_float_ptr(X,i,0);
    for (j=0; j<nt; j++) {
      u = (double)arr[j];
      if (sd > tiny) arr[j] = (float)((u-mean)/sd);
      else arr[j] = 0;
    }
  }

  for (i=0; i<X->size1; i++) {
    for (j=0; j<X->size2; j++) {
      double u = (double)gsl_matrix_float_get(X,i,j);
      if (gsl_isnan(u) || gsl_isinf(u)) gsl_matrix_float_set(X,i,j,0);
    }
  }
}



VImage VoxelMap(VAttrList list,VImage mask,size_t nvox)
{
  int b,r,c;
  
  /* read image dims */
  int nrows=0,ncols=0,ntimesteps=0;
  int nslices = VAttrListNumImages(list);
  VImage *src = VAttrListGetImages(list,nslices);
  VImageDimensions(src,nslices,&ntimesteps,&nrows,&ncols);
  
 
  /*  voxel addresses */
  VImage map = VCreateImage(1,5,nvox,VIntegerRepn);
  if (map == NULL) VError(" error allocating addr map");
  VFillImage(map,VAllBands,0);
  VPixel(map,0,3,0,VInteger) = nslices;
  VPixel(map,0,3,1,VInteger) = nrows;
  VPixel(map,0,3,2,VInteger) = ncols;

  int i = 0;
  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	if (VGetPixel(mask,b,r,c) < 0.1) continue;
	VPixel(map,0,0,i,VInteger) = b;
	VPixel(map,0,1,i,VInteger) = r;
	VPixel(map,0,2,i,VInteger) = c;
	i++;
      }
    }
  }
  return map;
}

/* count number of voxels */
size_t NumVoxels(VImage mask)
{
  int b,r,c;
  size_t nvox = 0;
  for (b=0; b<VImageNBands(mask); b++) {
    for (r=0; r<VImageNRows(mask); r++) {
      for (c=0; c<VImageNColumns(mask); c++) {
	if (VGetPixel(mask,b,r,c) < 0.1) continue;
	nvox++;
      }
    }
  }
  fprintf(stderr," nvoxels: %ld\n",(long)nvox);
  return nvox;
}


gsl_matrix_float *ReadDataECM(VAttrList list,VImage mask,VImage map,VShort first,VShort length,int type,size_t nvox)
{
  size_t i,j;
  int b,r,c;

  /* read image dims */
  int nrows=0,ncols=0,ntimesteps=0,nt=0;
  int nslices = VAttrListNumImages(list);
  VImage *src = VAttrListGetImages(list,nslices);
  VImageDimensions(src,nslices,&ntimesteps,&nrows,&ncols);

  if (VImageNRows(mask) != nrows) 
    VError(" inconsistent image dims, mask has %d rows, image has %d rows",VImageNRows(mask),nrows);
  if (VImageNColumns(mask) != ncols) 
    VError(" inconsistent image dims, mask has %d columns, image has %d columns",VImageNColumns(mask),ncols);
  if (VImageNBands(mask) != nslices) 
    VError(" inconsistent image dims, mask has %d slices, image has %d slices",VImageNBands(mask),nslices);

  
   /* get time steps to include */
  if (length < 1) length = ntimesteps-1;
  int last = first + length;
  if (last >= ntimesteps) last = ntimesteps-1;
  if (first < 0) first = 0;
  nt = last - first + 1;
  if (nt < 2) VError(" not enough timesteps, nt= %d",nt);
  fprintf(stderr," ntimesteps: %d, first= %d, last= %d, nt= %d\n",
	  (int)ntimesteps,(int)first,(int)last,(int)nt);

  
  /* alloc and fill  matrix X */
  size_t dim = nt;
  if (type == 0 || type == 7) dim = 2*nt;
  gsl_matrix_float *X = gsl_matrix_float_calloc(nvox,dim);
  for (i=0; i<nvox; i++) {

    b = VPixel(map,0,0,i,VInteger);
    r = VPixel(map,0,1,i,VInteger);
    c = VPixel(map,0,2,i,VInteger);

    float *ptr = gsl_matrix_float_ptr(X,i,0);
    int k;
    j = 0;
    for (k=first; k<=last; k++) {
      if (j >= X->size2) VError(" j= %d %d",j,X->size2);
      if (k >= VImageNBands(src[b])) VError(" k= %d %d",k, VImageNBands(src[b]));
      *ptr++ = (float) VGetPixel(src[b],k,r,c);
      j++;
    }
  }


  /* normalize */
  NormalizeData(X,nt);

  return X;
}
