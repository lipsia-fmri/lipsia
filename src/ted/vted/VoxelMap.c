/*
** read data matrix, map of voxel addresses
**
** G.Lohmann, Jan 2013
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort.h>

#include "viaio/Vlib.h"
#include "viaio/VImage.h"
#include "viaio/mu.h"

#define SQR(x) ((x) * (x))
#define ABS(x) ((x) > 0 ? (x) : -(x))


VImage *VImagePointer(VAttrList list,int *nt)
{
  VImage tmp;
  VAttrListPosn posn;
  int i,ntimesteps,nrows,ncols,nslices;
 
  /* get image dimensions, read functional data */
  nslices = 0;
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    VGetAttrValue (& posn, NULL,VImageRepn, & tmp);
    if (VPixelRepn(tmp) != VShortRepn) continue;
    nslices++;
  }
  /* VDestroyImage(tmp); */
  if (nslices < 1) VError(" no slices");

  VImage *src = (VImage *) VCalloc(nslices,sizeof(VImage));
  i = ntimesteps = nrows = ncols = 0;
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    VGetAttrValue (& posn, NULL,VImageRepn, & src[i]);
    if (VPixelRepn(src[i]) != VShortRepn) continue;
    if (VImageNBands(src[i]) > ntimesteps) ntimesteps = VImageNBands(src[i]);
    if (VImageNRows(src[i])  > nrows)  nrows = VImageNRows(src[i]);
    if (VImageNColumns(src[i]) > ncols) ncols = VImageNColumns(src[i]);
    i++;
  }
  nslices = i;
  if (nslices < 1) return NULL;
  if (ntimesteps < 2) return NULL;
  *nt = ntimesteps;
  return src;
}


/* read functional data */
void VReadImagePointer(VAttrList list,VImage *src)
{
  VAttrListPosn posn;
 
  int i = 0;
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    VGetAttrValue (& posn, NULL,VImageRepn, & src[i]);
    if (VPixelRepn(src[i]) != VShortRepn) continue;
    i++;
  }
}

/* check if mask covers data. If not, set surplus voxels to zero */
long VMaskCoverage(VAttrList list,VImage mask)
{
  VImage src=NULL;
  VAttrListPosn posn;
  int b,r,c;
  long count=0;

  b = 0; 
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    VGetAttrValue (& posn, NULL,VImageRepn, & src);
    if (VPixelRepn(src) != VShortRepn) continue;
    if (VImageNRows(src) < 3) continue;

    if (VImageNRows(src) != VImageNRows(mask)) 
      VError(" inconsistent nrows: %d,  mask: %d",VImageNRows(src),VImageNRows(mask));
    if (VImageNColumns(src) != VImageNColumns(mask)) 
      VError(" inconsistent ncols: %d,  mask: %d",VImageNColumns(src),VImageNColumns(mask));

    for (r=0; r<VImageNRows(src); r++) {
      for (c=0; c<VImageNColumns(src); c++) {
	if (b >= VImageNBands(mask)) VError("VMaskCoverage, illegal addr, band= %d",b);
	float u = VGetPixel(mask,b,r,c);
	int j = (int)VPixel(src,0,r,c,VShort);
	if (j == 0 && u > 0.3) {
	  VSetPixel(mask,b,r,c,0);
	  count++;
	}
      }
    }
    b++;
  }
  return count;
}


/* check if data matrix contains zero voxels */
void VCheckMatrix(gsl_matrix_float *X)
{
  long i,j,nvox=X->size1,nt=X->size2;
  double sum1,sum2,nx,mean,var,tiny=1.0e-4;

  long count = 0;
  nx = (double)nt;
  for (i=0; i<nvox; i++) {
    const float *arr1 = gsl_matrix_float_const_ptr(X,i,0);
    sum1 = sum2 = 0;
    for (j=0; j<nt; j++) {
      const double u = (double)arr1[j];
      sum1 += u;
      sum2 += u*u;
    }
    mean = sum1/nx;
    var = (sum2 - nx * mean * mean) / (nx - 1.0);
    if (var < tiny) count++;
  }
  if (count > 0)
    VWarning(" number of empty voxels: %ld",count);
}


VImage VoxelMap(VImage mask,size_t *nvoxels)
{
  int b,r,c;
  int nslices = VImageNBands(mask);
  int nrows = VImageNRows(mask);
  int ncols = VImageNColumns(mask);

  /* count number of non-zero voxels */
  size_t nvox = 0;  
  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	float u = VGetPixel(mask,b,r,c);
	if (ABS(u) < 0.5) continue;
	nvox++;
      }
    }
  }
  (*nvoxels) = nvox;

  /* voxel addresses */
  VImage map = VCreateImage(1,4,nvox,VShortRepn);
  if (map == NULL) VError(" error allocating addr map");
  VFillImage(map,VAllBands,0);
  VImage xmap=NULL;
  VExtractAttr (VImageAttrList(mask),"map",NULL,VImageRepn,&xmap,FALSE);
  VCopyImageAttrs (mask,map);
  VSetAttr(VImageAttrList(map),"nvoxels",NULL,VLongRepn,(VLong)nvox);
  VSetAttr(VImageAttrList(map),"nslices",NULL,VLongRepn,(VLong)nslices);
  VSetAttr(VImageAttrList(map),"nrows",NULL,VLongRepn,(VLong)nrows);
  VSetAttr(VImageAttrList(map),"ncols",NULL,VLongRepn,(VLong)ncols);

  VPixel(map,0,3,0,VShort) = nslices;
  VPixel(map,0,3,1,VShort) = nrows;
  VPixel(map,0,3,2,VShort) = ncols;

  size_t i = 0;
  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	float u = VGetPixel(mask,b,r,c);
	if (ABS(u) < 0.5) continue;
	VPixel(map,0,0,i,VShort) = b;
	VPixel(map,0,1,i,VShort) = r;
	VPixel(map,0,2,i,VShort) = c;
	i++;
      }
    }
  }
  return map;
}

/* two-pass formula, correct for round-off error */
void VNormalize(float *data,int nt,VBoolean stddev)
{
  int j;
  float ave,var,sd,nx,s,u,tiny=1.0e-6;

  nx = (float)nt;
  ave = 0;
  for (j=0; j<nt; j++) ave += data[j];
  ave /= nx;

  var = u = 0;
  for (j=0; j<nt; j++) {
    s = data[j]-ave;
    u   += s;
    var += s*s;
  }
  var=(var-u*u/nx)/(nx-1);
  sd = sqrt(var);
  if (sd < tiny) {
    for (j=0; j<nt; j++) data[j] = 0;
  }
  else {
    if (stddev==FALSE) sd = 1.0;  /* only subtract mean !!!! */
    for (j=0; j<nt; j++) {
      u = data[j];
      data[j] = (u-ave)/sd;
    }
  }
}


/*  copy data to input matrix X, contains fMRI time series */
void VDataMatrix(VImage *src,int first,int len,VImage map,gsl_matrix_float *X)
{
  long i,j,k,nvox;
  int b,r,c;
  float *data=NULL;

  data = (float *) VCalloc(len,sizeof(float));

  gsl_matrix_float_set_zero (X);
  nvox = VImageNColumns(map);
  if (nvox != (long)X->size1) VError(" err, %ld %ld",(long)X->size1,nvox);
  if (len  != (long)X->size2) VError(" err, %ld %ld",(long)X->size2,len);

  for (i=0; i<nvox; i++) {
    
    b = VPixel(map,0,0,i,VShort);
    r = VPixel(map,0,1,i,VShort);
    c = VPixel(map,0,2,i,VShort);

    if ((first+len) > VImageNBands(src[b])) VError(" illegal len addr, %d %d",first+len,VImageNBands(src[b]));
    if (r >= VImageNRows(src[b])) VError(" illegal row addr");
    if (c >= VImageNColumns(src[b])) VError(" illegal column addr");

    k = 0;
    for (j=first; j<first+len; j++) {
      data[k] = (float) VGetPixel(src[b],j,r,c);
      k++;
    }
    VNormalize(data,len,TRUE);

    float *ptr = gsl_matrix_float_ptr(X,i,0);
    for (j=0; j<len; j++) *ptr++ = data[j];
  }
  VFree(data);
}

