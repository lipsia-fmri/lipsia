#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/* read voxel value */
double VGetVoxel(VImage *src,int nimages,int b,int r,int c,int j)
{
  double u=0;
  if (nimages <= 1) {  /* 3D */
    u = VGetPixel(src[j],b,r,c);
  }
  else {  /* 4D */
    u = VGetPixel(src[b],j,r,c);
  }
  return u;
}

/* write a voxel value */
int VSetVoxel(VImage *src,int nvolumes,int b,int r,int c,int j,double value)
{
  if (j >= nvolumes) return -1;
  if (nvolumes <= 1) {  /* 3D */
    VSetPixel(src[j],b,r,c,value);
  }
  else {  /* 4D */
    VSetPixel(src[b],j,r,c,value);
  }
  return 1;
}

/* read a single image from a file */
VImage VReadImageFile(VString filename)
{
  VImage image=NULL;
  VAttrList list=NULL;
  if (strlen(filename) > 1) {
    list = VReadAttrList(filename,0L,TRUE,FALSE);
    image = VReadImage(list);
  }
  return image;
}


/* apply mask or threshold */
void VMaskMinval(VAttrList list,VImage mask,double minval)
{
  int b,r,c,j;
  double u;
  if (mask == NULL && minval < NO_MINVAL+0.001) return;

  int nslices=0,nrows=0,ncols=0,nvolumes=0;
  int nimages = VAttrListNumImages(list);
  VImage *src = VAttrListGetImages(list,nimages);
  VDimensions(src,nimages,&nslices,&nrows,&ncols,&nvolumes);

  if (mask != NULL) {
    if (nslices != VImageNBands(mask))
      VError(" Number of slices in mask inconsistent with input image (%d %d)",(int)VImageNBands(mask),nslices);
    if (nrows != VImageNRows(mask))
      VError(" Number of rows in mask inconsistent with input image (%d %d)",(int)VImageNRows(mask),nrows);
    if (ncols != VImageNColumns(mask))
      VError(" Number of columns in mask inconsistent with input image (%d %d)",(int)VImageNColumns(mask),ncols);
  }
  
  for (j=0; j<nvolumes; j++) {
    for (b=0; b<nslices; b++) {
      for (r=0; r<nrows; r++) {
	for (c=0; c<ncols; c++) {
	  
	  if (mask != NULL) {
	    u = VGetPixel(mask,b,r,c);
	    if (fabs(u) < TINY) {
	      VSetVoxel(src,nvolumes,b,r,c,j,0.0);
	    }
	  }
	  else {
	    u = VGetVoxel(src,nvolumes,b,r,c,j);
	    if (u < minval) {
	      VSetVoxel(src,nvolumes,b,r,c,j,0.0);
	    }
	  }
	}
      }
    }
  }
}

/* apply mask or threshold to a single list, read mask file first */
void VMinval(VAttrList list,VString mask_filename,double minval)
{
  VImage mask = VReadImageFile(mask_filename);
  if (mask == NULL && minval < NO_MINVAL+0.001) return;
  VMaskMinval(list,mask,minval);
}


/* apply mask or threshold to multiple lists */
void VMultMinval(VAttrList *list,int nlists,VString mask_filename,double minval)
{
  int k;
  VImage mask = VReadImageFile(mask_filename);
  if (mask == NULL && minval < NO_MINVAL+0.001) return;
  for (k=0; k<nlists; k++) {
    VMaskMinval(list[k],mask,minval);
  }
}


/* 
** check if time course is a foreground voxel,
** normalize if needed
*/
int Foreground(double *data,size_t n,int norm)
{
  size_t i=0;
  double u=0,nx=0,sum1=0,sum2=0,mean=0,var=0,sd=0;

  for (i=0; i<n; i++) {
    u = data[i];
    sum1 += u;
    sum2 += u*u;
  }
  nx = (double)n;
  mean = sum1/nx;
  var = (sum2 - nx * mean * mean) / (nx - 1.0);
  if (var < TINY) return -1;

  /* subtract mean */
  if (norm == 1) {
    for (i=0; i<n; i++) {
      u = data[i];
      data[i] = (u-mean);
    }
  }

  /* z-scoring */
  if (norm == 2) {
    sd = sqrt(var);
    for (i=0; i<n; i++) {
      u = data[i];
      data[i] = (u-mean)/sd;
    }
  }
  return 1;
}


int VCheckMinval(VImage *src,int b,int r,int c,int nt)
{
  int j;
  double s=0;
  for (j=0; j<nt; j++) {
    s += VGetPixel(src[b],j,r,c);
  }
  if (fabs(s) < TINY) return -1;
  return 1;
}
