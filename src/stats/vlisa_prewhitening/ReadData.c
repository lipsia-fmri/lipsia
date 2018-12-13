/*
** Read input data into array
**
** G.Lohmann, MPI-KYB, 2018
*/
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>

#define TINY 1.0e-8


/*
** get global mean regressor
*/
void GlobalMean(gsl_matrix *Data,gsl_matrix *covariates,int column)
{
  int i=0,j=0;
  int nt = Data->size2;

  double sum=0,nx=(double)Data->size1;
  for (i=0; i<nt; i++) {
    for (j=0; j<Data->size1; j++) {
      sum += gsl_matrix_get(Data,j,i);
    }
    gsl_matrix_set(covariates,i,column,sum/nx);
  }

  /* normalize */
  double u=0,s1=0,s2=0,mx=(double)Data->size2;
  for (i=0; i<nt; i++) {
    u = gsl_matrix_get(covariates,i,column);
    s1 += u;
    s2 += u*u;
  }
  double mean = s1/mx;
  double var = (s2 - mx * mean * mean) / (mx - 1.0);
  double sd = sqrt(var);

  for (i=0; i<nt; i++) {
    u = gsl_matrix_get(covariates,i,column);
    gsl_matrix_set(covariates,i,column,(u-mean/sd));
  }
}



/* normalize time courses */
void VRowNormalize(gsl_matrix *Data)
{
  size_t i,j;
  double u=0,s1=0,s2=0,mean=0,sd=0;
  double nx=(double)Data->size2;

  for (i=0; i<Data->size1; i++) {
    s1 = s2 = 0;
    for (j=0; j<Data->size2; j++) {
      u = gsl_matrix_get(Data,i,j);
      s1 += u;
      s2 += u*u;
    }
    mean = s1/nx;
    sd = sqrt((s2 - nx * mean * mean) / (nx - 1.0));
    
    for (j=0; j<Data->size2; j++) {
      u = gsl_matrix_get(Data,i,j);
      gsl_matrix_set(Data,i,j,(u-mean)/sd);
    }
  }
}


/* get timing information */
void VGetTimeInfos(VAttrList *list,int nlists,double *mtr,float *run_duration)
{
  VImage src=NULL;
  VAttrListPosn posn;
  int k,ntimesteps=0;
  VLong itr=0;
  double tr=-1.0;

  int ntt=0;
  for (k=0; k<nlists; k++) {
    ntt = 0;
    for (VFirstAttr (list[k], & posn); VAttrExists (& posn); VNextAttr (& posn)) {
      if (VGetAttrRepn (& posn) != VImageRepn) continue;
      VGetAttrValue (& posn, NULL,VImageRepn, & src);
      if (tr < 0) {
	if (VGetAttr(VImageAttrList(src),"repetition_time",NULL,VLongRepn,&itr) != VAttrFound) 
	  VError(" TR info missing in header");
	tr = (double) itr / 1000.0;
      }
      ntt = VImageNBands(src);
      run_duration[k] = tr*(float)ntt;
    }
    ntimesteps += ntt;
  }
  if (tr < 0) VError(" TR info missing");
  *mtr = tr;
}



/* read data block */
gsl_matrix *VReadImageData(VAttrList *list,int nlists)
{
  VAttrListPosn posn;
  int i=0,k,ntt=0,slice,row,col,nrows=0,ncols=0,ntimesteps=0;
  double u=0,s1=0,s2=0,nx=0,sd=0,mean=0,tiny=1.0e-8;

  int nslices = VAttrListNumImages(list[0]);

  /* alloc src image */
  VImage **src = (VImage **) VCalloc(nlists,sizeof(VImage *));
  int *nt = (int *) VCalloc(nlists,sizeof(int));
  for (k=0; k<nlists; k++) {
    src[k] = (VImage *) VCalloc(nslices,sizeof(VImage));
    i = 0;
    for (VFirstAttr (list[k], & posn); VAttrExists (& posn); VNextAttr (& posn)) {
      if (VGetAttrRepn (& posn) != VImageRepn) continue;
      if (i >= nslices) VError(" inconsistent number of slices, %d %d %d",k,i,nslices);
      VGetAttrValue (& posn, NULL,VImageRepn, & src[k][i]);
      i++;
    }
  }

  /* get image dimensions */
  nrows = VImageNRows(src[0][0]);
  ncols = VImageNColumns(src[0][0]);
  ntimesteps = 0;
  for (k=0; k<nlists; k++) {
    ntt = VImageNBands(src[k][0]);
    ntimesteps += ntt;
    nt[k] = ntt;
    if (nrows != VImageNRows(src[k][0])) VError(" inconsistent number of rows ");
    if (ncols != VImageNColumns(src[k][0])) VError(" inconsistent number of columns ");
  }

  /* number of nonzero voxels */
  size_t nvox=0;
  for (slice=0; slice<nslices; slice++) {
    for (row=0; row<nrows; row++) {
      for (col=0; col<ncols; col++) {
	u = VGetPixel(src[0][slice],0,row,col);
	if (gsl_isnan(u) || gsl_isinf(u)) continue;
	if (fabs(u) < TINY) continue;
	nvox++;
      }
    }
  }

  /* alloc new data struct */
  gsl_matrix *Data = gsl_matrix_calloc(nvox,ntimesteps);
  if (!Data) VError(" error allocating Data");

  
  /* fill map and data */
  for (k=0; k<nlists; k++) {

    ntt = 0;
    for (i=0; i<k; i++) ntt += nt[i];	  

    nvox=0;
    for (slice=0; slice<nslices; slice++) {
      for (row=0; row<nrows; row++) {
	for (col=0; col<ncols; col++) {
	  u = VGetPixel(src[0][slice],0,row,col);
	  if (gsl_isnan(u) || gsl_isinf(u)) continue;
	  if (fabs(u) < TINY) continue;

	  s1=s2=nx=0;
	  for (i=0; i<nt[k]; i++) {
	    u = VGetPixel(src[k][slice],i,row,col);
	    s1 += u;
	    s2 += u*u;	    
	    nx++;
	  }
	  mean = s1/nx;
	  sd = sqrt((s2 - nx * mean * mean) / (nx - 1.0));
	  if (sd < tiny) {
	    u = 0.0;
	    for (i=0; i<nt[k]; i++) {
	      gsl_matrix_set(Data,nvox,ntt+i,u);
	    }
	  }
	  else {
	    for (i=0; i<nt[k]; i++) {
	      u = VGetPixel(src[k][slice],i,row,col);
	      u = (u-mean)/sd;
	      gsl_matrix_set(Data,nvox,ntt+i,u);
	    }
	  }
	  nvox++;
	}
      }
    }
  }
  VFree(nt);
  return Data;
}


/* get map of voxel addresses */
VImage VoxelMap(VAttrList list)
{ 
  VAttrListPosn posn;
  int i=0,slice,row,col,nrows=0,ncols=0;

  int nslices = VAttrListNumImages(list);

  /* alloc src image */
  i = nrows = ncols = 0;
  VImage *src = (VImage *) VCalloc(nslices,sizeof(VImage));
  i = 0;
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    if (i >= nslices) VError(" inconsistent number of slices, %d %d",i,nslices);
    VGetAttrValue (& posn, NULL,VImageRepn, & src[i]);
    if (VImageNRows(src[i])  > nrows)  nrows = VImageNRows(src[i]);
    if (VImageNColumns(src[i]) > ncols) ncols = VImageNColumns(src[i]);
    i++;
  }


  /* number of nonzero voxels */
  double u=0;
  size_t nvox=0;
  for (slice=0; slice<nslices; slice++) {
    for (row=0; row<nrows; row++) {
      for (col=0; col<ncols; col++) {
	u = VGetPixel(src[slice],0,row,col);
	if (fabs(u) < TINY) continue;
	nvox++;
      }
    }
  }
  if (nvox < 1) VError(" No non-zero voxels in input image");


  /* voxel addresses */
  VImage map = VCreateImage(1,4,nvox,VShortRepn);
  if (map == NULL) VError(" error allocating addr map");
  VFillImage(map,VAllBands,0);
  VCopyImageAttrs (src[0],map);
  VSetAttr(VImageAttrList(map),"nvoxels",NULL,VLongRepn,(VLong)nvox);
  VSetAttr(VImageAttrList(map),"nslices",NULL,VLongRepn,(VLong)nslices);
  VSetAttr(VImageAttrList(map),"nrows",NULL,VLongRepn,(VLong)nrows);
  VSetAttr(VImageAttrList(map),"ncols",NULL,VLongRepn,(VLong)ncols);
  VPixel(map,0,3,0,VShort) = nslices;
  VPixel(map,0,3,1,VShort) = nrows;
  VPixel(map,0,3,2,VShort) = ncols;


  /* fill map */
  nvox=0;
  for (slice=0; slice<nslices; slice++) {
    for (row=0; row<nrows; row++) {
      for (col=0; col<ncols; col++) {
	u = VGetPixel(src[slice],0,row,col);
	if (fabs(u) < TINY) continue;
	VPixel(map,0,0,nvox,VShort) = slice;
	VPixel(map,0,1,nvox,VShort) = row;
	VPixel(map,0,2,nvox,VShort) = col;
	nvox++;
      }
    }
  }
  return map;
}
