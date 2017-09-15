#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/* estimate threshold for image background, multiple lists */
float VGetMinvalNlists(VAttrList *list,int nlists)
{
  VImage src=NULL;
  VAttrListPosn posn;
  int row,col,i,k;
  double umin = 1.0e+12;
  double sum=0,nx=0,mx=0,u=0,tiny=1.0e-8;

  for (k=0; k<nlists; k++) {
    for (VFirstAttr (list[k], & posn); VAttrExists (& posn); VNextAttr (& posn)) {
      if (VGetAttrRepn (& posn) != VImageRepn) continue;
      VGetAttrValue (& posn, NULL,VImageRepn, & src);
      for (i=0; i<VImageNBands(src); i++) {
	for (row=0; row<VImageNRows(src); row++) {
	  for (col=0; col<VImageNColumns(src); col++) {
	    u = VGetPixel(src,i,row,col);
	    if (fabs(u) > tiny) {
	      if (u < umin) umin = u;
	      sum += u;
	      nx++;
	    }
	    mx++;
	  }
	}
      }
    }
  }
  double mean = sum/nx;
  float minval = (float)(mean*0.875);

  /* many zero voxels, probably pre-segmented */
  if (nx/mx < 0.875) minval = (float)(umin-tiny);

  return minval;
}


/* estimate threshold for image background, single list */
float VGetMinval(VAttrList list)
{
  VImage src=NULL;
  VAttrListPosn posn;
  int row,col,i;
  double sum=0,nx=0,mx=0,u=0,tiny=1.0e-8;
  double umin = 1.0e+12;

  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    VGetAttrValue (& posn, NULL,VImageRepn, & src);
    for (i=0; i<VImageNBands(src); i++) {
      for (row=0; row<VImageNRows(src); row++) {
	for (col=0; col<VImageNColumns(src); col++) {
	  u = VGetPixel(src,i,row,col);
	  if (fabs(u) > tiny) {
	    if (u < umin) umin = u;
	    sum += u;
	    nx++;
	  }
	  mx++;
	}
      }
    }
  }
  double mean = sum/nx;
  float minval = (float)(mean*0.875);

  /* many zero voxels, probably pre-segmented */
  if (nx/mx < 0.875) minval = (float)(umin-tiny);

  return minval;
}


/* apply minval to list of functional slices */
void VApplyMinval(VAttrList list,VFloat minval)
{
  VImage src=NULL;
  VAttrListPosn posn;
  int r,c,j,n,m;
  double sum,nx,u;


  /* apply minval */
  m = 0;
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    VGetAttrValue (& posn, NULL,VImageRepn, & src);
    if (VImageNRows(src) < 2) continue;

    n = VImageNBands(src);

    for (r=0; r<VImageNRows(src); r++) {
      for (c=0; c<VImageNColumns(src); c++) {

	sum = nx = 0;
	for (j=0; j<n; j++) {
	  u = (double)VGetPixel(src,j,r,c);
	  sum += u;
	  nx++;
	}

	if (sum/nx < minval) {
	  for (j=0; j<n; j++) VSetPixel(src,j,r,c,0.0);
	}
	m++;
      }
    }
  }
  if (m < 1) VError("No voxels above threshold 'minval', use lower threshold");
}


void VApplyMinvalNlists(VAttrList *list,int nlists,float minval)
{
  VAttrListPosn posn;
  int slice,row,col,i,k,nrows=0,ncols=0,ntimesteps=0;

  /* alloc src image */
  int nslices = VAttrListNumImages(list[0]);
  VImage **src = (VImage **) VCalloc(nlists,sizeof(VImage *));
  int *nt = (int *) VCalloc(nlists,sizeof(int));
  int ntt=0;
  i = nrows = ncols = 0;
  for (k=0; k<nlists; k++) {
    src[k] = (VImage *) VCalloc(nslices,sizeof(VImage));
    ntt = 0;
    i = 0;
    for (VFirstAttr (list[k], & posn); VAttrExists (& posn); VNextAttr (& posn)) {
      if (VGetAttrRepn (& posn) != VImageRepn) continue;
      if (i >= nslices) VError(" inconsistent number of slices, %d %d %d",k,i,nslices);
      VGetAttrValue (& posn, NULL,VImageRepn, & src[k][i]);
      if (VImageNBands(src[k][i]) > ntt) ntt = VImageNBands(src[k][i]);
      if (VImageNRows(src[k][i])  > nrows)  nrows = VImageNRows(src[k][i]);
      if (VImageNColumns(src[k][i]) > ncols) ncols = VImageNColumns(src[k][i]);
      i++;
    }
    ntimesteps += ntt;
    nt[k] = ntt;
  }

  double sum=0,nx=0,mean=0;
  for (slice=0; slice<nslices; slice++) {
    for (row=0; row<nrows; row++) {
      for (col=0; col<ncols; col++) {

	sum = nx = 0;
	for (k=0; k<nlists; k++) {
	  for (i=0; i<nt[k]; i++) {
	    sum += VGetPixel(src[k][slice],i,row,col);
	    nx++;
	  }
	}
	mean = sum/nx;

	if (mean < minval) {
	  for (k=0; k<nlists; k++) {
	    for (i=0; i<nt[k]; i++) {
	      VSetPixel(src[k][slice],i,row,col,0.0);
	    }
	  }
	}
      }
    }
  }
}
