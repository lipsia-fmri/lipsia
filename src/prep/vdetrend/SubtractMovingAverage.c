/*
** detrending, subtract moving average
** G.Lohmann, 2018, MPI-KYB
*/
#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>


extern void VGetStats(double *data,int nt,int i0,double *ave,double *sigma);

void VSubtractMovingAverage(VAttrList list,float window,int del)
{ 
  int slice,row,col,i,j,k;
  double ave=0,sigma=0,sum=0,nx=0;

  
  /* get image dimensions */
  int nslices = VAttrListNumImages(list);
  VImage *src = VAttrListGetImages(list,nslices);
  int nrows  = VImageNRows (src[0]);
  int ncols  = VImageNColumns (src[0]);
  int ntimesteps = VImageNBands (src[0]);
  if (nslices < 1) VError(" no slices");
  double smax = VPixelMaxValue(src[0]);
  double smin = VPixelMinValue(src[0]);

  double tr=0;
  VFloat xtr=-1;
  if (VGetAttr(VImageAttrList(src[0]),"repetition_time",NULL,VFloatRepn,&xtr) == VAttrFound) {
    tr = (double)xtr;
    if (tr > 100) tr /= 1000.0;  /* convert to seconds */
  }
  else {
    VError(" TR information missing");
  }
  if (tr < TINY) VError("FreqFilter: implausible or missing TR" );
  int wn = (int)(0.5 * window / tr);
  fprintf(stderr," Subtract moving mean,  window size: %.4f sec (%d time points)\n",window,2*wn+1);


  double *x = (double *)VCalloc(ntimesteps,sizeof(double));
  double *y = (double *)VCalloc(ntimesteps,sizeof(double));


  /* ini dest image */
  int flag=0;
  VImage *dst = (VImage *)VCalloc(nslices,sizeof(VImage));
  if (VPixelRepn(src[0]) != VFloatRepn) {
    flag=1;
    for (i=0; i<nslices; i++) {
      dst[i] = VCreateImage(ntimesteps,nrows,ncols,VFloatRepn);
      VFillImage(dst[i],VAllBands,0);
      VCopyImageAttrs (src[i], dst[i]);
    }
  }
  
  
  /* subtract moving average */
  for (slice=0; slice<nslices; slice++) {
    if (slice%5==0)fprintf(stderr," slice  %5d  of  %d\r",slice,nslices);
    for (row=0; row<nrows; row++) {
      for (col=0; col<ncols; col++) {
    
	for (j=0; j<ntimesteps; j++) {
	  x[j] = VGetPixel(src[slice],j,row,col);
	}


	/* get mean,sigma over entire time series */
	VGetStats(x,ntimesteps,del,&ave,&sigma);
	if (sigma < TINY) continue;


	/* subtract moving median */
	for (j=0; j<ntimesteps; j++) {
	  sum = nx = 0;
	  for (k=j-wn; k<=j+wn; k++) {
	    if (k < 0 || k >= ntimesteps) continue;
	    if (del > 0 && k <= del) continue;
	    sum += x[k];
	    nx++;
	  }
	  if (nx < 1.0) continue;

	  if (del == 0 || j > del) {
	    y[j] = x[j] - sum/nx + ave;
	  }
	  else if (del > 0 && j <= del) {
	    y[j] = ave;
	  }
	  else
	    fprintf(stderr," j=%d, del= %d\n",j,del);

	  if (y[j] < smin) y[j] = smin;
	  if (y[j] >= smax) y[j] = smax;
	  if (flag == 0) VSetPixel(src[slice],j,row,col,y[j]);
	  else VPixel(dst[slice],j,row,col,VFloat) = y[j];
	}
      }
    }
  }
  fprintf(stderr,"\n");
  VFree(x);
  VFree(y);
  if (flag == 0) return;


  /* update attr list */
  VAttrListPosn posn;
  i=0;
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    VDestroyImage(src[i]);
    if (i >= nslices) break;
    VSetAttrValue (& posn, NULL,VImageRepn,dst[i]);
    i++;
  } 
}
