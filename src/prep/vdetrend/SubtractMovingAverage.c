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
  int slice,row,col,j,k;
  double ave=0,sigma=0,sum=0,nx=0,tiny=1.0e-6;


  /* get image dimensions */
  int nslices = VAttrListNumImages(list);
  VImage *src = VAttrListGetImages(list,nslices);
  int nrows  = VImageNRows (src[0]);
  int ncols  = VImageNColumns (src[0]);
  int ntimesteps = VImageNBands (src[0]);
  if (nslices < 1) VError(" no slices");
  double smax = VPixelMaxValue(src[0]);
  double smin = VPixelMinValue(src[0]);


  VAttrList geolist = VGetGeoInfo(list);
  if (!geolist) VError(" no geoinfo");
  double *D = VGetGeoPixdim(geolist,NULL);
  double tr = D[4]/1000.0;
  int wn = (int)(0.5 * window / tr);
  fprintf(stderr," Subtract moving mean,  window size: %.4f sec (%d time points)\n",window,2*wn+1);


  double *x = (double *)VCalloc(ntimesteps,sizeof(double));
  double *y = (double *)VCalloc(ntimesteps,sizeof(double));


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
	if (sigma < tiny) continue;


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
	  VSetPixel(src[slice],j,row,col,y[j]);
	}
      }
    }
  }
  fprintf(stderr,"\n");
}
