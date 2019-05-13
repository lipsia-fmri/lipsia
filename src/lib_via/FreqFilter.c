#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>


int MeanSigma(double *data,size_t n,double *ave,double *sigma)
{
  size_t i=0;
  double u=0,nx=0,sum1=0,sum2=0,mean=0,var=0;


  for (i=0; i<n; i++) {
    u = data[i];
    sum1 += u;
    sum2 += u*u;
  }
  nx = (double)n;
  mean = sum1/nx;
  var = (sum2 - nx * mean * mean) / (nx - 1.0);
  if (var > TINY) {
    *ave = mean;
    *sigma = sqrt(var);
    return 1;
  }
  else {
    return -1;
  }
}


int NormalizeTimecourse(double *data,size_t n,double *ave,double *sigma)
{
  size_t i=0;
  double u=0,mean=0,sd=0;

  if (MeanSigma(data,n,&mean,&sd) < 0) return -1;

  for (i=0; i<n; i++) {
    u = data[i];
    data[i] = (u-mean)/sd;
  }
  *ave = mean;
  *sigma = sd;
  return 1;
}


void VFreqFilter(VAttrList list,VFloat highpass,VFloat lowpass,VBoolean stop)
{
  int i,j,k,b,r,c;
  double x=0,freq=0;

  
  /* read functional data */
  int nrows=0,ncols=0,nt=0;
  int nslices = VAttrListNumImages(list);
  VImage *src = VAttrListGetImages(list,nslices);
  VImageDimensions(src,nslices,&nt,&nrows,&ncols);


  /* read header info */
  double *D = (double *)VCalloc(8,sizeof(double));
  for (i=0; i<8; i++) D[i] = 1.0;
  VAttrList geolist = VGetGeoInfo(list);
  if (geolist != NULL) D = VGetGeoPixdim(geolist,D);
  double tr = D[4];

  VFloat xtr=-1;
  if (VGetAttr(VImageAttrList(src[0]),"repetition_time",NULL,VFloatRepn,&xtr) == VAttrFound) {
    tr = (double)xtr;
    if (tr > 100) tr /= 1000.0;  /* convert to seconds */
    if (tr < 0.01) VError("FreqFilter: implausible TR (%f seconds)",tr);
  }
  fprintf(stderr," image dimensions: %d x %d x %d,  nt= %d,  TR= %.3f secs\n",nslices,nrows,ncols,nt,tr);  


  /* ini dest image */
  VImage *dst = (VImage *)VCalloc(nslices,sizeof(VImage));
  for (i=0; i<nslices; i++) {
    dst[i] = VCreateImage(nt,nrows,ncols,VFloatRepn);
    VFillImage(dst[i],VAllBands,0);
    VCopyImageAttrs (src[i], dst[i]);
  }
  double smin = VPixelMinValue(dst[0])+TINY;
  double smax = VPixelMaxValue(dst[0])-TINY;

  

  /* ini fft, zero-padding, add tails */
  size_t tail = 40;
  if (nt < tail) tail = 0;
  size_t ntt = nt + 2*tail;
  
  size_t nfft = 1;
  while (nfft < ntt) nfft *= 2;
  size_t nc  = (nfft / 2) + 1;
  double alpha = 1.0 / ((double)nfft * tr);

  double *data = (double *)VCalloc(nt,sizeof(double));
  double *data_fft = (double *)VCalloc(nfft,sizeof(double));

  
  
  /* bandpass filter function */
  double s=0,xs=1,xf=1;
  double *bandpass = (double *)VCalloc(nc,sizeof(double));
  for (i=1; i<nc; i++) {
    freq = (double)i * alpha;
    s = 1.0/freq;
    xs = xf = 1.0;
    if (lowpass > 0) xf = 1.0 / (1.0 + exp(lowpass-s));
    if (highpass > 0) xs = 1.0 / (1.0 + exp(s-highpass));
    bandpass[i] = xs*xf;
    if (stop == TRUE) bandpass[i] = 1.0-bandpass[i];
  }


  /* main process */
  double ave=0,sigma=0;
  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	
	/* get data */	
	for (j=0; j<nt; j++) data[j] = VGetPixel(src[b],j,r,c);
	if (NormalizeTimecourse(data,nt,&ave,&sigma) < 0) {
	  for (j=0; j<nt; j++) VSetPixel(src[b],j,r,c,0.0);
	  continue;
	}

	
	/* zero padding and tail added */
	for (j=0; j<nfft; j++) data_fft[j] = 0;
	k=tail;
	for (j=0; j<tail; j++) {
	  data_fft[j] = data[k--];
	}
	k=0;
	for (j=tail; j<nt+tail; j++) {
	  data_fft[j] = data[k++];
	}
	k=nt;
	for (j=tail+nt; j<nt+2*tail; j++) {
	  data_fft[j] = data[k--];
	}

		
	/*  fft */
	gsl_fft_real_radix2_transform(data_fft,1,(size_t)nfft);
	
	/* apply bandpass filter */
	for (i=1; i<nc; i++) {
	  x = bandpass[i];
	  data_fft[i] *= x;
	  data_fft[nfft-i] *= x;
	}
	
	/* inverse fft */
	gsl_fft_halfcomplex_radix2_inverse(data_fft,1,(size_t)nfft);
	k=0;
	for (j=0; j<nt; j++) {
	  data[k++] = data_fft[j+tail];
	}
	
	for (j=0; j<nt; j++) {
	  x = data[j]*sigma + ave;
	  if (x < smin) x = smin+TINY;
	  if (x >= smax) x = smax-TINY;
	  VPixel(dst[b],j,r,c,VFloat) = x;
	}
      }
    }
  }
  VFree(data);
  VFree(data_fft);


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
