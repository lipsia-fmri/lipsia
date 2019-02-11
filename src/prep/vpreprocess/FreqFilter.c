/****************************************************************
 *
 * vpreprocess: FreqFilter.c
 *
 * Copyright (C) Max Planck Institute 
 * for Human Cognitive and Brain Sciences, Leipzig
 *
 * <lipsia@cbs.mpg.de>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * 
 * $Id: FreqFilter.c 3181 2008-04-01 15:19:44Z karstenm $
 *
 *****************************************************************/

#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <fftw3.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>

extern gsl_matrix *DataMap(VImage *src,VImage mask,size_t nvox,int);

gsl_matrix *DataMap(VImage *src,VImage mask,size_t nvox,int tail)
{
  int j,b,r,c;
  int nrows=0,ncols=0,nt0=0;
  double tiny=1.0e-6;

  int nslices = VImageNBands(mask);
  VImageDimensions(src,nslices,&nt0,&nrows,&ncols);

  int nt = nt0 + 2*tail;
  gsl_matrix *A = gsl_matrix_calloc(nvox,nt);

  
  size_t n=0;
  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {

	double u = VGetPixel(src[b],0,r,c);
	if (fabs(u) < tiny) continue;

	for (j=0; j<tail; j++) {
	  u = (double)VGetPixel(src[b],tail-j,r,c); 
	  gsl_matrix_set(A,n,j,u);
	}
	for (j=tail; j<nt-tail; j++) {
	  u = (double)VGetPixel(src[b],j-tail,r,c);
	  gsl_matrix_set(A,n,j,u);
	}
	for (j=nt-tail; j<nt; j++) {
	  u = (double)VGetPixel(src[b],2*nt-3*tail-2-j,r,c); 
	  gsl_matrix_set(A,n,j,u);
	}
	n++;
      }
    }
  }
  return A;
}


void VFreqFilter(VAttrList list,VFloat high,VFloat low,VBoolean stop,VFloat sharp)
{
  VString str1=NULL;
  double *in=NULL;
  fftw_complex *out;
  fftw_plan p1,p2;
  int i,j,nc,b,r,c;
  float freq=0,alpha=0;
  double *highp=NULL, *lowp=NULL;
  double x_high, x_low, x;
  double u=0,v=0,tiny=1.0e-6;


  /* dialog messages */
  str1 = (VString)VMalloc(sizeof(char)*4);
  if (stop) strcpy(str1, "stop");
  else strcpy(str1, "pass");
  if (low>0 && high>0) 
    fprintf(stderr," band %s filter: 1/%.1f Hz and 1/%.1f Hz\n",str1,high,low);
  else {
    if (high>0) fprintf(stderr," high %s filter: 1/%.1f Hz\n",str1,high);
    if (low>0)  fprintf(stderr," low %s filter: 1/%.1f Hz\n",str1,low);
  }


  
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
  }  
  if (tr > 100) tr /= 1000.0;  /* convert to seconds */
  if (tr < 0.01) VError("FreqFilter: implausible TR (%f seconds)",tr);
  fprintf(stderr," image dimensions: %d x %d x %d,  nt= %d,  TR= %.3f secs\n",nslices,nrows,ncols,nt,tr);  


  /* alloc memory */
  double nx = (double)nt;
  int tail = 50;
  if (nt<tail) tail=nt-2;
  nt  += 2 * tail;
  nc  = (nt / 2) + 1;
  in  = (double *)fftw_malloc(sizeof(double) * nt);
  out = fftw_malloc (sizeof (fftw_complex ) * nc);
  for (i=0; i<nt; i++) in[i] = 0;


  /* make plans */
  p1 = fftw_plan_dft_r2c_1d (nt,in,out,FFTW_ESTIMATE);
  p2 = fftw_plan_dft_c2r_1d (nt,out,in,FFTW_ESTIMATE);


  /* repetition time */
  alpha = (double)nt * tr;


  /* filter function */
  if (sharp > 0) {
    highp = (double *)malloc(sizeof(double) * nc);
    lowp  = (double *)malloc(sizeof(double) * nc);
    for (i=1; i <nc; i++) {
      highp[i] = 1.0 / (1.0 +  exp( (alpha/high -(double)i)*sharp )   );
      lowp[i]  = 1.0 / (1.0 +  exp( ((double)i - alpha/low)*sharp )   );
    }
  }

  /* create ROI mask */
  VImage mask = VCreateImage(nslices,nrows,ncols,VBitRepn);
  VFillImage(mask,VAllBands,0);
  size_t nvox=0;
  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	double u = VGetPixel(src[b],0,r,c);
	if (fabs(u) < tiny) continue;
	VPixel(mask,b,r,c,VBit) = 1;
	nvox++;
      }
    }
  }

  /* map to double values */
  gsl_matrix *A = DataMap(src,mask,nvox,tail);
  double zmax = gsl_matrix_max(A);
  double zmin = gsl_matrix_min(A);


  /* main process */
  size_t n=0;
  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	
	if (VPixel(mask,b,r,c,VBit) == 0) continue;
	for (i=0; i<nt; i++) {
	  in[i] = gsl_matrix_get(A,n,i);
	}

	/*  fft */
	fftw_execute(p1);
	  
	/* remove specified frequencies */
	for (i=1; i <nc; i++) {
	    
	  if (sharp > 0) {
	      
	    /* highpass */
	    if (high > 0) x_high = highp[i];
	    else x_high = 1.0;
	      
	    /* lowpass */
	    if (low > 0) x_low = lowp[i];
	    else x_low = 1.0;
	      
	    x = x_high + x_low - 1.0;
	    if (stop) x = fabs(1.0-x);
	      
	    out[i][0] *= x;
	    out[i][1] *= x;
	      
	  }
	  else {
	      
	    /* hard thresholding */	    
	    freq = 1.0/(double)i * alpha;     /* 1/frequency */
	      
	    if ((!stop && (freq < low || (freq > high && high>0)))
		|| (stop && !(freq < low || (freq > high && high>0))))
	      out[i][0] = out[i][1] = 0;
	  }
	}
	  
	/* inverse fft */
	fftw_execute(p2);
	for (i=0; i<nt; i++) {
	  gsl_matrix_set(A,n,i,in[i]/nx);
	}
	n++;
      }
    }
  }

  
  /* rescale */
  VRepnKind repn = VPixelRepn(src[0]);
  double smax = VPixelMaxValue(src[0]);
  double smin = VPixelMinValue(src[0]);

  double ymax = 32000;
  if (repn == VFloatRepn || repn == VDoubleRepn) ymax = 1000.0;
  
  for (b=0; b<nslices; b++) {
    VFillImage(src[b],VAllBands,0);
  }
  zmax = gsl_matrix_max(A);
  zmin = gsl_matrix_min(A);


  /* write output */
  n=0;
  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	if (VPixel(mask,b,r,c,VBit) == 0) continue;
	
	/* map to output image */
	for (j=tail; j<nt-tail; j++) {
	  u = gsl_matrix_get(A,n,j);
	  if (fabs(u) > 0.0001) u = (u-zmin)/(zmax-zmin);
	  else u = 0;
	  if (u > 1.0) u = 1.0;
	  if (u < 0) u = 0;
	  
	  v = u*ymax;
	  if (v > smax) v = smax;
	  if (v < smin) v = smin;
	  VSetPixel(src[b],j-tail,r,c,v);
	}
	n++;
      }
    }
  }
  VDestroyImage(mask);
  gsl_matrix_free(A);
}

