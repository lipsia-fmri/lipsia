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


void VFreqFilter(VAttrList list,VFloat high,VFloat low,VBoolean stop,VFloat sharp)
{
  VString str1=NULL;
  double *in=NULL;
  fftw_complex *out;
  fftw_plan p1,p2;
  int i,j,nc,b,r,c,tail=0;
  float freq=0,alpha=0;
  double *highp=NULL, *lowp=NULL;
  double x_high, x_low, x;
  double tiny=1.0e-6;


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
  if (tr > 100) tr /= 1000.0;  /* convert to seconds */
  if (tr < 0.01) VError("FreqFilter: implausible TR (%f seconds)",tr);
  fprintf(stderr," image dimensions: %d x %d x %d,  nt= %d,  TR= %.3f secs\n",nslices,nrows,ncols,nt,tr);  


  /* alloc memory */
  double nx = (double)nt;
  tail = 50;
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
  

  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {

	double u = VGetPixel(src[b],0,r,c);
	if (fabs(u) < tiny) continue;

	for (j=0; j<tail; j++) {
	  in[j] = (double)VGetPixel(src[b],tail-j,r,c); 
	}
	for (j=tail; j<nt-tail; j++) {
	  in[j] = (double)VGetPixel(src[b],j-tail,r,c); 
	}
	for (j=nt-tail; j<nt; j++) {
	  in[j] = (double)VGetPixel(src[b],2*nt-3*tail-2-j,r,c); 
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
	for (i=0; i<nt; i++) in[i] /= nx;

	for (j=tail; j<nt-tail; j++) {
	  VSetPixel(src[b],j-tail,r,c,in[j]);
	}
      }
    }
  }
}

