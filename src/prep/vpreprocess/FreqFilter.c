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

#define ABS(x) ((x) < 0 ? -(x) : (x))


void VFreqFilter(VAttrList list,VFloat high,VFloat low,VBoolean stop,VFloat sharp)
{
  VAttrListPosn posn;
  VImage src=NULL;
  VString str=NULL, str1=NULL;
  double *in=NULL, sum=0;
  fftw_complex *out;
  fftw_plan p1,p2;
  int i,j,k,n,nc,r,c,nrows,ncols,not=0,tail=0;
  float freq,alpha,tr=0;
  double *highp=NULL, *lowp=NULL;
  double x_high, x_low, x;
  double eps=0.00000001;


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


  /* get image dimensions */
  n = nrows = ncols = 0;
  str = VMalloc(100);
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    VGetAttrValue (& posn, NULL,VImageRepn, & src);
    if (VImageNBands(src) > n) n = VImageNBands(src);
    if (VImageNRows(src) > nrows) nrows = VImageNRows(src);
    if (VImageNColumns(src) > ncols) ncols = VImageNColumns(src);
    
    /* Get TR */
    if (VGetAttr(VImageAttrList(src),"repetition_time",NULL,VFloatRepn,&tr) == VAttrFound)
      tr /= 1000.0;
    else {
      if (VGetAttr (VImageAttrList (src), "MPIL_vista_0", NULL,VStringRepn, (VPointer) & str) == VAttrFound) {
	sscanf(str," repetition_time=%f",&tr);
	tr /= 1000.0;
      } else
	VError(" TR unknown");
    }
  }
  if (tr < 0.01 || tr > 100) VError(" illegal TR %f",tr);
 

  /* alloc memory */
  tail = 50;
  if (n<tail) tail=n-2;
  n  += 2 * tail;
  nc  = (n / 2) + 1;
  in  = (double *)fftw_malloc(sizeof(double) * n);
  out = fftw_malloc (sizeof (fftw_complex ) * nc);
  for (i=0; i<n; i++) in[i] = 0;


  /* make plans */
  p1 = fftw_plan_dft_r2c_1d (n,in,out,FFTW_ESTIMATE);
  p2 = fftw_plan_dft_c2r_1d (n,out,in,FFTW_ESTIMATE);

  /* repetition time */
  alpha = (double)n * tr;

  /* filter function */
  if (sharp > 0) {
    highp = (double *)malloc(sizeof(double) * nc);
    lowp  = (double *)malloc(sizeof(double) * nc);
    for (i=1; i <nc; i++) {
      highp[i] = 1.0 / (1.0 +  exp( (alpha/high -(double)i)*sharp )   );
      lowp[i]  = 1.0 / (1.0 +  exp( ((double)i - alpha/low)*sharp )   );
      /* Butterworth example
	 lowp[i]  = ( (gain*gain) / (1 + pow(((double)i*low/alpha),2*(double)order)) );
	 highp[i] = ( (gain*gain) / (1 + pow((alpha/high/(double)i),2*(double)order)) ); */
    }
  }
  
  k = -1;
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    VGetAttrValue (& posn, NULL,VImageRepn, & src);
    k++;

    if (VImageNRows(src) < 2) continue;

    if (k%2 == 0) fprintf(stderr," slice %4d\r",k);
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {


	/* get data */
	sum=0;

	for (j=0; j<tail; j++) {
	  in[j] = (double)VGetPixel(src,tail-j,r,c); 
	  sum += in[j];
	}
	for (j=tail; j<n-tail; j++) {
	  in[j] = (double)VGetPixel(src,j-tail,r,c); 
	  sum += in[j];
	}
	for (j=n-tail; j<n; j++) {
	  in[j] = (double)VGetPixel(src,2*n-3*tail-2-j,r,c); 
	  sum += in[j];
	}

	/*  fft */
	if ( ABS(sum)>eps ) {
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
	      if (stop) x = ABS(1-x);
	      
	      out[i][0] *= x;
	      out[i][1] *= x;
	      
	    } else {
	      
	      /* hard thresholding */	    
	      freq = 1.0/(double)i * alpha;     /* 1/frequency */
	      
	      if ((!stop && (freq < low || (freq > high && high>0)))
		  || (stop && !(freq < low || (freq > high && high>0))))
		out[i][0] = out[i][1] = 0;
	    }
	  }
	  
	  
	  /* inverse fft */
	  fftw_execute(p2);
	  for (i=0; i<n; i++) in[i] /= (double)n;

	  for (j=tail; j<n-tail; j++) 
	    VSetPixel(src,j-tail,r,c,in[j]);
	  
	}
	else not++;
      }
    }
  }
  /* fprintf(stderr,"\n %d voxels excluded (mean == 0)\n",not); */
}

