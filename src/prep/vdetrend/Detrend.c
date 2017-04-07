/*
** detrending, fit polynomial or order 1 (linear), 3 or 5
** G.Lohmann, Aug 2016, MPI-KYB
*/
#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>


#define ABS(x) ((x) < 0 ? -(x) : (x))


void VGetStats(double *data,int nt,int i0,double *ave,double *sigma)
{
  int i;
  double sum1,sum2,nx;
  sum1 = sum2 = 0;
  for (i=i0; i<nt; i++) {
    double u = data[i];
    sum1 += u;
    sum2 += u*u;
  }
  nx = (double)nt;
  double mean = sum1/nx;
  (*ave) = mean;
  (*sigma) = sqrt((sum2 - nx * mean * mean) / (nx - 1.0));
}


void VDetrend(VAttrList list,VFloat minval,VShort type,VShort i0)
{ 
  VImage tmp=NULL;
  VAttrListPosn posn;
  int row,col,i,j,ntimesteps,nrows,ncols,nslices;


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
  if (nslices < 1) VError(" no slices in image data");
  if (ntimesteps < 3) VError(" no timesteps in image data");


  int p = 4;
  if (type == 2) p = 6;
  gsl_vector *y = gsl_vector_calloc(ntimesteps);
  gsl_vector *c = gsl_vector_calloc(p);
  gsl_matrix *X = gsl_matrix_calloc(ntimesteps,p);
  for (j=0; j<ntimesteps; j++) {
    double u = (double)j;
    gsl_matrix_set(X,j,0,1);
    gsl_matrix_set(X,j,1,u);
    gsl_matrix_set(X,j,2,u*u);
    gsl_matrix_set(X,j,3,u*u*u);
    if (type == 2) {
      gsl_matrix_set(X,j,4,u*u*u*u);
      gsl_matrix_set(X,j,5,u*u*u*u*u);
    }
  }
  gsl_multifit_linear_workspace *workspace = gsl_multifit_linear_alloc ((size_t)ntimesteps,(size_t)p);
  gsl_matrix *cov = gsl_matrix_calloc(p,p);


  double *x = (double *) VCalloc(ntimesteps,sizeof(double));
  for (j=0; j<ntimesteps; j++) x[j] = (float)j;

  gsl_vector *w = gsl_vector_calloc(ntimesteps);
  for (j=0; j<ntimesteps; j++) w->data[j] = (float)1.0;
  for (i=0; i<=i0; i++) w->data[i] = 0.01; /* ignore initial time steps */


  double smin = VRepnMinValue(VShortRepn);
  double smax = VRepnMaxValue(VShortRepn);


  /* remove baseline drift */
  double c0=0,c1=0,c2=0,c3=0,c4=0,c5=0,cov00,cov01,cov11,chisq;
  for (i=0; i<nslices; i++) {
    if (i%5==0)fprintf(stderr," slice  %5d  of  %d\r",i,nslices);
    for (row=0; row<nrows; row++) {
      for (col=0; col<ncols; col++) {
    
	for (j=0; j<ntimesteps; j++) {
	  y->data[j] = (double) VPixel(src[i],j,row,col,VShort);
	}
	if (y->data[0] < minval) {
	  for (j=0; j<ntimesteps; j++) VPixel(src[i],j,row,col,VShort) = 0;
	  goto skip;
	}
	double mean=0,sigma=0;
	VGetStats(y->data,ntimesteps,i0,&mean,&sigma);

	/* remove linear trend */
	if (type == 0) {
	  gsl_fit_wlinear (x,1,w->data,1,y->data,1,j,&c0,&c1,&cov00,&cov01,&cov11,&chisq);
	  for (j=0; j<i0; j++) y->data[j] = 0;
	  for (j=i0; j<ntimesteps; j++) {
	    y->data[j] = y->data[j] - (c0 + c1*x[j]);
	  }
	}

	/* remove cubic polynomial */
	if (type == 1) {
	  gsl_multifit_wlinear (X,w,y,c,cov,&chisq,workspace);
	  c0 = c->data[0];
	  c1 = c->data[1];
	  c2 = c->data[2];
	  c3 = c->data[3];	  

	  for (j=0; j<i0; j++) y->data[j] = 0;
	  for (j=i0; j<ntimesteps; j++) {
	    double u = x[j];
	    y->data[j] = y->data[j] - (c0 + c1*u + c2*u*u + c3*u*u*u);
	  }
	}

	/* remove polynomial of order 5 */
	if (type == 2) {
	  gsl_multifit_wlinear (X,w,y,c,cov,&chisq,workspace);
	  c0 = c->data[0];
	  c1 = c->data[1];
	  c2 = c->data[2];
	  c3 = c->data[3];
	  c4 = c->data[4];
	  c5 = c->data[5];
	  for (j=0; j<i0; j++) y->data[j] = 0;
	  for (j=i0; j<ntimesteps; j++) {
	    double u = x[j];
	    double u2 = u*u;
	    double u3 = u2*u;
	    y->data[j] = y->data[j] - (c0 + c1*u + c2*u2 + c3*u3 + c4*u2*u2 + c5*u2*u3);
	  }
	}

	double a = 1.25;
	double sum=0,nx=0;
	for (j=i0; j<ntimesteps; j++) {
	  sum += a*y->data[j] + mean;
	  nx++;
	}
	if (nx > 0) sum /= nx;
	
	for (j=0; j<i0; j++) VPixel(src[i],j,row,col,VShort) = sum;
	for (j=i0; j<ntimesteps; j++) {
	  double u = a*y->data[j] + mean;
	  VShort v = (int)(u+0.5);
	  if (v > smax) v = smax;
	  if (v < smin) v = smin;
	  VPixel(src[i],j,row,col,VShort) = v;
	}
      skip: ;
      }
    }
  }
  fprintf(stderr,"\n");
}
