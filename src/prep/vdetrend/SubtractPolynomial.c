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
  double u,sum1,sum2,nx;
  sum1 = sum2 = 0;
  for (i=i0; i<nt; i++) {
    u = data[i];
    sum1 += u;
    sum2 += u*u;
  }
  nx = (double)nt;
  double mean = sum1/nx;
  (*ave) = mean;
  (*sigma) = sqrt((sum2 - nx * mean * mean) / (nx - 1.0));
}


void VSubtractPolynomial(VAttrList list,VShort type,VShort i0)
{ 
  int slice,row,col,i,j;
  double tiny=1.0e-6;


  /* get image dimensions */
  int nslices = VAttrListNumImages(list);
  VImage *src = VAttrListGetImages(list,nslices);
  int nrows  = VImageNRows (src[0]);
  int ncols  = VImageNColumns (src[0]);
  int ntimesteps = VImageNBands (src[0]);
  if (nslices < 1) VError(" no slices");

  double smax = VPixelMaxValue(src[0]);
  double smin = VPixelMinValue(src[0]);
  
  if (type == 1)
    fprintf(stderr," Subtract linear drift\n");
  if (type == 2)
    fprintf(stderr," Subtract cubic polynmial drift\n");
  if (type == 3)
    fprintf(stderr," Subtract polynmial of order 5\n");


  /* ini regression matrix */
  int p = 4;
  if (type == 3) p = 6;
  gsl_vector *y = gsl_vector_calloc(ntimesteps);
  gsl_vector *c = gsl_vector_calloc(p);
  gsl_matrix *X = gsl_matrix_calloc(ntimesteps,p);
  for (j=0; j<ntimesteps; j++) {
    double u = (double)j;
    gsl_matrix_set(X,j,0,1);
    gsl_matrix_set(X,j,1,u);
    gsl_matrix_set(X,j,2,u*u);
    gsl_matrix_set(X,j,3,u*u*u);
    if (type == 3) {
      gsl_matrix_set(X,j,4,u*u*u*u);
      gsl_matrix_set(X,j,5,u*u*u*u*u);
    }
  }
  gsl_multifit_linear_workspace *workspace = gsl_multifit_linear_alloc ((size_t)ntimesteps,(size_t)p);
  gsl_matrix *cov = gsl_matrix_calloc(p,p);


  double *x = (double *) VCalloc(ntimesteps,sizeof(double));
  for (j=0; j<ntimesteps; j++) x[j] = (double)j;

  gsl_vector *w = gsl_vector_calloc(ntimesteps);
  for (j=0; j<ntimesteps; j++) w->data[j] = (double)1.0;
  for (i=0; i<=i0; i++) w->data[i] = 0.01; /* ignore initial time steps */



  /* remove baseline drift */
  double c0=0,c1=0,c2=0,c3=0,c4=0,c5=0,cov00,cov01,cov11,chisq;
  for (slice=0; slice<nslices; slice++) {
    if (slice%5==0)fprintf(stderr," slice  %5d  of  %d\r",slice,nslices);
    for (row=0; row<nrows; row++) {
      for (col=0; col<ncols; col++) {
    
	for (j=0; j<ntimesteps; j++) {
	  y->data[j] = VGetPixel(src[slice],j,row,col);
	}

	double mean=0,sigma=0;
	VGetStats(y->data,ntimesteps,i0,&mean,&sigma);
	if (sigma < tiny) continue;


	/* remove linear trend */
	if (type == 1) {
	  gsl_fit_wlinear (x,1,w->data,1,y->data,1,j,&c0,&c1,&cov00,&cov01,&cov11,&chisq);
	  for (j=0; j<i0; j++) y->data[j] = 0;
	  for (j=i0; j<ntimesteps; j++) {
	    y->data[j] = y->data[j] - (c0 + c1*x[j]);
	  }
	}

	/* remove cubic polynomial */
	if (type == 2) {
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
	if (type == 3) {
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
	
	for (j=0; j<i0; j++) VSetPixel(src[slice],j,row,col,sum);
	for (j=i0; j<ntimesteps; j++) {
	  double u = a*y->data[j] + mean;
	  if (u > smax) u = smax;
	  if (u < smin) u = smin;
	  VSetPixel(src[slice],j,row,col,u);
	}
      }
    }
  }
  fprintf(stderr,"\n");
}
