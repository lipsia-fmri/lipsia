/*
** LISA, 2nd level using design matrix
**
** G.Lohmann, July 2018
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
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_math.h>


#define ABS(x) ((x) > 0 ? (x) : -(x))
#define SQR(x) ((x)*(x))

extern void SubtractMean(gsl_matrix *X);
extern gsl_matrix *PseudoInv(gsl_matrix *A,gsl_matrix *B);
extern gsl_matrix *XRead2ndLevel(VString);
extern double t2z(double t,double df);


double EstimateVariance(gsl_matrix *XInv,gsl_vector *con)
{
  int nbeta = XInv->size1;
  double var=0;

  gsl_matrix *bcov = gsl_matrix_calloc(nbeta,nbeta);
  gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,XInv,XInv,0.0,bcov);

  gsl_vector *tmp  = gsl_vector_alloc(nbeta);
  gsl_blas_dgemv(CblasTrans,1.0,bcov,con,0.0,tmp);
  gsl_blas_ddot (tmp,con,&var);
  double sigma = sqrt(var);
  gsl_matrix_free(bcov);
  gsl_vector_free(tmp);

  return sigma;
}


/* 
** get effective degrees of freedom 
*/
double DoF(gsl_matrix *X,gsl_matrix *XInv,double *xtrace)
{
  int i,nimages = X->size1;

  gsl_matrix *R = gsl_matrix_calloc (nimages,nimages);
  gsl_matrix *P = gsl_matrix_calloc (nimages,nimages);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,X,XInv,0.0,P);
  gsl_matrix_set_identity(R);
  gsl_matrix_sub(R,P);

  double trace = 0;
  for (i=0; i<nimages; i++)
    trace += gsl_matrix_get(R,i,i);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,R,R,0.0,P);
  gsl_matrix_free(R);

  double trace2 = 0;
  for (i=0; i<nimages; i++) {
    trace2 += gsl_matrix_get(P,i,i);
  }
  gsl_matrix_free(P);

  (*xtrace) = trace;
  double df = (trace*trace) / trace2;
  return df;
}



/*
** general linear regression, 2nd level
*/
void GLM2(VImage *src,gsl_matrix *X,gsl_vector *contrast,int *permtable,int *signtable,int signswitch,VImage dest)
{
  int i,j,ip,b,r,c;
  double t=0,z=0,sum=0,tsigma=0,d=0,err=0,var=0,u=0,tiny=1.0e-6;

  int ncols   = VImageNColumns(src[0]);
  int nrows   = VImageNRows(src[0]);
  int nslices = VImageNBands(src[0]);
  int nimages = X->size1;
  VFillImage(dest,VAllBands,0);


  /* permutation or sign switching applied to design matrix */
  gsl_matrix *XP = gsl_matrix_calloc(X->size1,X->size2);
  gsl_matrix_memcpy(XP,X);

  if (signswitch >= 0) {   /* sign switching  (one-sample test) */
    for (i=0; i<X->size1; i++) {
      if (signtable[i] < 0) {
	u = gsl_matrix_get(X,i,signswitch);
	gsl_matrix_set(XP,i,signswitch,-u);
      }
    }
  }

  else {       /*  permutations  */
    for (i=0; i<X->size1; i++) {
      ip = permtable[i];
      for (j=0; j<X->size2; j++) {
	if (fabs(contrast->data[j]) > 0) {  /* only permute columns with nonzero contrast */
	  u = gsl_matrix_get(X,i,j);
	  gsl_matrix_set(XP,ip,j,u);
	}
      }
    }
  }

  
  /* pseudo inverse */
  gsl_matrix *XInv = PseudoInv(XP,NULL);


  /* get DoF and variance */
  double trace = 0;
  double df = DoF(XP,XInv,&trace);
  double sigma = EstimateVariance(XInv,contrast);
  if (df < tiny) VError(" Zero degrees of freedom");
  if (sigma < tiny) VError(" No variance in design/contrast");


  /* main loop */
  gsl_set_error_handler_off();
  gsl_vector *y    = gsl_vector_alloc (nimages);
  gsl_vector *yz   = gsl_vector_alloc (nimages);
  gsl_vector *beta = gsl_vector_alloc (X->size2);

  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
  
	int nonzero=0;
	for (j=0; j<nimages; j++) {
	  y->data[j] = VGetPixel(src[j],b,r,c);
	  if (fabs(y->data[j]) > 0) nonzero++;
	}
	if (nonzero < nimages-2) continue;


	/* compute beta's */
	gsl_blas_dgemv(CblasNoTrans,1.0,XInv,y,0.0,beta);


	/* residuals */
	gsl_blas_dgemv(CblasNoTrans,1.0,XP,beta,0.0,yz);
	err = 0;
	for (j=0; j<nimages; j++) {
	  d = y->data[j] - yz->data[j];
	  err += d*d;
	}
	

	/* get z-value */
	t = z = sum = 0;
	gsl_blas_ddot (beta,contrast,&sum);
	if (fabs(sum) < tiny) continue;
	var = err / trace;
	tsigma = sqrt(var) * sigma;
	if (tsigma > tiny) {
	  t = sum / tsigma;
	  z = t2z(t,df);
	  if (sum < 0) z = -z;
	}
	VPixel(dest,b,r,c,VFloat) = z;
      }
    }
  }

  gsl_vector_free(y);
  gsl_vector_free(yz);
  gsl_vector_free(beta);
  gsl_matrix_free(XP);
  gsl_matrix_free(XInv);
}
