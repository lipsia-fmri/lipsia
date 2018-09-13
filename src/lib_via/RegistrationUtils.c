/* 
** get transformation matrix for co-registration 
** Translation, rotation and one scaling parameter (7 degrees of freedom)
**
** G.Lohmann, 2015
*/
#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <viaio/os.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

#define ABS(x) ((x) > 0 ? (x) : -(x))

#define PI 3.14159265358979323846


double **MatrixAlloc(int n,int m)
{
  int i,j;
  double **A = (double **)VCalloc(n,sizeof(double));
  for(i=0; i<n; i++) {
    A[i]=(double *)VCalloc(m,sizeof(double));
  }
  for (i=0; i<n; i++) {
    for (j=0; j<m; j++) {
      A[i][j] = 0.0;
    }
  }
  return A;
}

void MatrixFree(double **A,int n,int m)
{
  int i;
  for(i=0; i<n; i++) VFree(A[i]);
  VFree(A);
}

/* matrix mult of 3x3 rotation matrices */
void MatrixMult3x3(double **A, double **B, double **C)
{
  int i,j,k,n=3;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      double sum = 0;
      for (k=0; k<n; k++) {
        sum += A[i][k] * B[k][j];
      }
      C[i][j] = sum;
    }
  }
}


/* matrix ini */
void MatrixIni3x3(double **A)
{
  int i,j,n=3;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      A[i][j] = 0.0;
    }
    for (i=0; i<n; i++) A[i][i] = 1.0;
  }
}


/* 3x3 matrix times vector */
void MatrixVector3x3(double **A, double *x,double *y)
{
  int i,j,n=3;
  for (i=0; i<n; i++) {
    double sum = 0;
    for (j=0; j<n; j++) {
      sum += A[i][j]*x[j];
    }
    y[i] = sum;
  }
}


void MatrixCopy3x3(double **A,double **B)
{
  int i,j,n=3;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      B[i][j] = A[i][j];
    }
  }
}


void MatrixScale3x3(double **R,double sx,double sy,double sz)
{
  int i;
  for (i=0; i<3; i++) {
    R[i][0] *= sx;
    R[i][1] *= sy;
    R[i][2] *= sz;
  }
}

void printmat(VImage src)
{
  int i,j;
  fprintf(stderr,"printmat:\n");
  for (i=0; i<VImageNRows(src); i++) {
    for (j=0; j<VImageNColumns(src); j++) {
      float u = VGetPixel(src,0,i,j);
      fprintf(stderr," %8.3f",u);
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");
}



void xprintmat(double **A,int n,int m)
{
  int i,j;
  fprintf(stderr,"printmat:\n");
  for (i=0; i<n; i++) {
    for (j=0; j<m; j++) {
      fprintf(stderr," %8.3f",A[i][j]);
    }
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");
}

void VCheckImage(VImage src)
{
  VAttrList out_list = VCreateAttrList();
  VAppendAttr(out_list,"image",NULL,VImageRepn,src);
  FILE *out_file = fopen("test.v","w");
  VWriteFile (out_file, out_list);
}


double clamp(double v) 
{
  v = fmod(v, 2.0*PI);
  return v>=0 ? v : 2.0*PI+v;
}

void EulerAngles(double **R,double *ea)
{
  int i;

  double u = -R[0][2];
  if (u < -1.0) u = u + ceil(u);
  if (u >  1.0) u = u - floor(u);
  if (u < -1.0) u = -1.0;
  if (u >  1.0) u =  1.0;
  ea[1] = asin(u);

  double xcos = cos(ea[1]);
  if (ABS(xcos)<0.001*PI/180.0) {   /* Gimball lock? */
    ea[0] = atan2(R[1][0], R[1][1]);
    ea[2] = 0;
  } 
  else {
    ea[0] = atan2(R[1][2], R[2][2]);
    ea[2] = atan2(R[0][1], R[0][0]);
  }
  for (i=0; i<3; i++) ea[i] = fmod(ea[i],2.0*PI);
  /*   for (i=0; i<3; i++) ea[i] = clamp(ea[i]); */
}

void Euler2RotationMatrix(double **R,double *ea)
{
  double c1 = cos(ea[0]);
  double c2 = cos(ea[1]);
  double c3 = cos(ea[2]);
  double s1 = sin(ea[0]);
  double s2 = sin(ea[1]);
  double s3 = sin(ea[2]);

  R[0][0] = c2*c1;
  R[1][0] = -c3*s1+s3*s2*c1;
  R[2][0] = s3*s1+c3*s2*c1;
  R[0][1] = c2*s1;
  R[1][1] = c3*c1+s3*s2*s1;
  R[2][1] = -s3*c1+c3*s2*s1;
  R[0][2] = -s2;
  R[1][2] = s3*c2;
  R[2][2] = c3*c2;

  double qfac = -1.0;
  R[0][2] *= qfac;
  R[1][2] *= qfac;
  R[2][2] *= qfac;
}



void Normalize(VImage src)
{
  int b,r,c;
  double tiny = 1.0e-6;
  if (VPixelRepn(src) == VUByteRepn || VPixelRepn(src) == VShortRepn 
      || VPixelRepn(src) == VSByteRepn)tiny = 1.0;

  double sum1=0,sum2=0,nx=0;
  for (b=0; b<VImageNBands(src); b++) {
    for (r=0; r<VImageNRows(src); r++) {
      for (c=0; c<VImageNColumns(src); c++) {
	double u = VGetPixel(src,b,r,c);
	if (u >= tiny) {
	  sum1 += u;
	  sum2 += u*u;
	  nx++;
	}
      }
    }
  }
  double mean = sum1/nx;
  double sigma = sqrt((sum2 - nx * mean * mean) / (nx - 1.0));

  double xmax = VPixelMaxValue(src);
  double xmin = VPixelMinValue(src);
  double mean0 = (xmax + xmin)/2.0;
  double sig0 = xmax*0.2;
  fprintf(stderr," normalize, mean,sd= %f %f\n",mean0,sig0);

  for (b=0; b<VImageNBands(src); b++) {
    for (r=0; r<VImageNRows(src); r++) {
      for (c=0; c<VImageNColumns(src); c++) {
	double u = VGetPixel(src,b,r,c);
	if (u >= tiny) {
	  u = (u-mean)/sigma;
	  u = u*sig0 + mean0;
	  if (u > xmax) u = xmax;
	  if (u < xmin) u = xmin;
	  VSetPixel(src,b,r,c,u);
	}
      }
    }
  }
}


/* set up transformation info */
void ReadTransformation(VImage transimage,double **trans,double *shift,
			double *ref_pixdim,double *ref_dim,double *src_pixdim,double *src_dim)			
{
  int i,j;
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      trans[i][j] = VGetPixel(transimage,0,i,j);
    }
  }
  for (i=0; i<3; i++)
    shift[i] = VGetPixel(transimage,0,i,3);

  for (i=0; i<3; i++)
    ref_pixdim[i] = VGetPixel(transimage,0,i,4);

  for (i=0; i<3; i++) 
    ref_dim[i] = VGetPixel(transimage,0,i,5);

  for (i=0; i<3; i++)
    src_pixdim[i] = VGetPixel(transimage,0,i,6);

  for (i=0; i<3; i++) 
    src_dim[i] = VGetPixel(transimage,0,i,7);
}


void VSetResolution(VImage transimage,double **trans,double newres,double *shift,double *scale)
{
  /* resolution of reference image */
  double resx0 = VGetPixel(transimage,0,0,6);
  double resy0 = VGetPixel(transimage,0,1,6);
  double resz0 = VGetPixel(transimage,0,2,6);

  /* change resolution */
  scale[0] = resx0/newres;
  scale[1] = resy0/newres;
  scale[2] = resz0/newres;

  int i;
  MatrixScale3x3(trans,scale[0],scale[1],scale[2]);
  for (i=0; i<3; i++) shift[i] *= scale[i];
}
