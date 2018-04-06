#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/file.h>
#include <viaio/mu.h>

extern void xprintmat(double **A,int n,int m);
extern void printmat(VImage);
extern double **MatrixAlloc(int n,int m);

/* inverse of 3x3 matrix */
double MatrixInverse3x3(double **a,double **ainv)
{
  ainv[0][0] =  a[1][1]*a[2][2] - a[1][2]*a[2][1];
  ainv[1][0] = -a[1][0]*a[2][2] + a[1][2]*a[2][0];
  ainv[2][0] =  a[1][0]*a[2][1] - a[1][1]*a[2][0];

  ainv[0][1] = -a[0][1]*a[2][2] + a[0][2]*a[2][1];
  ainv[1][1] =  a[0][0]*a[2][2] - a[0][2]*a[2][0];
  ainv[2][1] = -a[0][0]*a[2][1] + a[0][1]*a[2][0];

  ainv[0][2] =  a[0][1]*a[1][2] - a[0][2]*a[1][1];
  ainv[1][2] = -a[0][0]*a[1][2] + a[0][2]*a[1][0];
  ainv[2][2] =  a[0][0]*a[1][1] - a[0][1]*a[1][0];

  /* determinant */
  double detA = a[0][0]*ainv[0][0] + a[0][1]*ainv[1][0] + a[0][2]*ainv[2][0]; 
  if (detA == 0) VError(" matrix is singular");

  int i,j;
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      ainv[i][j] /= detA;
    }
  }
  return detA;
}


/* determinant of 3x3 matrix */
double Determinant3x3(double **a)
{
  double ainv[3][3];
  ainv[0][0] =  a[1][1]*a[2][2] - a[1][2]*a[2][1];
  ainv[1][0] = -a[1][0]*a[2][2] + a[1][2]*a[2][0];
  ainv[2][0] =  a[1][0]*a[2][1] - a[1][1]*a[2][0];

  ainv[0][1] = -a[0][1]*a[2][2] + a[0][2]*a[2][1];
  ainv[1][1] =  a[0][0]*a[2][2] - a[0][2]*a[2][0];
  ainv[2][1] = -a[0][0]*a[2][1] + a[0][1]*a[2][0];

  ainv[0][2] =  a[0][1]*a[1][2] - a[0][2]*a[1][1];
  ainv[1][2] = -a[0][0]*a[1][2] + a[0][2]*a[1][0];
  ainv[2][2] =  a[0][0]*a[1][1] - a[0][1]*a[1][0];

  /* determinant */
  double detA = a[0][0]*ainv[0][0] + a[0][1]*ainv[1][0] + a[0][2]*ainv[2][0];
  return detA;
}


VImage VTriLinearResample(VImage src,VImage dest,double **transform,double *shift,
			  int dst_nbands,int dst_nrows,int dst_ncolumns)
{
  int   b,r,c,bb,rr,cc;
  double c000,c100,c001,c101,c010,c110,c011,c111;
  double c00,c01,c10,c11,c0,c1;
  double x,y,z,xd,yd,zd,x0,y0,z0;
  double bx,rx,cx;
  double val;

  VRepnKind repn = VPixelRepn(src);
  int src_nrows  = VImageNRows(src);
  int src_ncols  = VImageNColumns(src);
  int src_nbands = VImageNBands(src);


  /* get inverse transform : */
  double **ainv = MatrixAlloc(3,3);
  double det = MatrixInverse3x3(transform,ainv);
  if (det == 0) VError(" singular matrix");

  /* create output image */
  if (dest == NULL) {
    dest = VCreateImage(dst_nbands,dst_nrows,dst_ncolumns,repn);
  }
  if (! dest) return NULL;
  VFillImage(dest,VAllBands,0);


  /* Determines the value of each pixel in the destination image: */
  for (b=0; b<dst_nbands; b++) {
    bx = (double) b - shift[2];

    for (r=0; r<dst_nrows; r++) {
      rx = (double) r - shift[1];

      for (c=0; c<dst_ncolumns; c++) {
	cx = (double) c - shift[0];


	/* interpolation point */
	x = ainv[0][0] * cx + ainv[0][1] * rx + ainv[0][2] * bx;
	y = ainv[1][0] * cx + ainv[1][1] * rx + ainv[1][2] * bx;
	z = ainv[2][0] * cx + ainv[2][1] * rx + ainv[2][2] * bx;
	

	/* check range */
	cc = (int)x;
	rr = (int)y;
	bb = (int)z;
	if (cc < 1 || cc >= src_ncols-1) continue;
	if (rr < 1 || rr >= src_nrows-1) continue;
	if (bb < 1 || bb >= src_nbands-1) continue;


	/* cube */
	x0 = (int) (x);
	y0 = (int) (y);
	z0 = (int) (z);


	/* distance to corners of cube */
	xd = (x-x0);
	yd = (y-y0);
	zd = (z-z0);


	/* function values at corners of cube */
	c000 = VGetPixel(src,bb,rr,cc);
	c100 = VGetPixel(src,bb,rr,cc+1);
	c001 = VGetPixel(src,bb+1,rr,cc);
	c101 = VGetPixel(src,bb+1,rr,cc+1);
	c010 = VGetPixel(src,bb,rr+1,cc);
	c110 = VGetPixel(src,bb,rr+1,cc+1);
	c011 = VGetPixel(src,bb+1,rr+1,cc);
	c111 = VGetPixel(src,bb+1,rr+1,cc+1);


	/* interpolate */
	c00 = c000*(1.0-xd) + c100*xd;
	c01 = c001*(1.0-xd) + c101*xd;
	c10 = c010*(1.0-xd) + c110*xd;
	c11 = c011*(1.0-xd) + c111*xd;

	c0 = c00*(1.0-yd) + c10*yd;
	c1 = c01*(1.0-yd) + c11*yd;
	val = c0*(1.0-zd) + c1*zd;

	VSetPixel(dest,b,r,c,(double)val);
      }
    }
  }

  return dest;
}
