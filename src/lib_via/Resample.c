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


VImage VTriLinearResample(VImage src,VImage dest,double **transform,double *shift,double *fixpoint,
			  int dst_nbands,int dst_nrows,int dst_ncolumns)
{
  int   b,r,c;
  float bp,rp,cp,bx,rx,cx;
  int   sx, sy, sz;   /* origin of subcube    */
  float px, py, pz;   /* fractions of subcube */
  float qx, qy, qz;   /* fractions of subcube */
  int   lx, ly, lz;     /* lengths */
  int   ox, oy, oz;     /* offsets */

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


#define GetValues(type) \
{ \
  type *src_pp; \
  src_pp = (type *) VPixelPtr (src, sz, sy, sx); \
  val += (float) pz * py * px * (*src_pp); src_pp += ox; \
  val += (float) pz * py * qx * (*src_pp); src_pp += oy; \
  val += (float) pz * qy * px * (*src_pp); src_pp += ox; \
  val += (float) pz * qy * qx * (*src_pp); src_pp += oz; \
  val += (float) qz * py * px * (*src_pp); src_pp += ox; \
  val += (float) qz * py * qx * (*src_pp); src_pp += oy; \
  val += (float) qz * qy * px * (*src_pp); src_pp += ox; \
  val += (float) qz * qy * qx * (*src_pp); \
  VPixel(dest,b,r,c,type) = (type)val;		\
}


  /* Determines the value of each pixel in the destination image: */
  for (b=0; b<dst_nbands; b++) {
    bx = (float) b - shift[2];

    for (r=0; r<dst_nrows; r++) {
      rx = (float) r - shift[1];

      for (c=0; c<dst_ncolumns; c++) {
	cx = (float) c - shift[0];

	cp = ainv[0][0] * cx + ainv[0][1] * rx + ainv[0][2] * bx;
	rp = ainv[1][0] * cx + ainv[1][1] * rx + ainv[1][2] * bx;
	bp = ainv[2][0] * cx + ainv[2][1] * rx + ainv[2][2] * bx;

	cp += fixpoint[0];
	rp += fixpoint[1];
	bp += fixpoint[2];

	/* fprintf(stderr," %3d %3d %3d,  %f %f %f\n",b,r,c,bp,rp,cp); */

	if (bp < 0 || bp > src_nbands) continue;
	if (rp < 0 || rp > src_nrows) continue;
	if (cp < 0 || cp > src_ncols) continue;

	/* compute origin of subcube */
	sx = (int) (cp);
	sy = (int) (rp);
	sz = (int) (bp);

	/* check subcube */
	if ((sx < -1) || (sx > src_ncols  - 1)) continue;
	if ((sy < -1) || (sy > src_nrows  - 1)) continue;
	if ((sz < -1) || (sz > src_nbands - 1)) continue;

	/* compute fractions of subcube */
	qx = cp - sx; px = 1 - qx;
	qy = rp - sy; py = 1 - qy;
	qz = bp - sz; pz = 1 - qz;

	/* compute lengths and offsets */
	lx = 1;
	ly = src_ncols;
	lz = src_nrows * src_ncols;
	if (sx == -1) {sx = 0; lx = 0;};
	if (sy == -1) {sy = 0; ly = 0;};
	if (sz == -1) {sz = 0; lz = 0;};
	if (sx == src_ncols  - 1) lx = 0;
	if (sy == src_nrows  - 1) ly = 0;
	if (sz == src_nbands - 1) lz = 0;
	ox = lx;
	oy = ox + ly - 2 * lx;
	oz = oy + lz - 2 * ly;
	double val = 0;

	switch(repn) {
	case VShortRepn:
	  GetValues(VShort);
	  break;
	case VUByteRepn:
	  GetValues(VUByte);
	  break;
	case VFloatRepn:
	  GetValues(VFloat);
	  break;
	case VSByteRepn:
	  GetValues(VSByte);
	  break;
	default:
	  VError(" illegal pixel repn");
	}
      }
    }
  }

  return dest;
}
