#include <viaio/VImage.h>
#include <viaio/Vlib.h>

/* From the standard C library: */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "../cylutils/cyl.h"

extern void XWriteOutput(VImage image,VAttrList geolist,char *filename);

/* first order spatial autocorrelation */
double XSpatialAutocorrelation(VImage src,float xthr)
{
  size_t nn=0;
  int b,r,c;
  int nslices = VImageNBands(src);
  int nrows = VImageNRows(src);
  int ncols = VImageNColumns(src);
  double u=0,v0=0,v1=0,v2=0;

   
  double sum=0,rho=0;
  for (b=0; b<nslices-1; b++) {
    for (r=0; r<nrows-1; r++) {
      for (c=0; c<ncols-1; c++) {
	u = VPixel(src,b,r,c,VFloat);
	if (u < xthr) continue;
	v0 = VPixel(src,b+1,r,c,VFloat);
	v1 = VPixel(src,b,r+1,c,VFloat);
	v2 = VPixel(src,b,r,c+1,VFloat);
	if (v0 > xthr) { rho += u*v0; 	sum += u*u;  nn++; }
	if (v1 > xthr) { rho += u*v1; 	sum += u*u;  nn++; }
	if (v2 > xthr) { rho += u*v2; 	sum += u*u;  nn++; }
      }
    }
  }
  if (nn > 2) rho /= sum;
  if (rho > 1) rho = 1;
  return rho;
}

double SpatialAutocorrelation(Cylinders *cyl,size_t cid,gsl_vector *residuals)
{
  int i,k,b,r,c,bmin,rmin,cmin,bmax,rmax,cmax;
  size_t n=cyl->addr[cid]->size;

  bmin = rmin = cmin = 99999;
  bmax = rmax = cmax = 0;
  for (i=0; i<n; i++) {
    k = cyl->addr[cid]->data[i];
    b = gsl_matrix_int_get(cyl->xmap,k,0);
    r = gsl_matrix_int_get(cyl->xmap,k,1);
    c = gsl_matrix_int_get(cyl->xmap,k,2);
    if (b < bmin) bmin = b;
    if (r < rmin) rmin = r;
    if (c < cmin) cmin = c;
    if (b > bmax) bmax = b;
    if (r > rmax) rmax = r;
    if (c > cmax) cmax = c;
  }

  size_t nslices = bmax-bmin+2;
  size_t nrows = rmax-rmin+2;
  size_t ncols = cmax-cmin+2;
  VImage tmp = VCreateImage(nslices,nrows,ncols,VFloatRepn);
  float xmin = -9999;
  float xthr = -9998;
  VFillImage(tmp,VAllBands,xmin);
 
  for (i=0; i<n; i++) {
    k = cyl->addr[cid]->data[i];
    b = gsl_matrix_int_get(cyl->xmap,k,0);
    r = gsl_matrix_int_get(cyl->xmap,k,1);
    c = gsl_matrix_int_get(cyl->xmap,k,2);
    VPixel(tmp,b-bmin,r-rmin,c-cmin,VFloat) = residuals->data[i];
  }
  /*
  XWriteOutput(tmp,NULL,"rho.v");
  exit(0);
  */
 
  double rho = XSpatialAutocorrelation(tmp,xthr);
  VDestroyImage(tmp);
  return rho;
}
