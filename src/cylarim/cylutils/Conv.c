#include <viaio/VImage.h>
#include <viaio/Vlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>

#include "../cylutils/cyl.h"


double gaussian(double x,int i)
{
  double mean[3] = {0.2,0.5,0.8};
  double sigma=0.1;
  double z = exp(-(x-mean[i])*(x-mean[i])/(2.0*sigma*sigma));
  return z;
}


VImage Convert2Repn(VImage src,VImage dest,VRepnKind repn)
{
  size_t b,r,c;
  double u;
  size_t nslices = VImageNBands(src);
  size_t nrows = VImageNRows(src);
  size_t ncols = VImageNColumns(src);

  if (VPixelRepn(src) == repn) { 
    dest = VCopyImage(src,dest,VAllBands);
    return dest;
  }

  dest = VCreateImage(nslices,nrows,ncols,repn);
  VFillImage(dest,VAllBands,0);
  VCopyImageAttrs (src, dest);

  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	
	u = VGetPixel(src,b,r,c);
	
	if (repn == VBitRepn) {
	  if (u < 0.1) VPixel(dest,b,r,c,VBit) = 0;
	  else VPixel(dest,b,r,c,VBit) = 1;
	}
	else if (repn == VUByteRepn) {
	  VPixel(dest,b,r,c,VUByte) = (int)u;
	}
	else {
	  VSetPixel(dest,b,r,c,u);
	}
      }
    }
  }
  return dest;
}



/* create 4D image from several 3D images */
VImage *VCreate4DImage(VImage *src,int nimages,VAttrList geoinfo,VRepnKind repn)
{
  int i;
  int nslices = VImageNBands(src[0]);
  int nrows = VImageNRows(src[0]);
  int ncols = VImageNColumns(src[0]);

  /* dim */
  double *dim = (double *) VCalloc(8,sizeof(double));
  dim[0] = 4;
  dim[1] = ncols;
  dim[2] = nrows;
  dim[3] = nslices;
  dim[4] = nimages;
  VSetGeoDim(geoinfo,dim);

  /* alloc rgb image */
  VImage *dst = (VImage *) VCalloc(nslices, sizeof(VImage));
  for (i=0; i<nslices; i++) {
    dst[i] = VCreateImage(nimages,nrows,ncols,repn);
    VFillImage(dst[i],VAllBands,0);
    VCopyImageAttrs (src[0],dst[i]);
  }
  return dst;
}


/*
  DT_UINT8
  NIFTI_INTENT_VECTOR  intent_code
*/
int Dbl2Int(double u)
{
  if (u < 0) return 0;
  int v = (int)(u*255.0 + 0.5);
  if (v < 0) v = 0;
  if (v > 255) v = 255;
  return v;
}

double NormBeta(VImage *betaimage,int nbeta,char *nonzero,double *bmin)
{
  size_t i,j,k,n;
  size_t npixels = VImageNPixels(betaimage[0]);
  VFloat *pb = NULL;
  VFloat *p0 = NULL;

  /* use intercept if present */
  if (nbeta > 3) p0 = VImageData(betaimage[3]);

  n=0;
  for (i=0; i<npixels; i++) {
    if (nonzero[i] == 1) n++;
  }
  double *tab = (double *)VCalloc(3*n,sizeof(double));

 
  k=0;
  for (i=0; i<npixels; i++) {
    if (nonzero[i] == 0) continue;
    for (j=0; j<3; j++) {
      pb = VImageData(betaimage[j]);
      if (p0 != NULL) tab[k++] = (double)(pb[i]+p0[i]);
      else tab[k++] = (double)(pb[i]);
    }
  }
  gsl_sort(tab,1,k);
  *bmin = gsl_stats_quantile_from_sorted_data(tab,1,k,0.05);
  double bmax = gsl_stats_quantile_from_sorted_data(tab,1,k,0.95);
  VFree(tab);
  return bmax;
}


/* convert beta image to RGB format */
VImage *ConvRGB(VImage *betaimage,int nbeta,VBoolean rgb,VBoolean invert,VAttrList geoinfo,VRepnKind repn)
{
  int j,b,r,c;
  double z[3],zx=0;
  VFloat *pb = NULL;
  size_t npixels = VImageNPixels(betaimage[0]);
  int nslices = VImageNBands(betaimage[0]);
  int nrows = VImageNRows(betaimage[0]);
  int ncols = VImageNColumns(betaimage[0]);

  
  char *nonzero = (char *)VCalloc(npixels,sizeof(char));
  size_t i,n=0;
  for (i=0; i<npixels; i++) {
    float s = 0;
    for (j=0; j<3; j++) {
      pb = VImageData(betaimage[j]);
      if (invert) pb[i] = -pb[i];
      s += fabs(pb[i]);
    }
    if (s > TINY) nonzero[i] = 1;
  }  

  double bmin = 0;
  double bmax = NormBeta(betaimage,nbeta,nonzero,&bmin);

  
  /* rescale to [0,255] for RGB output (unsigned char) */
  VImage *dst = VCreate4DImage(betaimage,3,geoinfo,repn);


  n=0;
  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {

	if (nbeta > 3) zx = VPixel(betaimage[3],b,r,c,VFloat);  /* add intercept if present */
	z[0] = VPixel(betaimage[0],b,r,c,VFloat)+zx;
	z[1] = VPixel(betaimage[1],b,r,c,VFloat)+zx;
	z[2] = VPixel(betaimage[2],b,r,c,VFloat)+zx;

	char ii = nonzero[n++];
	if (ii == 0) continue;

	for (j=0; j<3; j++) {
	  z[j] = (z[j] - bmin)/(bmax-bmin);
	  if (z[j] < 0) z[j] = 0;
	  if (z[j] > 1) z[j] = 1;
	  if (!rgb) z[j] = 1.0 - z[j];   /* CMY */
	  z[j] = z[j]*z[j];
	}
	
	if (VPixelRepn(dst[0]) == VUByteRepn) {
	  for (j=0; j<3; j++) {
	    VPixel(dst[b],j,r,c,VUByte) = Dbl2Int(z[j]);
	  }
	}
	else {
	  for (j=0; j<3; j++) {
	    VPixel(dst[b],j,r,c,VFloat) = z[j];
	  }
	}
      }
    }
  }

  VFree(nonzero);
  return dst;
}


