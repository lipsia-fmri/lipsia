/*! \file
  2D/3D median filter


\par Author:
Gabriele Lohmann, MPI-CBS
*/



/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/mu.h>


/* From the standard C library: */
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>

#define ABS(x) ((x) > 0 ? (x) : -(x))


/*!
\fn VImage VMedianImage3d (VImage src, VImage dest, int dim, VBoolean ignore)
\param src      input image
\param dest     output image 
\param dim      kernel size (3,5,7...)
\param ignore   whether to ignore zero voxels
*/
VImage 
VMedianImage3d (VImage src, VImage dest, int dim, VBoolean ignore)
{
  int nbands,nrows,ncols;
  int i,len,len2,b,r,c,bb,rr,cc,b0,b1,r0,r1,c0,c1,d=0;
  gsl_vector *vec=NULL;
  double u,tiny=1.0e-8;

  if (dim%2 == 0) VError("VMedianImage3d: dim (%d) must be odd",dim);

  nrows  = VImageNRows (src);
  ncols  = VImageNColumns (src);
  nbands = VImageNBands (src);
  if (nbands <= d) VError("VMedianImage3d: number of slices too small (%d)",nbands);


  d = dim/2;
  len = dim * dim * dim;
  len2 = len/2;
  
  vec = gsl_vector_calloc(len);

  dest = VCopyImage(src,dest,VAllBands);
  VFillImage(dest,VAllBands,0);

  for (b=d; b<nbands-d; b++) {
    for (r=d; r<nrows-d; r++) {
      for (c=d; c<ncols-d; c++) {

	u = VGetPixel(src,b,r,c);
	
	i = 0;
	b0 = b-d;
	b1 = b+d;
	for (bb=b0; bb<=b1; bb++) {

	  r0 = r-d;
	  r1 = r+d;
	  for (rr=r0; rr<=r1; rr++) {

	    c0 = c-d;
	    c1 = c+d;
	    for (cc=c0; cc<=c1; cc++) {
	      u = VGetPixel(src,bb,rr,cc);
	      if (ABS(u) > tiny && ignore == TRUE) {
		gsl_vector_set(vec,i,u);
		i++;
	      }
	    }
	  }
	  
	  if (i < len2) continue;
	  gsl_sort(vec->data,vec->stride,(size_t)i);
	  u = gsl_stats_median_from_sorted_data(vec->data,vec->stride,i);
	  VSetPixel(dest,b,r,c,u);
	}
      }
    }
  }
  gsl_vector_free(vec);
  return dest;
}




/*!
\fn VImage VMedianImage2d (VImage src, VImage dest, int dim, VBoolean ignore)
\param src      input image
\param dest     output image 
\param dim      kernel size
\param ignore   whether to ignore zero voxels
*/
VImage VMedianImage2d (VImage src, VImage dest, int dim, VBoolean ignore)
{
  int nbands,nrows,ncols;
  int i,len,len2,b,r,c,rr,cc,d=0;
  gsl_vector *vec=NULL;
  double u=0;

  if (dim%2 == 0) VError("VMedianImage2d: dim (%d) must be odd",dim);

  nrows  = VImageNRows (src);
  ncols  = VImageNColumns (src);
  nbands = VImageNBands (src);

  d = dim/2;
  len = dim * dim;
  len2 = len/2;
  vec = gsl_vector_calloc(len);

  dest   = VCopyImage(src,dest,VAllBands);

  for (b=0; b<nbands; b++) {
    for (r=d; r<nrows-d; r++) {
      for (c=d; c<ncols-d; c++) {
	gsl_vector_set_zero(vec);

	u = VGetPixel(src,b,r,c);

	i = 0;
	for (rr=r-d; rr<=r+d; rr++) {
	  for (cc=c-d; cc<=c+d; cc++) {
	    u = VGetPixel(src,b,rr,cc);
	    if (ABS(u) > 0 && ignore == TRUE) {
	      gsl_vector_set(vec,i,u);
	      i++;
	    }
	  }
	}

	if (i < len2) continue;
	gsl_sort(vec->data,vec->stride,(size_t)i);
	u = gsl_stats_median_from_sorted_data(vec->data,vec->stride,i);
	VSetPixel(dest,b,r,c,u);
      }
    }
  }
  gsl_vector_free(vec);
  return dest;
}
