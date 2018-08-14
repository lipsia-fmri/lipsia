/*
**
** G.Lohmann, MPI-KYB, May 2018
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "viaio/Vlib.h"
#include "viaio/VImage.h"
#include "viaio/mu.h"
#include "viaio/option.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

extern void VCheckImage(VImage src);
extern VImage VFilterGauss3d (VImage src,VImage dest,double *sigma);

void Gauss4d(VAttrList list,double *sigma)
{
  /* get pointers to image data */
  double u=0;
  int j,b,r,c;
  int nrows=0,ncols=0,nt=0;
  int nslices = VAttrListNumImages(list);
  VImage *src = VAttrListGetImages(list,nslices);
  VImageDimensions(src,nslices,&nt,&nrows,&ncols);
  fprintf(stderr," image dimensions:  %d x %d x %d,  nt: %d\n",nslices,nrows,ncols,nt);

  
  VImage tmp = VCreateImage(nslices,nrows,ncols,VFloatRepn);
  VFillImage(tmp,VAllBands,0);
  VImage dest = VCreateImage(nslices,nrows,ncols,VFloatRepn);
  VFillImage(dest,VAllBands,0);


  for (j=0; j<nt; j++) {
    if (j%5==0) fprintf(stderr," %5d  of  %d\r",j,nt);
    for (b=0; b<nslices; b++) {
      for (r=0; r<nrows; r++) {
	for (c=0; c<ncols; c++) {
	  u = VGetPixel(src[b],j,r,c);
	  VSetPixel(tmp,b,r,c,u);
	}
      }
    }
    dest = VFilterGauss3d (tmp,dest,sigma);
    
    for (b=0; b<nslices; b++) {
      for (r=0; r<nrows; r++) {
	for (c=0; c<ncols; c++) {
	  u = VGetPixel(dest,b,r,c);
	  VSetPixel(src[b],j,r,c,u);
	}
      }
    }
  }
  fprintf(stderr,"\n");
  VDestroyImage(tmp);
  VDestroyImage(dest);
}
