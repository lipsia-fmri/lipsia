/*
** GLM (general linear modeling)
**
** G.Lohmann
*/
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

extern gsl_matrix *PseudoInv(gsl_matrix *A,gsl_matrix *B);
extern void printmat(gsl_matrix *R,char *str);
extern void printvec(gsl_vector *x,char *str);

extern void prewhite(gsl_matrix_float** pinvX,gsl_matrix_float** invM,gsl_matrix_float* X,int);

extern void dfs(VFloat** Dfs,gsl_matrix_float* pinvX,gsl_matrix_float* X,gsl_matrix_float* con,int);

extern void whitecov(VImage rho_vol,gsl_matrix_float* Y,gsl_matrix_float* invM,
		     gsl_matrix_float* pinvX,gsl_matrix_float* X,int,int);

extern void whitecov2(VImage effect_image, VImage rho_vol,
		      gsl_matrix_float* Y, gsl_matrix_float* X, gsl_matrix_float* con,
		      VFloat* Dfs, int, int slice);


void CopySlice(gsl_matrix *Data,gsl_matrix_float *Y,VImage map,int slice)
{
  size_t nvox;
  size_t i=0,j=0;
  double u=0;
  int ncols = VPixel(map,0,3,2,VShort);

  gsl_matrix_float_set_zero(Y);

  for (nvox=0; nvox < Data->size1; nvox++) {
    int b = VPixel(map,0,0,nvox,VShort);
    if (b != slice) continue;
    int r = VPixel(map,0,1,nvox,VShort);
    int c = VPixel(map,0,2,nvox,VShort);
    i = c + ncols*r;
    if (i >= Y->size2) continue; 

    for (j=0; j<Data->size2; j++) {
      u = gsl_matrix_get(Data,nvox,j);
      gsl_matrix_float_set(Y,j,i,(float)u);
    }
  }
}

gsl_matrix_float *GslFloatCpy(gsl_matrix *X)
{
  int i,j;
  double u;
  gsl_matrix_float *Z = gsl_matrix_float_calloc(X->size2,X->size1);
  for (i=0; i<X->size1; i++) {
    for (j=0; j<X->size2; j++) {
      u = gsl_matrix_get(X,i,j);
      gsl_matrix_float_set(Z,j,i,(float)u);
    }
  }
  return Z;
}


/*
** general linear model (GLM) using whitening
*/
void VWhiteGLM(gsl_matrix *Data,VImage map,gsl_matrix *X0,gsl_vector *con,int numlags,VImage zmap)
{
  int i,slice;
  int numcon=1;
  int nslices = VPixel(map,0,3,0,VShort);
  int nrows = VPixel(map,0,3,1,VShort);
  int ncols = VPixel(map,0,3,2,VShort);
  int npix = (nrows*ncols);
  int ntimesteps = Data->size2;

  gsl_set_error_handler_off();
  VFillImage(zmap,VAllBands,0);


  /* contrast vector */
  gsl_matrix_float *contrast = gsl_matrix_float_alloc(1,con->size);
  for (i=0; i<con->size; i++) gsl_matrix_float_set(contrast,0,i,(float)con->data[i]);


  /* copy design matrix */
  gsl_matrix_float *X = GslFloatCpy(X0);


  /* initialize rho_vol */
  VImage rho_vol = VCreateImage(npix,nslices,numlags,VFloatRepn);
  VFillImage(rho_vol,VAllBands,0);


  /* prewhitening */
  gsl_matrix_float *pinvX=NULL,*invM=NULL;
  prewhite(&pinvX,&invM,X,numlags);


  /* degrees of freedom */
  VFloat *dfs_gsl = (VFloat *)VMalloc(sizeof(VFloat)* (numcon+1));
  dfs(&dfs_gsl,pinvX,X,contrast,numlags); 


  /* first pass */
  gsl_matrix_float *Y = gsl_matrix_float_calloc(ntimesteps,npix);
  for (slice=0; slice<nslices; slice++) {
    CopySlice(Data,Y,map,(int)slice);
    whitecov(rho_vol,Y,invM,pinvX,X,numlags,slice);
  }


  /* second pass */
  VImage dest = VCreateImage(nslices,nrows,ncols,VFloatRepn);
  VFillImage(dest,VAllBands,0);
  for (slice=0; slice<nslices; slice++) {
    CopySlice(Data,Y,map,(int)slice);
    whitecov2(dest,rho_vol,Y,X,contrast,dfs_gsl,numlags,slice);
  }


  /* output */
  zmap = VCopyImagePixels(dest,zmap,VAllBands);


  /* cleanup */
  gsl_matrix_float_free(contrast);
  gsl_matrix_float_free(X);
  gsl_matrix_float_free(Y);
  gsl_matrix_float_free(pinvX);
  gsl_matrix_float_free(invM);
  VFree(dfs_gsl);
  VDestroyImage(rho_vol);
  VDestroyImage(dest);
}
