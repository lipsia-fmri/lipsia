#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int DNorm(double *data,size_t n)
{
  size_t j;
  double s1=0,s2=0,mean,sd,u,tiny=1.0e-6;
  double nx = (double)n;
  
  for (j=0; j<n; j++) {
    u = data[j];
    s1 += u;
    s2 += u*u;
  }
  mean = s1/nx;
  sd = sqrt((s2 - nx * mean * mean) / (nx - 1.0));
  if (sd < tiny) return -1;
  for (j=0; j<n; j++) {
    u = data[j];
    data[j] = (u-mean)/sd;
  }
  return 1;
}

double DCorr(double *data1,double *data2,size_t n)
{
  size_t j;
  double sum=0,nx=(double)n;
  for (j=0; j<n; j++) {
    sum += data1[j]*data2[j];
  }
  return sum/nx;
}



void DMN(VImage map,float *A,size_t nvox)
{
  size_t i,j,k;
  int b,r,c;

  fprintf(stderr," DMN...\n");
  
  VAttrList list = VReadAttrList("dmn.v",0L,TRUE,FALSE);
  VImage dmn = VReadImage(list);
  if (dmn == NULL) VError(" no DMN found");
  VImage dst = VCreateImageLike(dmn);
  VFillImage(dst,VAllBands,0);

  /*
  size_t i0=0;
  for (i=0; i<nvox; i++) {
    b = VPixel(map,0,0,i,VInteger);
    r = VPixel(map,0,1,i,VInteger);
    c = VPixel(map,0,2,i,VInteger);
    if (c==31 && r==25 && b==33) i0 = i;
    if (c==20 && r==18 && b==32) i0 = i;
  }
  i = i0;
  
  for (j=0; j<nvox; j++) {
    b = VPixel(map,0,0,j,VInteger);
    r = VPixel(map,0,1,j,VInteger);
    c = VPixel(map,0,2,j,VInteger);
    if (j < i) k = j+i*(i+1)/2;
    else k=i+j*(j+1)/2;
    VPixel(dst,b,r,c,VFloat) = A[k];
  }
  */

  
  double *data1 = (double *)VCalloc(nvox,sizeof(double));
  double *data2 = (double *)VCalloc(nvox,sizeof(double));

  for (j=0; j<nvox; j++) {
    b = VPixel(map,0,0,j,VInteger);
    r = VPixel(map,0,1,j,VInteger);
    c = VPixel(map,0,2,j,VInteger);
    data1[j] = VGetPixel(dmn,b,r,c);
  }
  DNorm(data1,nvox);

  for (i=0; i<nvox; i++) {
    if (i%100 == 0) fprintf(stderr," i: %7lu  of %lu\r",i,nvox);
    k=0;
    for (j=0; j<nvox; j++) {
      b = VPixel(map,0,0,j,VInteger);
      r = VPixel(map,0,1,j,VInteger);
      c = VPixel(map,0,2,j,VInteger);

      if (j < i) k = j+i*(i+1)/2;
      else k=i+j*(j+1)/2;
      data2[j] = A[k];
    }
    if (DNorm(data2,nvox) < 0) continue;

    b = VPixel(map,0,0,i,VInteger);
    r = VPixel(map,0,1,i,VInteger);
    c = VPixel(map,0,2,i,VInteger);
    VPixel(dst,b,r,c,VFloat) = DCorr(data1,data2,nvox);
  }
  VFree(data1);
  VFree(data2);
  
  
  VAttrList out_list = VCreateAttrList();
  VAttrList geolist = VGetGeoInfo(list);
  VSetGeoInfo(geolist,out_list);
  VAppendAttr(out_list,"image",NULL,VImageRepn,dst);
  FILE *fp = fopen("dst.v","w");
  VWriteFile (fp, out_list);
  fclose(fp);
  fprintf(stderr,"\n DMN, ok.\n");
}
