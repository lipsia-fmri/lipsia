/*
** centrality maps using spectral coherence
**
** get cross-periodogram via multiplication of FFTs
**
** Ref: T. Subba Rao,
** On the Cross Periodogram of a Stationary Gaussian Vector Process,
** The Annals of Mathematical Statistics, Vol. 38, No. 2 (Apr., 1967),
**  pp. 593-597.  (2239172.pdf)
**
** Implementation follows SAS, see
**  http://www.tau.ac.il/cc/pages/docs/sas8/ets/chap17/
**
** G.Lohmann, September 2009
*/
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_cblas.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265

#define SQR(x) ((x) * (x))
#define ABS(x) ((x) > 0 ? (x) : -(x))


/* re-implementation of cblas_sspmv,  cblas_sspmv causes problems */
void my_sspmv(float *A,float *x,float *y,int n)
{
  size_t i,j,k,kk;
  float tmp1=0,tmp2=0;
  
  for (j=0; j<n; j++) y[j] = 0;
  kk = k = 0;
  for (j=0; j<n; j++) {
    tmp1 = x[j];
    tmp2 = 0;
    k = kk;
    for (i=0; i<j; i++) {
      y[i] += tmp1 * A[k];
      tmp2 += A[k] * x[i];
      k++;
    }
    y[j] += tmp1*A[kk+j] + tmp2;
    kk += (j+1);
  }
}


int Normalize(double *data,int n)
{
  int i;
  double sum,tiny=1.0e-6;
  
  sum = 0;
  for (i=0; i<n; i++) sum += data[i];
  sum /= (double)n;
  if (sum < tiny) return -1;
  
  for (i=0; i<n; i++) data[i] -= sum;
  return 1;
}

double Tukey(double x,double wx)
{
  if (ABS(x) >= wx) return 0;
  else return 0.5*(1.0 + cos((PI*x)/wx));
}

int GetFFtIndex(int n,double tr,double wavelength)
{
  int k0;
  double nx,cycle;

  nx = (double)n;
  cycle = wavelength / tr;
  k0 = (int)(nx/cycle + 0.5);
  return k0;
}


double *Weights(int n,int k0,int p)
{
  int k;
  double *w=NULL;
  double kx,wx,sum;
  
  w  = (double *) VCalloc(n,sizeof(double));
  wx = (double) p;
  kx = 0;
  for (k=0; k<p; k++) {
    w[k] = Tukey(kx,wx);
    kx++;
  }
  sum = 0;
  for (k=0; k<p; k++) sum += w[k];
  sum *= (2.0*PI);
  for (k=0; k<p; k++) w[k] /= sum;
  return w;
}


void Tables(gsl_matrix *tablecos,gsl_matrix *tablesin)
{
  int i,k,n;
  double kx,ix,nx,w;
  
  n = tablecos->size1;
  nx = (double)n;
  kx = 0;
  for (k=0; k<n; k++) {
    w = 2.0*PI*kx/nx;
    ix = 0;
    for (i=1; i<n; i++) {
      gsl_matrix_set(tablecos,k,i,cos(w*ix));
      gsl_matrix_set(tablesin,k,i,sin(w*ix));
      ix++;
    }
    kx++;
  }
}		

void FFT(double *in,double *xcos,double *xsin,gsl_matrix *tablecos,gsl_matrix *tablesin)
{
  int i,k,n;
  double sum1,sum2,kx,ix,nx;
  
  n  = tablecos->size1; 
  nx = 1;
  kx = 0;
  for (k=0; k<n; k++) {
    sum1 = sum2 = 0;
    /* w = 2.0*PI*kx/nx; */
    ix = 0;
    for (i=1; i<n; i++) {
      sum1 += in[i]*gsl_matrix_get(tablecos,k,i);
      sum2 += in[i]*gsl_matrix_get(tablesin,k,i);
      ix++;
    }
    xcos[k] = 2.0*sum1/nx;
    xsin[k] = 2.0*sum2/nx;
    kx++;
  }
}



void EigenvectorCentrality(float *A,float *ev,int n)
{
  int i,iter,maxiter;
  float sum,d,nx;
  static float *y=NULL;

  fprintf(stderr," power iteration, n= %d\n",n);
  if (y==NULL) y = (float *) VCalloc(n,sizeof(float));
  if (!y) VError(" err allocating vector ");
  nx = sqrt((double)n);
  for (i=0; i<n; i++) ev[i] = y[i] = 1.0/nx;

  maxiter=50;
  for (iter=0; iter < maxiter; iter++) {

    /* y = Ax,  A symmetric and lower-triangular */
    /* cblas_sspmv(CblasRowMajor,CblasLower,(int)n,1.0f,A,ev,1,1.0f,y,1); */
    my_sspmv(A,ev,y,n);

    sum = 0;
    for (i=0; i<n; i++) sum += y[i]*y[i];
    sum = sqrt(sum);

    d = 0;
    for (i=0; i<n; i++) {
      d += SQR(ev[i] - y[i]/sum);
      ev[i] = y[i]/sum;
    }
    fprintf(stderr," %5d  %f\n",iter,d);
    if (d < 1.0e-6) break;
  }
  for (i=0; i<n; i++) ev[i] *= 100.0;
  fprintf(stderr," ecm done..\n");
}


VImage WriteOutput(VImage src,VImage map,int nslices,int nrows, int ncols, float *ev, int n)
{
  VImage dest=NULL;
  int i,b,r,c;
  
  dest = VCreateImage(nslices,nrows,ncols,VFloatRepn);
  VFillImage(dest,VAllBands,0);
  VCopyImageAttrs (src, dest);
  VSetAttr(VImageAttrList(dest),"modality",NULL,VStringRepn,"conimg");

  for (i=0; i<n; i++) {
    b = VPixel(map,0,0,i,VFloat);
    r = VPixel(map,0,1,i,VFloat);
    c = VPixel(map,0,2,i,VFloat);
    VPixel(dest,b,r,c,VFloat) = ev[i];
  }
  return dest;
}



VAttrList VSpectralECM(VAttrList list,VImage mask,VDouble wavelength,VShort nlags,
		       VShort first,VShort length,VShort type)
{
  VImage dest0=NULL;
  int b,r,c,k0,k,kk,n,p=0;
  size_t i,j;
  size_t ii,m,nvox;
  int last;
  gsl_matrix *matcos=NULL,*matsin=NULL;
  float *A=NULL,*ev=NULL;
  float tr=-1,xtr=-1;
  double *in,*xsin=NULL,*xcos=NULL,*w=NULL;
  double nx,u,v,sumx,freq=0;
  double re,im,x,y,cs,qs,c1,s1,c2,s2;
  double tiny=1.0e-6;
  gsl_matrix *tablecos=NULL,*tablesin=NULL;
  double pi = 2.0*acos(0);


  /* get image dimensions */
  int nrows=0,ncols=0,nt=0;
  int nslices = VAttrListNumImages(list);
  VImage *src = VAttrListGetImages(list,nslices);
  VImageDimensions(src,nslices,&nt,&nrows,&ncols);
  size_t ntimesteps = nt;


  /* get time steps to include */
  if (length < 1) length = ntimesteps-2;
  last = first + length -1;
  if (last >= ntimesteps) last = ntimesteps-1;
  if (first < 0) first = 1;

  n = last - first + 1;
  if (n < 2) VError(" not enough timesteps, nt= %d",n);
  nx = (double)n;


  /* 
  ** get freq range 
  */
  if (VGetAttr (VImageAttrList (src[0]), "repetition_time", NULL,
		VFloatRepn, (VPointer) & xtr) == VAttrFound) {
    tr = (xtr/1000.0);
  }
  if (tr <= 0) VError(" illegal tr: %f",tr);

  k0   = GetFFtIndex((int)n,tr,wavelength); 
  freq = 1.0/wavelength;
  fprintf(stderr,"# TR= %g secs,  wavelength= %g secs,  freq= %.3f Hz,  k0= %d\n",
	  tr,wavelength,freq,k0);
  if (k0 < 1 || k0 >= n/2) VError(" illegal freq, index must be in [1,%d]\n",n/2);


  /* weights */
  if (nlags < 0)
    p = n/15;   /* default */
  else
    p = (int)nlags;
  if (p < 3) VError(" nlags too small, %d",p);
  w = Weights((int)n,k0,p);



  /* count number of voxels */
  nvox = 0;
  for (b=0; b<nslices; b++) {
    if (VImageNRows(src[b]) < 2) continue;
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	if (VGetPixel(mask,b,r,c) < 0.01) continue;
	nvox++;
      }
    }
  }
  fprintf(stderr," nvoxels: %ld,  n= %d,  nlags= %d\n",(long)nvox, (int)n, (int)p);



  /*
  ** voxel maps
  */
  VImage map = VCreateImage(1,4,nvox,VIntegerRepn);
  if (map == NULL) VError(" error allocating addr map");
  VFillImage(map,VAllBands,0);
  VPixel(map,0,3,0,VInteger) = nslices;
  VPixel(map,0,3,1,VInteger) = nrows;
  VPixel(map,0,3,2,VInteger) = ncols;

  double *xmap = (double *)VCalloc(nvox,sizeof(double));

  i = 0;
  for (b=0; b<nslices; b++) {
    if (VImageNRows(src[b]) < 2) continue;
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	if (VGetPixel(mask,b,r,c) < 0.01) continue;
	VPixel(map,0,0,i,VInteger) = b;
        VPixel(map,0,1,i,VInteger) = r;
        VPixel(map,0,2,i,VInteger) = c;
	i++;
      }
    }
  }


  /* alloc memory */
  in   = (double *)VCalloc(n,sizeof(double));
  xsin = (double *)VCalloc(n,sizeof(double));
  xcos = (double *)VCalloc(n,sizeof(double));
  tablecos = gsl_matrix_calloc((int)n,(int)n);
  if (tablecos == NULL) VError(" err alloc cos");
  tablesin = gsl_matrix_calloc((int)n,(int)n);
  if (tablesin == NULL) VError(" err alloc sin");
  Tables(tablecos,tablesin);


  /*
  ** precompute FFT
  */
  matcos = gsl_matrix_calloc(nvox,(int)n);
  if (!matcos) VError(" err alloc matcos");
  matsin = gsl_matrix_calloc(nvox,(int)n);
  if (!matsin) VError(" err alloc matsin");
  Tables(tablecos,tablesin);

   for (i=0; i<nvox; i++) {

    b = VPixel(map,0,0,i,VInteger);
    r = VPixel(map,0,1,i,VInteger);
    c = VPixel(map,0,2,i,VInteger);

    j = 0;
    for (k=first; k<=last; k++)
      in[j++] = (double) VGetPixel(src[b],k,r,c);
    if (Normalize(in,n) < 0) continue;
    
    FFT(in,xcos,xsin,tablecos,tablesin);
    nx = (double)n;
    sumx = 0;
  
    for (k=-p; k<=p; k++) {
      kk = k+k0;
      if (kk < 0) kk = -k;
      if (kk >= n) kk = n-k;
      u = 2.0*(xcos[kk]*xcos[kk] + xsin[kk]*xsin[kk])/nx;
      sumx += w[ABS(k)]*u;
    }
    xmap[i] = sumx;

    for (k=0; k<n; k++) {
      gsl_matrix_set(matcos,i,k,xcos[k]);
      gsl_matrix_set(matsin,i,k,xsin[k]);
    }
  }


  /*
  ** compute spectral coherence matrix
  */
  m = (nvox*(nvox+1))/2;
  fprintf(stderr," matrix computation, n= %ld...\n",(long)nvox);
  A = (float *) VCalloc(m,sizeof(float));
  if (!A) VError(" err allocating matrix");
  memset(A,0,m*sizeof(float));

  ii = 0;
  for (i=0; i<nvox; i++) {
    if (i%100 == 0) fprintf(stderr,"i= %7ld\r",(long)i);
    
    x = xmap[i];

    for (j=0; j<=i; j++) {
      if (i == j) {
	ii++;
	continue;
      }
      y = xmap[j];
      if (x*y < tiny) continue;
 
      cs = qs = 0;
      for (k=-p; k<=p; k++) {

	kk = k+k0;
	if (kk  < 0) kk = -k;
	if (kk >= n) kk = n-k;

	c1 = gsl_matrix_get(matcos,i,kk);
	s1 = gsl_matrix_get(matsin,i,kk);

	c2 = gsl_matrix_get(matcos,j,kk);
	s2 = gsl_matrix_get(matsin,j,kk);
	
	re = 2.0*(c1*c2 + s1*s2)/nx;
	im = 2.0*(c1*s2 - s1*c2)/nx;

	cs += w[ABS(k)] * re;
	qs += w[ABS(k)] * im;
      }
      
      if (type == 0) {  /* spectral coherence */
	v = (cs*cs + qs*qs) / (x*y);
	v = sqrt(v);
      }
      else {            /* phase sync */
	v = atan2(qs,cs)/pi;
	/* v = 1.0 - ABS(v); */
	v = cos(v);
      }

      if (v < tiny) v = tiny;  /* make matrix irreducible */
      A[ii] = v;
      ii++;
    }
  }

  /* free space */
  gsl_matrix_free(matcos);
  gsl_matrix_free(matsin);


  /*
  ** eigenvector centrality
  */
  ev = (float *) VCalloc(nvox,sizeof(float));
  EigenvectorCentrality(A,ev,nvox);
  dest0 = WriteOutput(src[0],map,nslices,nrows,ncols,ev,nvox);
  VSetAttr(VImageAttrList(dest0),"name",NULL,VStringRepn,"eigenvector_centrality");
  VAttrList out_list = VCreateAttrList();
  VAppendAttr(out_list,"image",NULL,VImageRepn,dest0);
  return out_list;
}


VDictEntry TYPDict[] = {
  { "freq", 0 },
  { "sync", 1 },
  { NULL }
};

int main (int argc,char *argv[])
{ 
  static VDouble  wavelength = 10;
  static VShort   nlags      = 10;
  static VShort   first      = 2;
  static VShort   length     = 0;
  static VShort   type       = 0;
  static VString  mask_filename   = "";
  static VOptionDescRec  options[] = {
    {"wavelength", VDoubleRepn,1,(VPointer) &wavelength,VRequiredOpt,NULL,"wavelength in secs"},
    {"first",VShortRepn,1,(VPointer) &first,VOptionalOpt,NULL,"first timestep to use"},
    {"length",VShortRepn,1,(VPointer) &length,VOptionalOpt,NULL,"length of time series to use"},
    {"mask",VStringRepn,1,(VPointer) &mask_filename,VRequiredOpt,NULL,"mask"},
    {"nlags",VShortRepn,1,(VPointer) &nlags,VOptionalOpt,NULL,
     "num lags in cross-corr, use 0 to get default value"},
    {"type",VShortRepn,1,(VPointer) &type,VOptionalOpt,TYPDict,"frequency or phase coherence"}
  };
  FILE *in_file,*out_file;
  VString in_filename=NULL;
  char *prg = "vspectralecm";


  /* parse command line */
  VParseFilterCmdZ (VNumber (options),options,argc,argv,&in_file,&out_file,&in_filename);
  if (type > 2) VError(" illegal type");


  /* read mask */
  VAttrList listm = VReadAttrList(mask_filename,0L,TRUE,FALSE);
  VImage mask = VReadImage(listm);
  if (mask == NULL) VError(" no ROI mask found");


  /* read functional data */
  VAttrList list = VReadAttrListZ(in_file,in_filename,0L,TRUE,FALSE);
  if (list == NULL) VError(" error reading input file %s",in_file);
  VAttrList geolist = VGetGeoInfo(list);


  /*
  ** process
  */
  VAttrList out_list = VSpectralECM(list,mask,wavelength,nlags,first,length,type);


  /* update geoinfo, 4D to 3D */
  if (geolist != NULL) {
    double *D = VGetGeoDim(geolist,NULL);
    D[0] = 3;  /* 3D */
    D[4] = 1;  /* just one timestep */
    VSetGeoDim(geolist,D);
    VSetGeoInfo(geolist,out_list);
  }

  VHistory(VNumber(options),options,prg,&list,&out_list);
  if (! VWriteFile (out_file, out_list)) exit (1);
  fprintf (stderr, "%s: done.\n", argv[0]);
  return 0;
}
