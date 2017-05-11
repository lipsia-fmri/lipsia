/*
** ECM - eigenvector centrality mapping using linar correlations
**
** G.Lohmann, Mar 2009
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
#include <assert.h>

#include <time.h>

#define NSLICES 10000

#define SQR(x) ((x) * (x))
#define ABS(x) ((x) > 0 ? (x) : -(x))
#define NEG(x) ((x) < 0 ? (-x) : 0)


#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/

/* re-implementation of cblas_sspmv,  cblas_sspmv causes problems */
void my_sspmv(float *A,float *x,float *y,int n)
{
  int i,j;
  size_t k,kk;
  float tmp1=0,tmp2=0;

  kk = k = 0;
  for (j=0; j<n; j++) {
    y[j] = 0;
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


double Correlation(const float * arr1,const float * arr2,int n)
{
  int i;
  double sx,sy,sxx,syy,sxy,rx;
  double tiny=1.0e-4;

  sxx = syy = sxy = sx = sy = 0;
  for (i=0; i<n; i++) {
    const double u = (double)arr1[i];
    const double v = (double)arr2[i];
    sx  += u;
    sy  += v;
    sxx += u*u;
    syy += v*v;
    sxy += u*v;
  }
  const double nx= n;
  const double u = nx*sxx - sx*sx;
  const double v = nx*syy - sy*sy;
  rx = 0;
  if (u*v > tiny)
    rx = (nx*sxy - sx*sy)/sqrt(u*v);

  return rx;
}

void EigenvectorCentrality(float *A,float *ev,int n)
{
  int i,iter,maxiter;
  float sum,d,nx;
  static float *y=NULL;

  fprintf(stderr," power iteration...\n");
  if (y==NULL) y = (float *) VCalloc(n,sizeof(float));
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
    fprintf(stderr," %5d  %f\n",(int)iter,d);
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
    if (b < 0 || r < 0 || c < 0) VError("  WriteOutput: illegal address %d %d %d",b,r,c);
    VPixel(dest,b,r,c,VFloat) = ev[i];
  }
  return dest;
}



VAttrList VECM(VAttrList list,VImage mask,VShort minval,VShort first,VShort length,VShort type)
{
  VAttrList out_list=NULL;
  VAttrListPosn posn;
  VImage src[NSLICES],map=NULL;
  VImage dest=NULL;
  int b,r,c,ntimesteps,nrows,ncols,nslices,nt,last;
  size_t i,j,i1,n,m;
  gsl_matrix_float *mat=NULL;
  float *A=NULL,*ev=NULL;
  float tiny=1.0e-6;

  /*
  ** get image dimensions
  */
  i = ntimesteps = nrows = ncols = 0;
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    VGetAttrValue (& posn, NULL,VImageRepn, & src[i]);
    if (VPixelRepn(src[i]) != VShortRepn) continue;
    if (VImageNBands(src[i]) > ntimesteps) ntimesteps = VImageNBands(src[i]);
    if (VImageNRows(src[i])  > nrows)  nrows = VImageNRows(src[i]);
    if (VImageNColumns(src[i]) > ncols) ncols = VImageNColumns(src[i]);
    i++;
    if (i >= NSLICES) VError(" too many slices");
  }
  nslices = i;

  if (VImageNBands(mask) != nslices) 
    VError(" number of slices inconsistent: %d, mask: %d",nslices,VImageNBands(mask));
  if (VImageNRows(mask) != nrows) 
    VError(" number of rows inconsistent: %d, mask: %d",nrows,VImageNRows(mask));
  if (VImageNColumns(mask) != ncols) 
    VError(" number of columns inconsistent: %d, mask: %d",ncols,VImageNColumns(mask));


  /* get time steps to include */
  if (length < 1) length = ntimesteps-1;
  last = first + length;
  if (last >= ntimesteps) last = ntimesteps-1;
  if (first < 0) first = 0;

  nt = last - first + 1;
  i1 = first;
  if (nt < 2) VError(" not enough timesteps, nt= %d",nt);
  fprintf(stderr,"# ntimesteps: %d, first= %d, last= %d, nt= %d\n",
	  (int)ntimesteps,(int)first,(int)last,(int)nt);


  /* count number of voxels */
  n = 0;
  for (b=0; b<nslices; b++) {
    if (VImageNRows(src[b]) < 2) continue;
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	if (VPixel(src[b],i1,r,c,VShort) < minval) continue;
	if (VGetPixel(mask,b,r,c) < 0.5) continue;
	n++;
      }
    }
  }
  fprintf(stderr," nvoxels: %ld\n",(long)n);



  /*
  ** voxel addresses
  */
  map = VCreateImage(1,5,n,VFloatRepn);
  if (map == NULL) VError(" error allocating addr map");
  VFillImage(map,VAllBands,0);
  VPixel(map,0,3,0,VFloat) = nslices;
  VPixel(map,0,3,1,VFloat) = nrows;
  VPixel(map,0,3,2,VFloat) = ncols;

  i = 0;
  for (b=0; b<nslices; b++) {
    if (VImageNRows(src[b]) < 2) continue;
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	if (VPixel(src[b],i1,r,c,VShort) < minval) continue;
	if (VGetPixel(mask,b,r,c) < 0.5) continue;

	VPixel(map,0,0,i,VFloat) = b;
	VPixel(map,0,1,i,VFloat) = r;
	VPixel(map,0,2,i,VFloat) = c;
	i++;
      }
    }
  }


  /*
  ** avoid casting to float, copy data to matrix
  */
  mat = gsl_matrix_float_calloc(n,nt);
  for (i=0; i<n; i++) {

    b = VPixel(map,0,0,i,VFloat);
    r = VPixel(map,0,1,i,VFloat);
    c = VPixel(map,0,2,i,VFloat);

    float *ptr = gsl_matrix_float_ptr(mat,i,0);
    int k;
    j = 0;
    for (k=first; k<=last; k++) {
      if (j >= mat->size2) VError(" j= %d %d",j,mat->size2);
      if (k >= VImageNBands(src[b])) VError(" k= %d %d",k, VImageNBands(src[b]));
      *ptr++ = (float) VPixel(src[b],k,r,c,VShort);
      j++;
    }
  }


  /*
  ** compute similarity matrix
  */
  m = (n*(n+1))/2;
  fprintf(stderr," matrix computation, n= %ld...\n",(long)n);
  A = (float *) calloc(m,sizeof(float));
  if (!A) VError(" err allocating correlation matrix");
  memset(A,0,m*sizeof(float));
  size_t progress=0;


#pragma omp parallel for shared(progress) private(j) schedule(guided) firstprivate(mat,A)
  for (i=0; i<n; i++) {
    if (i%100 == 0) fprintf(stderr," %d00\r",(int)(++progress));

    const float *arr1 = gsl_matrix_float_const_ptr(mat,i,0);
    for (j=0; j<=i; j++) {
      if (i == j) continue;
      const float *arr2 = gsl_matrix_float_const_ptr(mat,j,0);
      const double v = Correlation(arr1,arr2,nt);
      double u = 0;

      /* make positive */
      switch (type) {

      case 0:
	if (v > tiny) u = v;
	break;

      case 1:
	u = v + 1.0f;
	break;

      case 2:
	u = ABS(v);
	break;

      case 3:
	if (v < -tiny) u = -v;
	break;

      default:
	VError(" illegal type");
      }
      if (u < tiny) u = tiny;

      const size_t k=j+i*(i+1)/2;
      if (k >= m) VError(" illegal addr k= %d, m= %d",k,m);
      A[k] = u;
    }
  }
  fprintf(stderr," matrix done.\n");
  gsl_matrix_float_free(mat);

  /*
  ** eigenvector centrality
  */
  ev = (float *) VCalloc(n,sizeof(float));
  EigenvectorCentrality(A,ev,n);
  dest = WriteOutput(src[0],map,nslices,nrows,ncols,ev,n);
  VSetAttr(VImageAttrList(dest),"name",NULL,VStringRepn,"ECM");

  out_list = VCreateAttrList();
  VAppendAttr(out_list,"image",NULL,VImageRepn,dest);
  return out_list;
}

VDictEntry TYPDict[] = {
  { "pos", 0, 0,0,0,0 },
  { "add", 1, 0,0,0,0 },
  { "abs", 2, 0,0,0,0 },
  { "neg", 3, 0,0,0,0 },
  { NULL, 0,0,0,0,0 }
};

int main (int argc,char *argv[])
{
	
  static VString  filename = "";
  static VShort   first  = 0;
  static VShort   length = 0;
  static VShort   type   = 1;
  static VShort   minval = 0;
  static VShort   nproc  = 10;
  static VOptionDescRec  options[] = {
    {"mask",VStringRepn,1,(VPointer) &filename,VRequiredOpt,NULL,"mask"},
    {"minval",VShortRepn,1,(VPointer) &minval,VOptionalOpt,NULL,"Signal threshold"},
    {"first",VShortRepn,1,(VPointer) &first,VOptionalOpt,NULL,"First timestep to use"},
    {"length",VShortRepn,1,(VPointer) &length,VOptionalOpt,NULL,"Length of time series to use, '0' to use full length"},
    {"type",VShortRepn,1,(VPointer) &type,VOptionalOpt,TYPDict,"type of scaling to make correlations positive"},
    {"j",VShortRepn,1,(VPointer) &nproc,VOptionalOpt,NULL,"number of processors to use, '0' to use all"}
  };
  FILE *in_file,*out_file,*fp;
  VAttrList list=NULL,list1=NULL,out_list=NULL;
  VAttrListPosn posn;
  VImage mask=NULL;
  char *prg_name=GetLipsiaName("vecm");
  fprintf(stderr, "%s\n", prg_name);

  VParseFilterCmd (VNumber (options),options,argc,argv,&in_file,&out_file);
  if (type > 3) VError(" illegal type");


  /* read mask */
  fp = VOpenInputFile (filename, TRUE);
  list1 = VReadFile (fp, NULL);
  if (! list1) VError("Error reading mask file");
  fclose(fp);

  for (VFirstAttr (list1, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    VGetAttrValue (& posn, NULL,VImageRepn, & mask);
  }
  if (mask == NULL) VError(" no mask found");


  /* read functional data */
  if (! (list = VReadFile (in_file, NULL))) exit (1);
  fclose(in_file);



  /* omp-stuff */
#ifdef _OPENMP
  int num_procs=omp_get_num_procs();
  if (nproc > 0 && nproc < num_procs) num_procs = nproc;
  printf("using %d cores\n",(int)num_procs);
  omp_set_num_threads(num_procs);
#endif /* _OPENMP */


  /* main process */
  out_list = VECM(list,mask,minval,first,length,type );


  /* update geoinfo, 4D to 3D */
  VAttrList geolist = VGetGeoInfo(list);
  if (geolist != NULL) {
    double *D = VGetGeoDim(geolist,NULL);
    D[0] = 3;  /* 3D */
    D[4] = 1;  /* just one timestep */
    VSetGeoDim(geolist,D);
    VSetGeoInfo(geolist,out_list);
  }

  /* write output */
  VHistory(VNumber(options),options,prg_name,&list,&out_list);
  if (! VWriteFile (out_file, out_list)) exit (1);
  fprintf (stderr, "%s: done.\n", argv[0]);
  return 0;
}
