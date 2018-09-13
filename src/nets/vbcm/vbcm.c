/*
** BCM - bipartite connectivity mapping
**
** G.Lohmann, MPI-KYB, Aug 2017
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
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>

extern void gsl_sort_vector_index(gsl_permutation *,gsl_vector *);

#define SQR(x) ((x) * (x))
#define ABS(x) ((x) > 0 ? (x) : -(x))


#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/

extern gsl_matrix_float *DataMatrix(VImage *src,int nslices,VImage roi,int);
extern size_t NumVoxels(VImage src);
extern VImage VoxelMap(VImage roi);

typedef struct SpointStruct{
  VShort x;
  VShort y;
  VShort z;
} SPoint;

VDictEntry MetricDict[] = {
  { "pearson", 0 },
  { "spearman", 1 },
  { "MI", 2 },
  { NULL }
};


VImage VCallocImage(int nslices,int nrows,int ncols,VRepnKind repn,VImage ref)
{
  VImage image = VCreateImage(nslices,nrows,ncols,repn);
  if (image == NULL) VError(" error allocating image");
  VFillImage(image,VAllBands,0);
  if (ref != NULL) VCopyImageAttrs (ref,image);
  return image;
}

void ShowImage(gsl_matrix_float *A,char *filename)
{
  int i,j;
  int n=A->size1;
  int m=A->size2;
  VImage tstimage = VCallocImage(1,n,m,VFloatRepn,NULL);
  for (i=0; i <n; i++) {
    for (j=0; j<m; j++) {
      VPixel(tstimage,0,i,j,VFloat) = gsl_matrix_float_get(A,i,j);
    }
  }
  VAttrList out_list = VCreateAttrList();
  VAppendAttr(out_list,"image",NULL,VImageRepn,tstimage);
  FILE *out_file = fopen(filename,"w");
  VWriteFile (out_file, out_list);
}


VImage GetRoi(SPoint addr,int nslices,int nrows,int ncols,int radius)
{
  int i,b,r,c;
  int rad2 = radius*radius+1;

  int b0 = addr.z;
  int r0 = addr.y;
  int c0 = addr.x;

  i=0;
  for (b=b0-radius; b<=b0+radius; b++) {
    for (r=r0-radius; r<=r0+radius; r++) {
      for (c=c0-radius; c<=c0+radius; c++) {
	if (SQR(b-b0) + SQR(r-r0) + SQR(c-c0) > rad2) continue;
	i++;
      }
    }
  }
  VImage roi = VCreateImage(nslices,nrows,ncols,VBitRepn);
  VFillImage(roi,VAllBands,0);
  for (b=b0-radius; b<=b0+radius; b++) {
    if (b < 0 || b >= nslices) VError(" b %d",b);
    for (r=r0-radius; r<=r0+radius; r++) {
      if (r < 0 || r >= nrows) VError(" r %d",r);
      for (c=c0-radius; c<=c0+radius; c++) {
	if (c < 0 || c >= ncols) VError(" c %d",c);
	if (SQR(b-b0) + SQR(r-r0) + SQR(c-c0) > rad2) continue;
	VPixel(roi,b,r,c,VBit) = 1;
      }
    }
  }
  return roi;
}


double Pearson(const float *data1,const float *data2,int n)
{
  int i;
  double corr=0;
  for (i=0; i<n; i++) {   
    const double u = (double)data1[i];
    const double v = (double)data2[i];
    corr += u*v;
  }
  corr /= (double)n;
  return corr;
}

double Spearman(const float *data1,const float *data2,int n)
{
  int i;
  double nx = (double)n;
  double kx = nx*(nx*nx-1.0);
  double sxy=0.0;
  for (i=0; i<n; i++) {
    const double u = (double)data1[i];
    const double v = (double)data2[i];
    const double d = (u-v);
    sxy += d*d;
  }
  double rho = 1.0 - 6.0*sxy/kx;
  return rho;
}

/* if X,Y are normalized Gaussians, then I=-0.5*log(1-r^2), with r covariance */
double MutualInformation(const float *data1,const float *data2,int nt)
{
  int i;
  double covar,mi,u;
  double tiny=1.0e-10;

  covar = 0;
  for (i=0; i<nt; i++) {
    const double u = data1[i];
    const double v = data2[i];
    covar += u*v;
  }
  const double nx= nt;
  covar /= nx;

  mi = 0;
  u  = 1.0 - covar*covar;
  if (u > tiny) mi = -0.5 * log(u);
  if (mi < 0) {
    VWarning("MI: %f",mi);
    mi = 0;
  }
  return mi;
}


double Correlation(const float *data1,const float *data2,int n,int metric)
{
  double corr=0.0;
  if (metric == 0) corr = Pearson(data1,data2,n);
  if (metric == 1) corr = Spearman(data1,data2,n);
  if (metric == 2) corr = MutualInformation(data1,data2,n);
  if (gsl_isnan(corr) || gsl_isinf(corr)) corr = 0;
  return corr;
}


void VectorNormalize(gsl_vector_float *x)
{
  int i;
  double nx = (float)x->size;
  double sum1=0,sum2=0;
  for (i=0; i<x->size; i++) {
    sum1 += x->data[i];
    sum2 += x->data[i]*x->data[i];
  }
  double mean = sum1/nx;
  double sigma = sqrt((double)((sum2 - nx * mean * mean) / (nx - 1.0)));

  for (i=0; i<x->size; i++) {
    x->data[i] = (x->data[i] - mean)/sigma;
  }
}


/* bipartite eigenvector centrality mapping */
void VBiadjacencyECM(gsl_matrix_float *C,gsl_vector_float *xa,gsl_vector_float *xb)
{
  int na = C->size1;
  int nb = C->size2;
  int i,iter,maxiter=100;
  double d=0,sum=0,sum_old=0;

  /* ini */
  gsl_vector_float *ya = gsl_vector_float_calloc(na);
  gsl_vector_float *yb = gsl_vector_float_calloc(nb);  
  double nx = (double)(na+nb);
  for (i=0; i<na; i++) xa->data[i] = 1.0/nx;
  for (i=0; i<nb; i++) xb->data[i] = 1.0/nx;

  /* iterations */
  for (iter=0; iter<maxiter; iter++) {
    gsl_blas_sgemv (CblasNoTrans,1.0,C,xb,0.0,ya);
    gsl_blas_sgemv (CblasTrans,1.0,C,xa,0.0,yb);

    sum = 0;
    for (i=0; i<na; i++) sum += ya->data[i]*ya->data[i];
    for (i=0; i<nb; i++) sum += yb->data[i]*yb->data[i];
    sum = sqrt(sum);

    for (i=0; i<na; i++) xa->data[i] = ya->data[i]/sum;
    for (i=0; i<nb; i++) xb->data[i] = yb->data[i]/sum;

    d = fabs(sum-sum_old);
    if (iter > 0) fprintf(stderr," %5d   %f\n",(int)iter,d);
    if (d < 1.0e-8 && iter > 3) break;
    if (iter > 5 && sum > sum_old) break;
    sum_old = sum;
  }

  gsl_vector_float_free(ya);
  gsl_vector_float_free(yb);
}


/* bipartite degree centrality mapping */
void VBiadjacencyDegreeMap(gsl_matrix_float *C,gsl_vector_float *xa,gsl_vector_float *xb)
{
  size_t na = C->size1;
  size_t nb = C->size2;
  size_t i;

  gsl_vector_float *ya = gsl_vector_float_calloc(na);
  gsl_vector_float *yb = gsl_vector_float_calloc(nb);
  for (i=0; i<na; i++) xa->data[i] = 1.0/(double)na;
  for (i=0; i<nb; i++) xb->data[i] = 1.0/(double)nb;
  
  gsl_blas_sgemv (CblasNoTrans,1.0,C,xb,0.0,ya);
  gsl_blas_sgemv (CblasTrans,1.0,C,xa,0.0,yb);

  gsl_vector_float_memcpy(xa,ya);
  gsl_vector_float_memcpy(xb,yb);

  gsl_vector_float_free(ya);
  gsl_vector_float_free(yb);
}



/* project onto ROI */
void VNetworkProjection(gsl_matrix_float *C,gsl_vector_float *x)
{
  size_t na = C->size1;
  size_t nb = C->size2;
  size_t i,n1=0,n2=0;
 
  if (x->size == na) {
    n1 = na;
    n2 = nb;
  }
  else {
    n1 = nb;
    n2 = na;
  }
  float nx = (float)n2;

  fprintf(stderr," network projection\n");

  /* project bipartite graph onto ROI */
  gsl_matrix_float *A = gsl_matrix_float_calloc(n1,n1);
  if (!A) VError(" err allocating projection network");
  int progress=0;

#pragma omp parallel for schedule(guided) firstprivate(C,A)
  for (i=0; i<n1; i++) {
    size_t j=0,k=0;    
    if (i%10 == 0) fprintf(stderr," %d0  of  %lu\r",(int)(++progress),n1);
    for (j=0; j<n1; j++) {
      
      float w=0;
      for (k=0; k<n2; k++) {
	float u = gsl_matrix_float_get(C,i,k);
	float v = gsl_matrix_float_get(C,j,k);
	w += u*v;
      }
#pragma omp critical
      {
	gsl_matrix_float_set(A,i,j,w/nx);
      }
    }
  }
  fprintf(stderr,"\n");


  /* Hubs of projected network */
  float w=0;
  size_t j=0;
  for (i=0; i<n1; i++) {
    w=0;
    for (j=0; j<n1; j++) {
      w += gsl_matrix_float_get(A,i,j);
    }
    x->data[i] = w/nx;
  }
  gsl_matrix_float_free(A);
}



/* make matrix positive */
void VThresholdMatrix(gsl_matrix_float *C,float threshold)
{
  int i,j;
  size_t nneg=0,npos=0;
  float u=0,tiny=1.0e-6;
  double sum1=0,sum2=0;
  double nx = (double)(C->size1 * C->size2);
  for (i=0; i<C->size1; i++) {
    for (j=0; j<C->size2; j++) {
      u = gsl_matrix_float_get(C,i,j);
      sum1 += u;
      sum2 += u*u;
      if (u < threshold+tiny) nneg++;
      else npos++;
      if (u < threshold+tiny) gsl_matrix_float_set(C,i,j,0.0);
    }
  }
  if (npos < 1) VWarning(" Subthreshold correlations only");
  double mean = sum1/nx;
  double sigma = sqrt((double)((sum2 - nx * mean * mean) / (nx - 1.0)));
  fprintf(stderr," matrix mean,std: %f %f,  npos: %lu, nneg: %lu\n",mean,sigma,npos,nneg);
}


int main (int argc,char *argv[])
{ 
  static SPoint addr1;
  static SPoint addr2;  
  static VString roi1_filename = "";
  static VString roi2_filename = ""; 
  static VShort type=1;
  static VFloat threshold = 0.0;
  static VShort radius=5;
  static VShort   metric = 0;
  static VBoolean plotboth=FALSE;
  static VBoolean normalize=TRUE;
  static VOptionDescRec  options[] = {   
    {"roi1",VStringRepn,1,(VPointer) &roi1_filename,VOptionalOpt,NULL,"ROI 1"},
    {"roi2",VStringRepn,1,(VPointer) &roi2_filename,VOptionalOpt,NULL,"ROI 2"},
    {"seed1",VShortRepn,3,(VPointer) &addr1,VOptionalOpt,NULL,"Voxel address of seed point (x,y,z)"},
    {"seed2",VShortRepn,3,(VPointer) &addr2,VOptionalOpt,NULL,"Voxel address of seed point (x,y,z)"},
    {"type", VShortRepn,1,(VPointer) &type,VOptionalOpt,NULL,"Type, 0:ECM, 1: DCM, 2: project onto ROI1, 3: project onto ROI2"},    
    {"threshold", VFloatRepn,1,(VPointer) &threshold,VOptionalOpt,NULL,"Matrix threshold"},
    {"metric",VShortRepn,1,(VPointer) &metric,VOptionalOpt,MetricDict,"Correlation metric"},
    {"radius", VShortRepn,1,(VPointer) &radius,VOptionalOpt,NULL,"Radius around seed voxel"},
    {"both", VBooleanRepn,1,(VPointer) &plotboth,VOptionalOpt,NULL,"Whether to plot results of both ROIs"},
    {"normalize", VBooleanRepn,1,(VPointer) &normalize,VOptionalOpt,NULL,"Whether to normalize results"},
  };
  FILE *out_file=NULL;
  VString in_file=NULL;
  VImage roi1=NULL,roi2=NULL;
  int i;
  char *prg = GetLipsiaName("vbcm");


  /* Parse command line arguments and identify files: */
  VParseFilterCmdX (VNumber (options), options, argc, argv,&in_file,&out_file);
  fprintf(stderr," type= %d\n",type);

  /* read functional data */
  VAttrList list = VReadAttrList(in_file,0L,TRUE,FALSE);
  if (list == NULL) VError(" error reading input file %s",in_file);

  /* get pointers to image data */
  int nrows=0,ncols=0,nt=0;
  int nslices = VAttrListNumImages(list);
  VImage *src = VAttrListGetImages(list,nslices);
  VImageDimensions(src,nslices,&nt,&nrows,&ncols);
  fprintf(stderr," image dims:  %d  %d  %d,  nt: %d\n",nslices,nrows,ncols,nt);


  /* omp-stuff */
#ifdef _OPENMP
  int num_procs=omp_get_num_procs();
  printf("using %d cores\n",(int)num_procs);
  omp_set_num_threads(num_procs);
#endif /* _OPENMP */


  /* ROI 1 */
  if (strlen(roi1_filename) > 0) {
    VAttrList list1 = VReadAttrList(roi1_filename,0L,TRUE,FALSE);
    if (list1 == NULL) VError(" error reading roi file %s",roi1_filename);
    roi1 = VReadImage(list1);
    if (roi1 == NULL) VError(" err reading %s",roi1_filename);
  }

  /* ROI 2 */
  if (strlen(roi2_filename) > 0) {
    VAttrList list2 = VReadAttrList(roi2_filename,0L,TRUE,FALSE);
    if (list2 == NULL) VError(" error reading roi file %s",roi2_filename);
    roi2 = VReadImage(list2);
    if (roi2 == NULL) VError(" err reading %s",roi2_filename);
  }


  /* ROIs from seed voxels */
  if (roi1 == NULL) roi1 = GetRoi(addr1,nslices,nrows,ncols,radius);
  if (roi2 == NULL) roi2 = GetRoi(addr2,nslices,nrows,ncols,radius);


  /* voxel maps */
  VImage map1  = VoxelMap(roi1);
  VImage map2  = VoxelMap(roi2);
  size_t nvox1 = NumVoxels(roi1);
  size_t nvox2 = NumVoxels(roi2);
  fprintf(stderr," nvox1: %lu,  nvox2: %lu\n",nvox1,nvox2);


  /* read data matrices */
  gsl_matrix_float *X1 = DataMatrix(src,nslices,roi1,(int)metric);
  gsl_matrix_float *X2 = DataMatrix(src,nslices,roi2,(int)metric);
  

  /* correlation matrix */
  fprintf(stderr," CorrMatrix...\n");
  int rad2 = 3*3;     /* exclusion radius */
  gsl_matrix_float *A = gsl_matrix_float_calloc(nvox1,nvox2);
  if (A == NULL) VError(" err allocating corr matrix, %d x %d",nvox1,nvox2);
  int progress=0;

#pragma omp parallel for schedule(guided) firstprivate(A)
  for (i=0; i<nvox1; i++) {
    if (i%100 == 0) fprintf(stderr," %d00\r",(int)(++progress));
    int bi = VPixel(map1,0,0,i,VShort);
    int ri = VPixel(map1,0,1,i,VShort);
    int ci = VPixel(map1,0,2,i,VShort);

    int j=0;
    float u=0;
    float *tmp = (float *) VCalloc(nvox2,sizeof(float));
    for (j=0; j<nvox2; j++) tmp[j] = -1.0;
    const float *data1 = gsl_matrix_float_const_ptr(X1,i,0);
    for (j=0; j<nvox2; j++) {
      int bj = VPixel(map2,0,0,j,VShort);
      int rj = VPixel(map2,0,1,j,VShort);
      int cj = VPixel(map2,0,2,j,VShort);
      int di = SQR(bi-bj) + SQR(ri-rj) + SQR(ci-cj);
      if (di < rad2) continue;

      const float *data2 = gsl_matrix_float_const_ptr(X2,j,0);
      tmp[j] = Correlation(data1,data2,nt,(int)metric);
    }

#pragma omp critical
    {
      for (j=0; j<nvox2; j++) {
	u = tmp[j];
	gsl_matrix_float_set(A,i,j,u);
      }
    }
    VFree(tmp);
  }
  fprintf(stderr,"\n");

  /*
  ShowImage(A,"test.v");
  exit(0);
  */

  /* Biadjacency ECM */
  fprintf(stderr," Biadjacency Connectivity...\n");
  gsl_vector_float *xa = gsl_vector_float_calloc(nvox1);
  gsl_vector_float *xb = gsl_vector_float_calloc(nvox2);

  fprintf(stderr," type= %d\n",type);


  VThresholdMatrix(A,threshold);

  switch(type) {
  case 0:
    VBiadjacencyECM(A,xa,xb);
    break;
  case 1:
    VBiadjacencyDegreeMap(A,xa,xb);
    break;
  case 2:
    VNetworkProjection(A,xa);
    break;
  case 3:
    VNetworkProjection(A,xb);
    break;
  default:
    VError(" unknown type");
  }



  /* normalize */
  if (normalize) {
    VectorNormalize(xa);
    VectorNormalize(xb);
  }

  /* output */
  VImage dest = VCallocImage(nslices,nrows,ncols,VFloatRepn,src[0]);
  int b,r,c;
  for (i=0; i<nvox1; i++) {
    b = VPixel(map1,0,0,i,VShort);
    r = VPixel(map1,0,1,i,VShort);
    c = VPixel(map1,0,2,i,VShort);
    VPixel(dest,b,r,c,VFloat) = 100.0*xa->data[i];
  }
  if (plotboth) {
    for (i=0; i<nvox2; i++) {
      b = VPixel(map2,0,0,i,VShort);
      r = VPixel(map2,0,1,i,VShort);
      c = VPixel(map2,0,2,i,VShort);
      VPixel(dest,b,r,c,VFloat) = 100.0*xb->data[i];
    }
  }
  

  /* write output to disk */
  VAttrList out_list = VCreateAttrList();
  VAttrList geolist = VGetGeoInfo(list);
  if (geolist != NULL) {
    double *D = VGetGeoDim(geolist,NULL);
    D[0] = 3;  /* 3D */
    D[4] = 1;
    VSetGeoDim(geolist,D);
    VSetGeoInfo(geolist,out_list);
  }
  VHistory(VNumber(options),options,prg,&list,&out_list);
  VAppendAttr(out_list,"image",NULL,VImageRepn,dest);
  VWriteFile (out_file, out_list);
  fclose(out_file);
  fprintf (stderr, "%s: done.               \n", argv[0]);
  exit(0);
}
