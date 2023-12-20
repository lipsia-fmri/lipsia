/*
** read data matrix, map of voxel addresses
**
** G.Lohmann, Jan 2020
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>

#include "viaio/Vlib.h"
#include "viaio/VImage.h"
#include "viaio/option.h"

extern void ZNorm(gsl_vector *);
extern void VCheckImage(VImage src);

void YCheckImage(VImage src)
{
  VAttrList out_list = VCreateAttrList();
  VAppendAttr(out_list,"image",NULL,VImageRepn,src);
  FILE *fp = VOpenOutputFile ("test.v", TRUE);
  if (! VWriteFile (fp, out_list)) exit (1);
  fclose(fp);
  exit(0);
}

int test_ascii(int val)
{
  if (val >= 'a' && val <= 'z') return 1;
  if (val >= 'A' && val <= 'Z') return 1;
  if (val >= '0' && val <= '9') return 1;
  if (val ==  '-') return 1;
  if (val ==  ' ') return 1;
  if (val == '\0') return 1;
  if (val == '\n') return 1;
  if (val == '\r') return 1;
  if (val == '\t') return 1;
  if (val == '\v') return 1;
  return 0;
}

int line_empty(char *buf,int len)
{
  int i;
  for (i=0; i<len; i++) {
    if (buf[i] != ' ' && buf[i] != '\n' && buf[i] != 0) return 0;
  }
  return 1;
}

int CheckBuffer(char *buf,int len)
{
  int j;
  if(strlen(buf) < 1) return 0;
  if (buf[0] == '%' || buf[0] == '#' || buf[0] == '/' || buf[0] == '\n') return 0;
  for (j=0; j<len; j++) {
    if (buf[0] == '\t') buf[j] = ' ';
  }
  if (line_empty(buf,len) > 0) return 0;
  return 1;
}




VImage VoxelMap(VImage mask,size_t *nvoxels)
{
  int b,r,c;
  int nslices = VImageNBands(mask);
  int nrows = VImageNRows(mask);
  int ncols = VImageNColumns(mask);

  /* count number of non-zero voxels */
  size_t nvox = 0;  
  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	float u = VGetPixel(mask,b,r,c);
	if (u < TINY) continue;
	nvox++;
      }
    }
  }
  (*nvoxels) = nvox;

  /* voxel addresses */
  VImage map = VCreateImage(1,4,nvox,VIntegerRepn);
  if (map == NULL) VError(" error allocating addr map");
  VFillImage(map,VAllBands,0);
  VImage xmap=NULL;
  VExtractAttr (VImageAttrList(mask),"map",NULL,VImageRepn,&xmap,FALSE);
  VCopyImageAttrs (mask,map);
  VSetAttr(VImageAttrList(map),"nvoxels",NULL,VLongRepn,(VLong)nvox);
  VSetAttr(VImageAttrList(map),"nslices",NULL,VLongRepn,(VLong)nslices);
  VSetAttr(VImageAttrList(map),"nrows",NULL,VLongRepn,(VLong)nrows);
  VSetAttr(VImageAttrList(map),"ncols",NULL,VLongRepn,(VLong)ncols);

  VPixel(map,0,3,0,VInteger) = nslices;
  VPixel(map,0,3,1,VInteger) = nrows;
  VPixel(map,0,3,2,VInteger) = ncols;

  size_t i = 0;
  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	float u = VGetPixel(mask,b,r,c);
	if (u < TINY) continue;
	VPixel(map,0,0,i,VInteger) = b;
	VPixel(map,0,1,i,VInteger) = r;
	VPixel(map,0,2,i,VInteger) = c;
	i++;
      }
    }
  }
  return map;
}




gsl_vector *ReadRegressor(VString filename)
{
  double x=0;
  int n,len=4096;
  char *buf = (char *)VCalloc(len,sizeof(char));

  FILE *fp = fopen(filename,"r");
  if (!fp) VError(" Error opening file %s",filename);

  n=0;
  while (!feof(fp)) {
    memset(buf, 0, len);
    if (fgets(buf,len,fp) == NULL) continue;
    if (CheckBuffer(buf,len) < 1) continue;
    if (! test_ascii((int)buf[0])) VError(" Input file must be a text file");
    n++;
  }
  gsl_vector *vec = gsl_vector_calloc(n);
  rewind(fp);

  n=0;
  while (!feof(fp)) {
    memset(buf, 0, len);
    if (fgets(buf,len,fp) == NULL) continue;
    if (CheckBuffer(buf,len) < 1) continue;
    sscanf(buf,"%lf",&x);
    vec->data[n] = x;
    n++;
  }
  fclose(fp);

  return vec;
}


/* get nedges */
size_t GetNumEdges(VArgVector in_files)
{
  VAttrListPosn posn;
  VString str;
  VLong nedges=0;

  VString in_filename = ((VString *) in_files.vector)[0];
  VAttrList list = VReadAttrList(in_filename,0L,TRUE,FALSE);
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    str = VGetAttrName(&posn);
    if (strcmp(str,"nedges") == 0) {
      VGetAttrValue (& posn, NULL, VLongRepn, & nedges);
      break;
    }
  }
  VDestroyAttrList(list);
  if (nedges < 1) VError(" nedges");
  return (size_t)nedges;
}


/* get number of components */
size_t GetNumComponents(VArgVector in_files)
{
  VAttrListPosn posn;
  VString str;
  VLong n=0;

  VString in_filename = ((VString *) in_files.vector)[0];
  VAttrList list = VReadAttrList(in_filename,0L,TRUE,FALSE);
  for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    str = VGetAttrName(&posn);
    if (strcmp(str,"ncomponents") == 0) {
      VGetAttrValue (& posn, NULL, VLongRepn, & n);
      break;
    }
  }
  VDestroyAttrList(list);
  if (n < 1) VWarning(" ncomponents");
  return (size_t)n;
}


/* sign flip for ricci curvature */
void ScaleRicci(float *A,size_t nedges)
{
  size_t i;
  float u,nx = sqrt((double)nedges);
  for (i=0; i<nedges; i++) {
    u = A[i]/nx;
    A[i] = -u;
  }
}


void VReadBundles(VArgVector in_files,size_t nedges,size_t numsubjects,
		  gsl_matrix_float *X0,gsl_matrix_float *X1,size_t *n0,size_t *n1)
{
  size_t i,j,l=0,m=0;
  VAttrListPosn posn;
  VString str;
  VBundle zbundle=NULL;
  float *A;

  size_t numfiles = numsubjects;
  if (numfiles != (size_t)in_files.number) VError(" numfiles: %lu %lu",numfiles,in_files.number);

  size_t mist=0,total=0;
  size_t k=0;
  for (i=0; i<numsubjects; i++) {
    VString in_filename = ((VString *) in_files.vector)[k];
    VAttrList list = VReadAttrList(in_filename,0L,TRUE,FALSE);
    k++;

    zbundle = NULL;
    for (VFirstAttr (list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
      if (VGetAttrRepn (& posn) != VBundleRepn) continue;
      VGetAttrValue (& posn, NULL, VBundleRepn, & zbundle);
      A = (float *)zbundle->data;

           
      str = VGetAttrName(&posn);
      if (X0 != NULL && strcmp(str,"Correlation") == 0) {
	m++;
	for (j=0; j<nedges; j++) {
	  if (gsl_finite(A[j])==0) { A[j] = 0;  mist++; }
	  total++;
	  gsl_matrix_float_set(X0,i,j,A[j]);
	}	
      }
      if (X1 != NULL && strcmp(str,"Ricci") == 0) {
	l++;
	ScaleRicci(A,nedges);
	for (j=0; j<nedges; j++) {
	  gsl_matrix_float_set(X1,i,j,A[j]);
	}	
      }      
    }
    VDestroyAttrList(list);
  }
  if (mist > 0) VWarning(" read bundles, mist: %lu  of %lu,  %f",mist,total,(double)mist/(double)total);
  if (k != numfiles) VError(" k: %lu %lu",k,numfiles);
  *n0 = m;
  *n1 = l;
}


/* read selection file */
VImage ReadSelectionMap(VAttrList sel_list)
{
  VAttrListPosn posn;
  VImage map;
  VString str=NULL;
  
  for (VFirstAttr (sel_list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    str = VGetAttrName(&posn);
    if (strcmp(str,"map") == 0) {
      VGetAttrValue (& posn, NULL, VImageRepn, & map);
    }
  }
  if (!map) VError(" map not found");
  return map;
}


/* read selection file */
long *ReadSelectionI(VAttrList sel_list)
{
  VAttrListPosn posn;
  VBundle rowindex=NULL;
  VString str=NULL;
  
  for (VFirstAttr (sel_list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    str = VGetAttrName(&posn);
    if (strcmp(str,"RowIndex") == 0) {
      VGetAttrValue (& posn, NULL, VBundleRepn, & rowindex);
    }
  }
  long *I = (long *)rowindex->data;
  return I;
}



/* read selection file */
long *ReadSelectionJ(VAttrList sel_list)
{
  VAttrListPosn posn;
  VBundle colindex=NULL;
  VString str=NULL;
  
  for (VFirstAttr (sel_list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    str = VGetAttrName(&posn);
    if (strcmp(str,"ColIndex") == 0) {
      VGetAttrValue (& posn, NULL, VBundleRepn, & colindex);
    }
  }
  long *J = (long *)colindex->data;
  return J;
}


/* read selection file */
VImage ReadSelectionImage(VAttrList sel_list)
{
  VAttrListPosn posn;
  VImage image=NULL;
  VString str=NULL;
  
  for (VFirstAttr (sel_list, & posn); VAttrExists (& posn); VNextAttr (& posn)) {
    if (VGetAttrRepn (& posn) != VImageRepn) continue;
    str = VGetAttrName(&posn);
    if (strcmp(str,"image") == 0) {
      VGetAttrValue (& posn, NULL, VImageRepn, &image);
    }
  }
  return image;
}



/* read input confound files */
gsl_matrix *ReadConfounds(VArgVector confound_files,int numsubjects,FILE *fp,int verbose)
{
  int i,j;
  int numconfounds = (size_t)confound_files.number;
  if (numconfounds == 0) return NULL;
  
  VString in_filename=NULL;
  int dim=0;
  for (j=0; j<numconfounds; j++) {
    in_filename = ((VString *) confound_files.vector)[j];
    if (strcmp(in_filename,"xxx") == 0) continue;
    if (access(in_filename, F_OK) != 0) continue;
    dim++;
  }

  gsl_vector *tmp = NULL;
  gsl_matrix *X = gsl_matrix_calloc(numsubjects,dim);
  gsl_matrix_set_all(X,1.0);

  dim=0;
  for (j=0; j<numconfounds; j++) {
    in_filename = ((VString *) confound_files.vector)[j];
    if (strcmp(in_filename,"xxx") == 0) continue;
    if (access(in_filename, F_OK) != 0) continue;
    tmp = ReadRegressor(in_filename);
    if (tmp==NULL) continue;
    /* ZNorm(tmp); */
    
    for (i=0; i<numsubjects; i++) {
      gsl_matrix_set(X,i,dim,tmp->data[i]);
    }
    
    char *str = strchr(in_filename,'_');
    if (fp != NULL && verbose > 0) fprintf(fp,"# %s\n",str);
    if (verbose) fprintf(stderr,"# %3d  %s\n",dim,str);
    dim++;
  }
  gsl_vector_free(tmp);
  if (fp != NULL && verbose > 1) fprintf(fp,"#\n");
  if (verbose) fprintf(stderr,"\n");
  if (dim < 1) return NULL;
  
  return X;
}


void RowNormalizeBundle(gsl_matrix_float *X)
{
  if (X==NULL) return;
  size_t i,j;
  double u,z,s1=0,s2=0,mean=0,var=0,sd=0;
  double nx = (double)X->size2;

  for (i=0; i<X->size1; i++) {
    s1 = s2 = nx = 0;
    for (j=0; j<X->size2; j++) {
      u = (double)gsl_matrix_float_get(X,i,j);
      s1 += u;
      s2 += u*u;
      nx++;
    }
    mean = s1/nx;
    var = (s2 - nx * mean * mean) / (nx - 1.0);
    sd = sqrt(var);
    if (sd < TINY) VError(" sd=0, RowNormalizeBundle");

    for (j=0; j<X->size2; j++) {
      u = (double)gsl_matrix_float_get(X,i,j);
      z = (u-mean)/sd;
      gsl_matrix_float_set(X,i,j,(float)z);
    }
  }
}


void RowRankNormBundle(gsl_matrix_float *X)
{
  if (X==NULL) return;
  size_t i;
  size_t n = X->size2;
  double nx = (double)n;
  
#pragma omp parallel for schedule(guided)
  for (i=0; i<X->size1; i++) {
    size_t j;
    size_t *perm = (size_t *) VCalloc(n,sizeof(size_t));
    double *wx = (double *)VCalloc(n,sizeof(double));
    double *x = (double *)VCalloc(n,sizeof(double));
    for (j=0; j<n; j++) x[j] = (double)gsl_matrix_float_get(X,i,j);
    
    gsl_sort_index (perm,x,1,n);
    for (j=0; j<n; j++) wx[perm[j]] = (double)j/nx;
    for (j=0; j<n; j++) x[j] = 2.0*(wx[j]-0.5);

    for (j=0; j<X->size2; j++) {
      gsl_matrix_float_set(X,i,j,(float)x[j]);
    }
    VFree(x);
    VFree(wx);
    VFree(perm);
  }
}


VImage ReadMask(VString mask_filename,int verbose)
{
  VAttrList listm=NULL;
  VImage mask=NULL;
  size_t i,n=0;
  
  if (strlen(mask_filename) > 1 && strcmp(mask_filename,"xxx") != 0) {
    listm = VReadAttrList(mask_filename,0L,TRUE,FALSE);
    if (listm == NULL) VError(" error reading mask file %s",mask_filename);
    mask = VReadImage(listm);
    if (mask == NULL) VError(" err reading %s",mask_filename);
    if (VPixelRepn(mask) != VBitRepn && VPixelRepn(mask) != VShortRepn)
      VError(" mask must be in bit or short repn");

    if (VPixelRepn(mask) == VBitRepn) {
      VBit *pp = VImageData(mask);
      for (i=0; i<VImageNPixels(mask); i++) if (pp[i] > 0) n++;
    }
    if (VPixelRepn(mask) == VShortRepn) {
      VShort *pp = VImageData(mask);
      for (i=0; i<VImageNPixels(mask); i++) if (pp[i] > 0) n++;
    }    
    if (verbose) fprintf(stderr,"# mask: %s,  nvox: %lu\n",mask_filename,n);
  }
  return mask;
}


size_t ApplyMask(VImage mask,VImage map,long *I,long *J,char *table,size_t nedges,int type)
{
  size_t i,j,k,n=0;
  double u,v;

  if (map == NULL || mask == NULL) return nedges;

  int nslices = VPixel(map,0,3,0,VInteger);
  int nrows = VPixel(map,0,3,1,VInteger);
  int ncols = VPixel(map,0,3,2,VInteger);
  if (nslices != (int)VImageNBands(mask)) VError(" mask, inconsistent dims, slices");
  if (nrows != (int)VImageNRows(mask)) VError(" mask, inconsistent dims, rows");
  if (ncols != (int)VImageNColumns(mask)) VError(" mask, inconsistent dims, columns");
  
  for (k=0; k<nedges; k++) {
    i = I[k];
    j = J[k];
    int bi = (int)VPixel(map,0,0,i,VInteger);
    int ri = (int)VPixel(map,0,1,i,VInteger);
    int ci = (int)VPixel(map,0,2,i,VInteger);
    
    int bj = (int)VPixel(map,0,0,j,VInteger);
    int rj = (int)VPixel(map,0,1,j,VInteger);
    int cj = (int)VPixel(map,0,2,j,VInteger);
  
    u = VGetPixel(mask,bi,ri,ci);
    v = VGetPixel(mask,bj,rj,cj);

    if (type == 0) {  /* or */
      table[k] = 0;
      if (fabs(u) > 0.01 || fabs(v) > 0.01) table[k] = 1;
    }
    if (type == 1) {  /* and */
      table[k] = 0;
      if (fabs(u) > 0.01 && fabs(v) > 0.01) table[k] = 1;
    }
    if (type == 2) {   /* absolutely not */
      table[k] = 1;
      if (fabs(u) > 0.01 || fabs(v) > 0.01) table[k] = 0;
    }
    if (type == 3) {   /* not */
      table[k] = 1;
      if (fabs(u) > 0.01 && fabs(v) > 0.01) table[k] = 0;
    }
  }
  n=0;
  for (k=0; k<nedges; k++) {
    if (table[k] > 0) n++;
  }
  fprintf(stderr,"# mask,  non-zero: %lu of %lu,  type: %d\n",n,nedges,type);
  return n;
}


void ColumnNormalizeBundle(gsl_matrix_float *X)
{
  size_t i,j;
  double u,z,s1=0,s2=0,mean=0,var=0,sd=0;
  double nx = (double)X->size1;

  for (i=0; i<X->size2; i++) {
    s1 = s2 = 0;
    for (j=0; j<X->size1; j++) {
      u = (double)gsl_matrix_float_get(X,j,i);
      s1 += u;
      s2 += u*u;
    }
    mean = s1/nx;
    var = (s2 - nx * mean * mean) / (nx - 1.0);
    sd = sqrt(var);

    for (j=0; j<X->size1; j++) {
      u = (double)gsl_matrix_float_get(X,j,i);
      z = (u-mean)/sd;
      gsl_matrix_float_set(X,j,i,(float)z);
    }
  }
}


/* double centering */
void KNormBundle(gsl_matrix_float *X)
{
  size_t i,j;
  double u,v,s,sij,ni=(double)X->size1,nj=(double)X->size2;
  gsl_matrix_float *A = gsl_matrix_float_calloc(X->size1,X->size2);
  gsl_matrix_float_memcpy(A,X);


  sij=0;
  for (i=0; i<X->size1; i++) {
    for (j=0; j<X->size2; j++) {
      sij += gsl_matrix_float_get(X,i,j);
    }
  }
  sij /= (ni*nj);

  double *tmpj = (double *) VCalloc(X->size2,sizeof(double));
  for (j=0; j<X->size2; j++) {
    s=0;
    for (i=0; i<X->size1; i++) s += gsl_matrix_float_get(X,i,j);
    tmpj[j] = s/ni;
  }

  double *tmpi = (double *) VCalloc(X->size1,sizeof(double));
  for (i=0; i<X->size1; i++) {
    s=0;
    for (j=0; j<X->size2; j++) s += gsl_matrix_float_get(X,i,j);
    tmpi[i] = s/nj;
  }
  
  for (i=0; i<X->size1; i++) {
    for (j=0; j<X->size2; j++) {
      u = gsl_matrix_float_get(X,i,j);
      v = u - tmpi[i] - tmpj[j] + sij;
      gsl_matrix_float_set(A,i,j,v);
    }
  }

  gsl_matrix_float_memcpy(X,A);
  gsl_matrix_float_free(A);
  VFree(tmpi);
  VFree(tmpj);
}


void FisherZBundle(gsl_matrix_float *X)
{
  size_t i,n;
  double u;

  n = X->size1*X->size2;
  float *pp = gsl_matrix_float_ptr(X,0,0);
  for (i=0; i<n; i++) {
    u = (double)pp[i];
    pp[i] = (float)gsl_atanh(u);
  }
}
