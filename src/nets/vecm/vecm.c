/*
** ECM - eigenvector centrality mapping
**
** G.Lohmann, MPI-KYB, updated Oct 2018
*/

#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/


#define SQR(x) ((x) * (x))
#define ReLU(x) ((x) > 0 ? (x) : 0)

extern void   VMatrixECM(gsl_matrix_float *,float *,int,int,VImage);
extern void   VMatrixProjection(gsl_matrix_float *X,float *ev,int,int seed);
extern VImage VoxelMap(VAttrList list,VImage mask,size_t nvox);
extern size_t NumVoxels(VImage);
extern gsl_matrix_float *ReadDataECM(VAttrList,VImage,VImage,VShort,VShort,int,size_t);
extern VImage WriteOutput(VAttrList,VImage,float *,size_t);
extern int Irreducible(gsl_matrix_float *X,int type);


void NormVec(float *x,size_t nvox)
{
 size_t i;
 float sum=0;
 sum = 0;
 for (i=0; i<nvox; i++) sum += (x[i]*x[i]);
 sum = sqrt(sum);
 for (i=0; i<nvox; i++) x[i] /= sum;
}

/* needed for ECM-RLC, negative */
void SignFlip(gsl_matrix_float *X)
{
  float u=0;
  size_t nt = X->size2;
  size_t i,j,n=nt/2;
  for (i=0; i<X->size1; i++) {
    for (j=n; j<nt; j++) {
      u = gsl_matrix_float_get(X,i,j);
      gsl_matrix_float_set(X,i,j,-u);
    }
  }
}


/* fast, memory efficient power iteration, used in ECM-RLC, ECM-Add1 */
void PowerIteration(gsl_matrix_float *X,float *ev,int type,int maxiter)
{
  int iter;
  float vsum=0,d=0;
  size_t i;
  size_t nvox = X->size1;
  size_t nt = X->size2;
  float scale = 1.0/(float)nt;
  
  fprintf(stderr," power iteration...\n");
  gsl_vector_float *y = gsl_vector_float_calloc(nvox);
  gsl_vector_float *z = gsl_vector_float_calloc(nt);

  Irreducible(X,type);  /* check if corr matrix is irreducible */
  
  NormVec(ev,nvox);
  for (i=0; i<nvox; i++) y->data[i] = ev[i];


  /* power iterations */
  for (iter=0; iter<maxiter; iter++) {

    /* R = (X^t X)/n */
    gsl_blas_sgemv (CblasTrans,1.0,X,y,0.0,z);
    if (type == 7) SignFlip(X);
    gsl_blas_sgemv (CblasNoTrans,1.0,X,z,0.0,y);
    if (type == 7) SignFlip(X);
    gsl_vector_float_scale(y,scale);

    /* ECM-Add1,  (R+E)v = Rv + Ev */
    if (type == 1) {
      vsum = 0;
      for (i=0; i<nvox; i++) vsum += ev[i];
      for (i=0; i<nvox; i++) y->data[i] += vsum;
    }    

    /* normalize */
    NormVec(y->data,nvox);
    
    /* check convergence */
    d = 0;
    for (i=0; i<nvox; i++) {
      d += SQR(ev[i] - y->data[i]);
      ev[i] = y->data[i];
    }
    fprintf(stderr," %5d   %.6f\n",(int)iter,d);

    if (iter > 2 && d < 1.0e-6) break;
  }
  gsl_vector_float_free(y);
  gsl_vector_float_free(z);
}



float *VECM(gsl_matrix_float *X,int type,VBoolean project,int seed,int maxiter,VImage map)
{
  size_t i,j;
  size_t nvox = X->size1;
  size_t nt = X->size2;
  if (type == 0 || type == 7) nt /= 2;

  
  /* ini eigenvector */
  float u=0;
  float *eigvec = (float *) VCalloc(nvox,sizeof(float));
  for (i=0; i<nvox; i++) eigvec[i] = 1.0/(float)nt;


  /* no matrix projection */
  if (!project) {

    switch (type) {
      
    case 0:   /*  ECM-RLC positive */
      for (i=0; i<X->size1; i++) {
	for (j=0; j<nt; j++) {
	  u = gsl_matrix_float_get(X,i,j);
	  gsl_matrix_float_set(X,i,j+nt,fabs(u));
	}
      }
      PowerIteration(X,eigvec,type,maxiter);
      break;


    case 1:   /* ECM-Add1 */
      PowerIteration(X,eigvec,type,maxiter);
      break;

      
    case 6:   /* no filter, risky */
      PowerIteration(X,eigvec,1,maxiter);  /* initialize using ecm-add1 */
      PowerIteration(X,eigvec,0,maxiter);
      break;

      /* ECM-RLC negative */
    case 7:    
      for (i=0; i<X->size1; i++) {
	for (j=0; j<nt; j++) {
	  u = gsl_matrix_float_get(X,i,j);
	  gsl_matrix_float_set(X,i,j,fabs(u));
	  gsl_matrix_float_set(X,i,j+nt,u);
	}
      }
      PowerIteration(X,eigvec,type,maxiter);
      break;

      
    default:  /* dense matrix representation */
      VMatrixECM(X,eigvec,type,maxiter,map);
      break;
      
    }
  }

  /* Matrix projection, ECM-project */
  else {
    VMatrixProjection(X,eigvec,type,seed);
  }

  /* rescale for better visualization */
  float nx = (float)nvox;
  nx = sqrt(nx);
  for (i=0; i<nvox; i++) eigvec[i] *= nx;
 
  return eigvec;
}



VDictEntry TYPDict[] = {
  { "rlc",     0, 0,0,0,0 },
  { "add",     1, 0,0,0,0 },
  { "pos",     2, 0,0,0,0 },
  { "abs",     3, 0,0,0,0 },
  { "neg",     4, 0,0,0,0 },
  { "gauss",   5, 0,0,0,0 },
  { NULL, 0,0,0,0,0 }
};


int main (int argc,char *argv[])
{	
  static VString  mask_filename = "";
  static VShort   first  = 0;
  static VShort   length = 0;
  static VShort   type   = 0;
  static VShort   maxiter= 20;
  static VBoolean project = FALSE;
  static VInteger  seed   = 99402622;
  static VShort   nproc  = 0;
  static VOptionDescRec  options[] = {
    {"mask",VStringRepn,1,(VPointer) &mask_filename,VRequiredOpt,NULL,"Region of interest mask"},
    {"first",VShortRepn,1,(VPointer) &first,VOptionalOpt,NULL,"First timestep to use"},
    {"length",VShortRepn,1,(VPointer) &length,VOptionalOpt,NULL,"Length of time series to use, '0' to use full length"}, 
    {"metric",VShortRepn,1,(VPointer) &type,VOptionalOpt,TYPDict,"Type of correlation metric"},
    {"iterations",VShortRepn,1,(VPointer) &maxiter,VOptionalOpt,NULL,"Maximum number of power iterations"},
    {"project",VBooleanRepn,1,(VPointer) &project,VOptionalOpt,NULL,"Whether to do matrix projection"},
    {"seed",VIntegerRepn,1,(VPointer) &seed,VOptionalOpt,NULL,"Seed for random number generation (only for ECM-project)"},   
    {"j",VShortRepn,1,(VPointer) &nproc,VOptionalOpt,NULL,"Number of processors to use, '0' to use all"}
  };
  FILE *in_file=NULL,*out_file=NULL;
  VString in_filename=NULL;
  char *prg_name=GetLipsiaName("vecm");
  fprintf(stderr, "%s\n", prg_name);

  
  /* parse command line */
  VParseFilterCmdZ (VNumber (options),options,argc,argv,&in_file,&out_file,&in_filename);
  if (type < 0 || type > 8) VError(" unknown metric");
  if (type < 6) fprintf(stderr," Correlation metric %d:  %s\n",type,TYPDict[type].keyword);
  if (project && (type < 2 || type == 7))
    VError(" matrix projection not needed for metric %d (%s)\n",type,TYPDict[type].keyword);
  
  
  /* omp-stuff */
#ifdef _OPENMP
  int num_procs=omp_get_num_procs();
  if (nproc > 0 && nproc < num_procs) num_procs = nproc;
  if (type > 3 || project == TRUE) fprintf(stderr,"using %d cores\n",(int)num_procs);
  omp_set_num_threads(num_procs);
#endif /* _OPENMP */


  
  /* read functional data */
  VAttrList list = VReadAttrListZ(in_file,in_filename,0L,TRUE,FALSE);
  if (list == NULL) VError(" error reading input file %s",in_file);
  VAttrList geolist = VGetGeoInfo(list);


  /* read mask */
  VAttrList listm = VReadAttrList(mask_filename,0L,TRUE,FALSE);
  VImage mask = VReadImage(listm);
  if (mask == NULL) VError(" no mask found");
  size_t nvox = NumVoxels(mask);

  
  /* voxel map */
  VImage map = VoxelMap(list,mask,nvox);


  /* read data into X */
  gsl_matrix_float *X = ReadDataECM(list,mask,map,first,length,type,nvox);



  /* main process */
  float *eigvec = VECM(X,(int)type,project,(VLong)seed,(int)maxiter,map);

  
  
  /* create output image */
  VImage dest = WriteOutput(list,map,eigvec,nvox);
  VSetAttr(VImageAttrList(dest),"name",NULL,VStringRepn,"ECM");
  VAttrList out_list = VCreateAttrList();
  VAppendAttr(out_list,"image",NULL,VImageRepn,dest);
  

  /* update geoinfo, 4D to 3D */
  if (geolist != NULL) {
    double *D = VGetGeoDim(geolist,NULL);
    D[0] = 3;  /* 3D */
    D[4] = 1;  /* just one timestep */
    VSetGeoDim(geolist,D);
    VSetGeoInfo(geolist,out_list);
  }
  
  /* write to disk */
  VHistory(VNumber(options),options,prg_name,&list,&out_list);
  if (! VWriteFile (out_file, out_list)) exit (1);
  fprintf (stderr, "%s: done.\n", argv[0]);
  return 0;
}
