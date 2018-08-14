/*
** Implementation of LISA algorithm
** for statistical inference of fMRI images
**
** single subject design
** no precoloring, no prewhitening
**
** G.Lohmann, April 2017
*/
#include "viaio/Vlib.h"
#include "viaio/file.h"
#include "viaio/mu.h"
#include "viaio/option.h"
#include "viaio/os.h"
#include <viaio/VImage.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_permutation.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/


#define MINVAL 1.0e+8

typedef struct TrialStruct {
  int   id;
  float onset;
  float duration;
  float height;
} Trial;

extern void VIsolatedVoxels(VImage src,float threshold);
extern void VHistogram(gsl_histogram *histogram,VString filename);
extern void VCheckImage(VImage src);
extern void FDR(VImage src,VImage dest,gsl_histogram *nullhist,gsl_histogram *realhist,double);
extern double ttest1(double *data1,int n);
extern void ImageStats(VImage src,double *,double *,double *hmin,double *hmax);
extern Trial *ReadDesign(VStringConst designfile,int *numtrials,int *nevents);
extern gsl_matrix *VCreateDesign(int ntimesteps,int nevents,int deriv,gsl_matrix *);
extern void VHemoModel(Trial *trial,int ntrials,int nevents,int ntimesteps,double tr,int deriv,gsl_matrix *X,gsl_matrix *);
extern Trial *CopyTrials(Trial *trial,int numtrials);
extern void VGLM(gsl_matrix *Data,gsl_matrix *X,gsl_matrix *XInv,gsl_vector *con,VImage map,VImage zmap);
extern void PlotDesign(gsl_matrix *X,double tr,VString filename);
extern Trial *ConcatenateTrials(Trial **trial,int *numtrials,float *run_duration,int dlists,int sumtrials);

extern double VImageVar(VImage src);
extern void VImageCount(VImage src);
extern void VBilateralFilter(VImage src,VImage,int radius,double var1,double var2,int);
extern void VGetHistRange(VImage src,double *hmin,double *hmax);
extern void VZScale(VImage src,float,float stddev);
extern float VGetMode(VImage src);
extern void GlobalMean(gsl_matrix *Data,gsl_matrix *covariates,int column);
extern gsl_matrix *VReadCovariates(VString cfile,VBoolean normalize);

extern VImage VoxelMap(VAttrList list);
extern gsl_matrix *VReadImageData(VAttrList *list,int nlists);
extern void VGetTimeInfos(VAttrList *list,int nlists,double *mtr,float *run_duration);
extern void VApplyMinvalNlists(VAttrList *list,int nlists,float minval);
extern void VRowNormalize(gsl_matrix *Data);
extern void CheckTrialLabels(Trial *trial,int numtrials);
extern void HistoUpdate(VImage,gsl_histogram *);


void XCheckImage(VImage src,char *filename)
{
  VAttrList out_list = VCreateAttrList();
  VAppendAttr(out_list,"image",NULL,VImageRepn,src);
  FILE *out_file = fopen(filename,"w");
  VWriteFile (out_file, out_list);
}




/* shuffle each run separately to ensure exchangebility, concatenate individual permtables */
int **genperm(gsl_rng *rx,int *numtrials,int sumtrials,int dlists,int numperm)
{
  int i,j,k;

  int **permtable = (int **) VCalloc(numperm,sizeof(int *));

  gsl_permutation **perm = (gsl_permutation **) VCalloc(dlists,sizeof(gsl_permutation *));
  for (k=0; k<dlists; k++) {
    perm[k] = gsl_permutation_alloc((size_t)numtrials[k]);
    gsl_permutation_init (perm[k]);
  }

  for (i = 0; i < numperm; i++) {
    permtable[i] = (int *) VCalloc(sumtrials,sizeof(int));
    int jj=0;
    for (k=0; k<dlists; k++) {
      gsl_ran_shuffle (rx, perm[k]->data,numtrials[k],sizeof(size_t));
      for (j=0; j<numtrials[k]; j++) {
	permtable[i][j+jj] = perm[k]->data[j] + jj;
      }
      jj += numtrials[k];
    }
  }

  for (k=0; k<dlists; k++) {
    gsl_permutation_free(perm[k]);
  }
  return permtable;
}




VDictEntry HemoDict[] = {
  { "gamma_0", 0 },
  { "gamma_1", 1 },
  { "gamma_2", 2 },
  { "gauss", 3 },
  { NULL }
};

int main (int argc, char *argv[])
{
  static VArgVector in_files;
  static VArgVector des_files;
  static VString  cova_filename="";
  static VString  out_filename="";
  static VFloat   minval = MINVAL;
  static VShort   hemomodel = 0;
  static VArgVector contrast;
  static VFloat   alpha = 0.05;
  static VShort   radius = 2;
  static VFloat   rvar = 2.0;
  static VFloat   svar = 2.0;
  static VShort   numiter = 2;
  static VBoolean cleanup = TRUE;
  static VBoolean verbose = FALSE;  
  static VBoolean globalmean = FALSE;
  static VShort   numperm = 2000;
  static VLong    seed = 99402622;
  static VShort   nproc = 0;
  static VOptionDescRec options[] = {
    {"in", VStringRepn, 0, & in_files, VRequiredOpt, NULL,"Input files" },
    {"design", VStringRepn, 0, & des_files, VRequiredOpt, NULL,"Design files (1st level)" },
    {"covariates", VStringRepn,  1, & cova_filename, VOptionalOpt, NULL,"Additional covariates (optional)" },
    {"out", VStringRepn, 1, & out_filename, VRequiredOpt, NULL,"Output file" },
    {"contrast", VFloatRepn, 0, (VPointer) &contrast, VRequiredOpt, NULL, "Contrast vector"},
    {"hemo", VShortRepn, 1, (VPointer) &hemomodel, VOptionalOpt, HemoDict,"Hemodynamic model" },
    {"alpha",VFloatRepn,1,(VPointer) &alpha,VOptionalOpt,NULL,"FDR significance level"},
    {"perm",VShortRepn,1,(VPointer) &numperm,VOptionalOpt,NULL,"Number of permutations"},
    {"seed",VLongRepn,1,(VPointer) &seed,VOptionalOpt,NULL,"Seed for random number generation"},    
    {"minval",VFloatRepn,1,(VPointer) &minval,VOptionalOpt,NULL,"Signal threshold"},  
    {"radius",VShortRepn,1,(VPointer) &radius,VOptionalOpt,NULL,"Bilateral parameter (radius in voxels)"},
    {"rvar",VFloatRepn,1,(VPointer) &rvar,VOptionalOpt,NULL,"Bilateral parameter (radiometric)"},
    {"svar",VFloatRepn,1,(VPointer) &svar,VOptionalOpt,NULL,"Bilateral parameter (spatial)"},
    {"filteriterations",VShortRepn,1,(VPointer) &numiter,VOptionalOpt,NULL,"Bilateral parameter (number of iterations)"},
    {"cleanup",VBooleanRepn,1,(VPointer) &cleanup,VOptionalOpt,NULL,"Whether to remove isloated voxels"},    
    {"gsr",VBooleanRepn,1,(VPointer) &globalmean,VOptionalOpt,NULL,"Global signal regression"},
    {"j",VShortRepn,1,(VPointer) &nproc,VOptionalOpt,NULL,"number of processors to use, '0' to use all"},
  };

  FILE *fp=NULL;
  VString in_filename;
  VAttrList out_list=NULL,geolist=NULL;
  int i;
  char *prg_name=GetLipsiaName("vslisa");
  fprintf (stderr, "%s\n", prg_name);


  /* parse command line */
  if (! VParseCommand (VNumber (options), options, & argc, argv)) {
    VReportUsage (argv[0], VNumber (options), options, NULL);
    exit (EXIT_FAILURE);
  }
  if (argc > 1) {
    VReportBadArgs (argc, argv);
    exit (EXIT_FAILURE);
  }


  /* omp-stuff */
#ifdef _OPENMP
  int num_procs=omp_get_num_procs();
  if (nproc > 0 && nproc < num_procs) num_procs = nproc;
  fprintf(stderr," using %d cores\n",(int)num_procs);
  omp_set_num_threads(num_procs);
#endif /* _OPENMP */




  /* read functional image data */
  int nlists = in_files.number;
  if (nlists < 1) VError(" no input");
  int dlists = des_files.number;
  if (dlists != nlists) VError(" number of input images and design files do not match: %d %d",nlists,dlists);

  VAttrList *list = (VAttrList *) VCalloc(nlists,sizeof(VAttrList));
  for (i=0; i<nlists; i++) {
    in_filename = ((VString *) in_files.vector)[i];
    fprintf(stderr," %3d:  %s\n",i,in_filename);
    list[i] = VReadAttrList(in_filename,0L,TRUE,FALSE);
    if (geolist == NULL) geolist = VGetGeoInfo(list[i]);
  }


  /* if default minval parameter is set, then compute new minval threshold for brain mask */
  if (fabs(minval - MINVAL) < 0.001) {
    minval = VGetMinvalNlists(list,nlists);
    fprintf(stderr," Brain mask,  minval: %f\n",minval);
  }
  VApplyMinvalNlists(list,nlists,minval);


  /* read data and voxel map */
  double tr=0;
  float *run_duration = (float *) VCalloc(nlists,sizeof(float));
  VGetTimeInfos(list,nlists,&tr,run_duration);
  gsl_matrix *Data = VReadImageData(list,nlists);
  VImage map  = VoxelMap(list[0]);
  int nslices = VPixel(map,0,3,0,VShort);
  int nrows   = VPixel(map,0,3,1,VShort);
  int ncols   = VPixel(map,0,3,2,VShort);
  int ntimesteps = Data->size2;


  /* additional regressors, no task labels, not included in permutations */
  gsl_matrix *ctmp1=NULL;
  gsl_matrix *ctmp2=NULL;
  gsl_matrix *covariates=NULL;
  int cdim = 1;
  if (strlen(cova_filename) > 1) {
    ctmp1 = VReadCovariates(cova_filename,TRUE);
    if (ctmp1->size1 != Data->size2) VError(" num timesteps in covariate file not consistent with data");
  }
  if (globalmean) {
    if (ctmp1 != NULL) cdim = ctmp1->size2+1;
    ctmp2 = gsl_matrix_calloc(Data->size2,cdim);
    GlobalMean(Data,ctmp2,(int)(cdim-1));
  }
  if (ctmp1 != NULL && ctmp2 == NULL) covariates = ctmp1;
  if (ctmp2 != NULL) covariates = ctmp2;


  /* design files with task labels */
  Trial **trial = (Trial **) VCalloc(dlists,sizeof(Trial *));
  int *numtrials = (int *) VCalloc(dlists,sizeof(int *));
  int nevents = 0;
  int sumtrials = 0;
  for (i=0; i<dlists; i++) {
    in_filename = ((VString *) des_files.vector)[i];
    fprintf(stderr," %3d:  %s\n",i,in_filename);
    int kk=0,jj=0;
    trial[i] = ReadDesign(in_filename,&kk,&jj);
    numtrials[i] = kk;
    if (nevents == 0) nevents = jj;
    if (nevents != jj) VError(" inconsistent number of event types in design %d",i);
    sumtrials += numtrials[i];
  }
  fprintf(stderr," number of trials: %d, number of event types: %d\n",sumtrials,nevents-1);
  Trial *alltrials = ConcatenateTrials(trial,numtrials,run_duration,nlists,sumtrials);
  CheckTrialLabels(alltrials,sumtrials);


  /* read contrast vector */
  gsl_vector *cont = gsl_vector_alloc(contrast.number+1);
  gsl_vector_set_zero(cont);
  for (i=0; i < contrast.number; i++) {
    double u = ((VFloat *)contrast.vector)[i];
    gsl_vector_set(cont,i+1,u);
  }



  /* alloc initial design matrix X */
  gsl_matrix *X = VCreateDesign(ntimesteps,nevents,(int)hemomodel,covariates);
  gsl_matrix *XInv = gsl_matrix_calloc(X->size2,X->size1);
  if (X->size2 != cont->size) 
    VError(" dimension of contrast vector does not match design matrix %ld %ld",X->size2,cont->size);

  

  /* ini random permutations */
  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *rx = gsl_rng_alloc(T);
  gsl_rng_set(rx,(unsigned long int)seed);
  if (verbose) fprintf(stderr," seed: %ld\n",(long)seed);
  int **permtable = genperm(rx,numtrials,sumtrials,dlists,(int)numperm);


  /* estimate null variance to adjust radiometric parameter, use first 30 permutations */
  int nperm=0;
  float stddev = 1.0;
  double meanvar = 0.0;
  if (numperm > 0) {
    int tstperm = 30;
    if (tstperm > numperm) tstperm = numperm;
    VImage zmap = VCreateImage(nslices,nrows,ncols,VFloatRepn);
    double varsum=0,nx=0;
    for (nperm = 0; nperm < tstperm; nperm++) {
      Trial *permtrials = CopyTrials(alltrials,sumtrials);
      int j=0;
      for (j=0; j<sumtrials; j++) {
	int j0 = permtable[nperm][j];
	permtrials[j].id = alltrials[j0].id;
      }
      gsl_matrix *X = VCreateDesign(ntimesteps,nevents,(int)hemomodel,covariates);
      gsl_matrix *XInv = gsl_matrix_calloc(X->size2,X->size1);
      VHemoModel(permtrials,sumtrials,nevents,ntimesteps,tr,(int)hemomodel,X,covariates);
      VGLM(Data,X,XInv,cont,map,zmap);
      varsum += VImageVar(zmap);
      nx++;

      gsl_matrix_free(X);
      gsl_matrix_free(XInv);
      VFree(permtrials);
    }
    meanvar = varsum/nx;
    stddev = (float)(sqrt(meanvar));  /* update stddev */
    VDestroyImage(zmap);
  }


  /* no permutation */
  VImage zmap1 = VCreateImage(nslices,nrows,ncols,VFloatRepn);
  VCopyImageAttrs (map,zmap1);
  VImage dst1 = VCreateImageLike (zmap1);
  VHemoModel(alltrials,sumtrials,nevents,ntimesteps,tr,(int)hemomodel,X,covariates);
  VGLM(Data,X,XInv,cont,map,zmap1);


  if (numperm == 0) {
    double z = VImageVar(zmap1);
    stddev = sqrt(z); /* update stddev */
  }
  float mode=0;
  if (numperm > 0) VZScale(zmap1,mode,stddev);
  VBilateralFilter(zmap1,dst1,(int)radius,(double)rvar,(double)svar,(int)numiter);


  /* ini histograms */
  double hmin=0,hmax=0;
  VGetHistRange(dst1,&hmin,&hmax);
  if (verbose) fprintf(stderr," Histogram range:  [%.3f, %.3f]\n",hmin,hmax);
  size_t nbins = 10000;
  gsl_histogram *hist0 = gsl_histogram_alloc (nbins);
  gsl_histogram_set_ranges_uniform (hist0,hmin,hmax);
  gsl_histogram *histz = gsl_histogram_alloc (nbins);
  gsl_histogram_set_ranges_uniform (histz,hmin,hmax);
  HistoUpdate(dst1,histz);


  /* random permutations */
#pragma omp parallel for shared(Data) schedule(dynamic)
  for (nperm = 0; nperm < numperm; nperm++) {
    if (nperm%5 == 0) fprintf(stderr," perm  %4d  of  %d\r",nperm,(int)numperm);

    /* randomly shuffle trial labels */
    Trial *permtrials = CopyTrials(alltrials,sumtrials);
    int j=0;
    for (j=0; j<sumtrials; j++) {
      int j0 = permtable[nperm][j];
      permtrials[j].id = alltrials[j0].id;
    }  

    /* hemodynamic model */
    gsl_matrix *X = VCreateDesign(ntimesteps,nevents,(int)hemomodel,covariates);
    gsl_matrix *XInv = gsl_matrix_calloc(X->size2,X->size1);
    VHemoModel(permtrials,sumtrials,nevents,ntimesteps,tr,(int)hemomodel,X,covariates);
    

    /* GLM */
    VImage zmap = VCreateImageLike(zmap1);
    VGLM(Data,X,XInv,cont,map,zmap);
    
    
    VZScale(zmap,mode,stddev);
    gsl_matrix_free(X);
    gsl_matrix_free(XInv);
    VFree(permtrials);


    /* bilateral filter */
    VImage dst = VCreateImageLike (zmap);
    VBilateralFilter(zmap,dst,(int)radius,(double)rvar,(double)svar,(int)numiter);

#pragma omp critical 
    {
      HistoUpdate(dst,hist0);
    }
    VDestroyImage(dst);
    VDestroyImage(zmap);
  }
  

  /* apply fdr */
  VImage fdrimage = VCopyImage (dst1,NULL,VAllBands);
  if (numperm > 0) {
    FDR(dst1,fdrimage,hist0,histz,(double)alpha);

    if (cleanup && alpha < 1.0) {
      VIsolatedVoxels(fdrimage,(float)(1.0-alpha));
    }
  }



  /*
  ** output
  */
  out_list = VCreateAttrList ();
  VHistory(VNumber(options),options,prg_name,&list[0],&out_list);

  /* update geoinfo, 4D to 3D */
  if (geolist != NULL) {
    double *D = VGetGeoDim(geolist,NULL);
    D[0] = 3;
    D[4] = 1;
    VSetGeoDim(geolist,D);
  }
  VSetGeoInfo(geolist,out_list);
  VAppendAttr (out_list,"image",NULL,VImageRepn,fdrimage);
  fp = VOpenOutputFile (out_filename, TRUE);
  if (! VWriteFile (fp, out_list)) exit (1);
  fclose(fp);
  fprintf (stderr, "\n%s: done.\n", argv[0]);
  exit(0);
}
