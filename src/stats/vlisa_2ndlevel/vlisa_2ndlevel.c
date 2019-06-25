/*
** Implementation of LISA algorithm
** for statistical inference of fMRI images
**
** 2nd level inference using a GLM design matrix, 
**   e.g. for anova designs
**
** G.Lohmann, 2018
*/
#include <viaio/Vlib.h>
#include <viaio/file.h>
#include <viaio/mu.h>
#include <viaio/option.h>
#include <viaio/os.h>
#include <viaio/VImage.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/

#define ABS(x) ((x) > 0 ? (x) : -(x))

extern void   VIsolatedVoxels(VImage src,float threshold);
extern void   VHistogram(gsl_histogram *histogram,VString filename);
extern void   VCheckImage(VImage src);
extern void   FDR(VImage src,VImage dest,gsl_histogram *nullhist,gsl_histogram *realhist,double);
extern double ttest1(double *data1,int n);
extern void   ImageStats(VImage src,double *,double *,double *hmin,double *hmax);
extern void   VBilateralFilter(VImage src,VImage dest,int radius,double var1,double var2,int);
extern double VImageVar(VImage src);
extern void   VImageCount(VImage src);
extern void   VGetHistRange(VImage src,double *hmin,double *hmax);
extern void   VZScale(VImage src,float mode,float stddev);
extern float  VGetMode(VImage src);

extern void XReadExchange(VString filename,int *exchange,int n);
extern void GLM2(VImage *src,gsl_matrix *X,gsl_vector *contrast,int *,int *permtable,int *signtable,int signswitch,VImage dest);
extern gsl_matrix *XRead2ndLevel(VString);
extern gsl_matrix *PseudoInv(gsl_matrix *A,gsl_matrix *B);
extern int SignSwitch (gsl_matrix *X,gsl_vector *contrast,int *);
extern void genperm(long seed,int *exchange,int,int **,int **,int,int,gsl_vector *);
extern void HistoUpdate(VImage,gsl_histogram *);
extern gsl_matrix *VReadCovariates(VString,VBoolean);



int main (int argc, char *argv[])
{
  static VArgVector in_files1;
  static VString  out_filename="";  
  static VString  design_filename="";
  static VString  exchange_filename="";
  static VString  nuisance_filename="";
  static VString  mask_filename="";
  static VArgVector cont;
  static VFloat   alpha = 0.05; 
  static VShort   radius = 2;
  static VFloat   rvar = 2.0;
  static VFloat   svar = 2.0;
  static VShort   numiter = 2;
  static VShort   numperm = 5000;
  static VLong    seed = 99402622;
  static VBoolean centering = FALSE;
  static VBoolean demean  = TRUE;
  static VBoolean cleanup = TRUE;
  static VShort   nproc = 0;
  static VOptionDescRec options[] = {
    {"in", VStringRepn, 0, & in_files1, VRequiredOpt, NULL,"Input files" },
    {"out", VStringRepn, 1, & out_filename, VRequiredOpt, NULL,"Output file" },   
    {"design",VStringRepn,1,(VPointer) &design_filename,VRequiredOpt,NULL,"Design file (2nd level)"}, 
    {"grp", VStringRepn, 1, & exchange_filename, VOptionalOpt, NULL,"Group Ids needed for exchangeability" },
    {"contrast", VFloatRepn, 0, (VPointer) &cont, VRequiredOpt, NULL, "contrast vector"},
    {"nuisance", VStringRepn, 1, & nuisance_filename, VOptionalOpt, NULL,"Nuisance regressors" },
    {"demean",VBooleanRepn,1,(VPointer) &demean,VOptionalOpt,NULL,"Whether to subtract mean in nuisance regressors"},
    {"alpha",VFloatRepn,1,(VPointer) &alpha,VOptionalOpt,NULL,"FDR significance level"},
    {"perm",VShortRepn,1,(VPointer) &numperm,VOptionalOpt,NULL,"Number of permutations"},
    {"mask", VStringRepn, 1, (VPointer) &mask_filename, VRequiredOpt, NULL, "Mask"},   
    {"seed",VLongRepn,1,(VPointer) &seed,VOptionalOpt,NULL,"Seed for random number generation"},
    {"radius",VShortRepn,1,(VPointer) &radius,VOptionalOpt,NULL,"Bilateral parameter (radius in voxels)"},
    {"rvar",VFloatRepn,1,(VPointer) &rvar,VOptionalOpt,NULL,"Bilateral parameter (radiometric)"},
    {"svar",VFloatRepn,1,(VPointer) &svar,VOptionalOpt,NULL,"Bilateral parameter (spatial)"},
    {"filteriterations",VShortRepn,1,(VPointer) &numiter,VOptionalOpt,NULL,"Bilateral parameter (number of iterations)"},
    {"cleanup",VBooleanRepn,1,(VPointer) &cleanup,VOptionalOpt,NULL,"Whether to remove isolated voxels"},
    {"j",VShortRepn,1,(VPointer) &nproc,VOptionalOpt,NULL,"Number of processors to use, '0' to use all"},
  };
  FILE *fp=NULL;
  VString in_filename;
  VAttrList list1=NULL,out_list=NULL,geolist=NULL;
  VImage *src1;
  int i,j,nimages,npix=0;
  double u=0;
  char *prg_name=GetLipsiaName("vlisa_2ndlevel");
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


  /* images  */
  VImage mask = VReadImageFile(mask_filename);
  if (mask==NULL) VError("Error reading mask file %s",mask_filename);

  nimages = in_files1.number;
  src1 = (VImage *) VCalloc(nimages,sizeof(VImage));
  for (i = 0; i < nimages; i++) {
    in_filename = ((VString *) in_files1.vector)[i];
    list1   = VReadAttrList(in_filename,0L,TRUE,FALSE);
    VMaskMinval(list1,mask,0.0);
    src1[i] = VReadImage(list1);
    if (src1[i] == NULL) VError(" no input image found");
    if (VPixelRepn(src1[i]) != VFloatRepn) VError(" input pixel repn must be float");
    if (i == 0) npix = VImageNPixels(src1[i]);
    else if (npix != VImageNPixels(src1[i])) VError(" inconsistent image dimensions");
    fprintf(stderr," %3d:  %s\n",i,in_filename);

    /* use geometry info from 1st file */
    if (geolist == NULL) geolist = VGetGeoInfo(list1);
  }


  /* get nuisance file */
  gsl_matrix *XN = NULL;
  if (strlen(nuisance_filename) > 1) {
    XN = VReadCovariates(nuisance_filename,demean);
    fprintf(stderr," nuisance covariates dimensions:  %d x %d\n",(int)XN->size1,(int)XN->size2);
    if (XN->size1 != nimages) 
      VError(" number of input images (%d) does not match number of rows in nuisance covariates file (%d)",nimages,(int)XN->size1);
  }


  
  /* get design file and its inverse */
  gsl_matrix *X0 = XRead2ndLevel(design_filename);
  fprintf(stderr," design file dimensions:  %d x %d\n",(int)X0->size1,(int)X0->size2);
  if (X0->size1 != nimages) 
    VError(" number of input images (%d) does not match number of rows in design file (%d)",nimages,(int)X0->size1);


  /* read contrast vector */
  if (cont.number != X0->size2) 
    VError(" Dimension of contrast vector (%d) does not match number of columns in design file (%d)",cont.number,X0->size2);
  gsl_vector *contrast0 = gsl_vector_alloc(cont.number);
  for (i=0; i < cont.number; i++) {
    double u = ((VFloat *)cont.vector)[i];
    gsl_vector_set(contrast0, i, u);
  }


  /* concatenate design matrix and nuisance covariates */
  gsl_matrix *X = NULL;
  gsl_vector *contrast = NULL;
  int *permflag = NULL;

  /* include nuisance covariates, if present */
  if (XN != NULL) {
    permflag = (int *) VCalloc(XN->size2 + X0->size2,sizeof(int));
    for (j=0; j<XN->size2 + X0->size2; j++) permflag[j] = 0;
    for (j=0; j<X0->size2; j++) permflag[j] = 1;
    
    X = gsl_matrix_calloc(X0->size1, X0->size2 + XN->size2);
    for (i=0; i<X->size1; i++) {
      for (j=0; j<X0->size2; j++) {
	u = gsl_matrix_get(X0,i,j);
	gsl_matrix_set(X,i,j,u);	
      }
      for (j=0; j<XN->size2; j++) {
	u = gsl_matrix_get(XN,i,j);
	gsl_matrix_set(X,i,j+X0->size2,u);
      }
    }
    contrast = gsl_vector_calloc(cont.number + XN->size2);
    gsl_vector_set_zero(contrast);
    for (j=0; j<contrast0->size; j++) gsl_vector_set(contrast,j,contrast0->data[j]);
  }

  /* no nuisance covariates */
  else {
    if (X0->size2 == 1) VError(" Please use 'vlisa_onesample' for onesample tests");
    permflag = (int *) VCalloc(X0->size2,sizeof(int));
    for (j=0; j<X0->size2; j++) permflag[j] = 1;
    X = gsl_matrix_calloc(X0->size1,X0->size2);
    gsl_matrix_memcpy(X,X0);
    contrast = gsl_vector_calloc(cont.number);
    gsl_vector_memcpy(contrast,contrast0);
  }

 
  
  /* read exchangeability information */
  int *exchange = (int *) VCalloc(X->size1,sizeof(int));
  for (i=0; i<X->size1; i++) exchange[i] = 1;  /* no restrictions w.r.t. exchangeability */
  if (strlen(exchange_filename) > 1) { 
    XReadExchange(exchange_filename,exchange,(int)X->size1);
  }


  /* ini random permutations and sign switching */
  int signswitch = SignSwitch(X,contrast,permflag);
  int nperm=0;
  int **permtable = (int **) VCalloc(numperm,sizeof(int *));
  int **signtable = (int **) VCalloc(numperm,sizeof(int *));
  genperm((long)seed,exchange,signswitch,permtable,signtable,(int)nimages,(int)numperm,contrast);



  /* estimate null variance based on first 30 permutations */
  double hmin=0,hmax=0;
  float stddev=1.0;
  if (numperm > 0) {
    int tstperm = 30;
    if (tstperm > numperm) tstperm = numperm;
    double varsum=0,nx=0;

#pragma omp parallel for shared(src1,permtable) schedule(dynamic)
    for (nperm = 0; nperm < tstperm; nperm++) {
      VImage zmap = VCreateImageLike(src1[0]);
      GLM2(src1,X,contrast,permflag,permtable[nperm],signtable[nperm],signswitch,zmap);
#pragma omp critical 
      {
	varsum += VImageVar(zmap);
	nx++;
      }
      VDestroyImage(zmap);
    }
    double meanvar = varsum/nx;
    stddev = sqrt(meanvar);
  }

  

  /* no permutation */
  int *nopermtable = (int *) VCalloc(nimages,sizeof(int));
  int *nosigntable = (int *) VCalloc(nimages,sizeof(int));
  for (i=0; i<nimages; i++) nopermtable[i] = i;
  for (i=0; i<nimages; i++) nosigntable[i] = 1;
  VImage dst1  = VCreateImageLike (src1[0]);
  VImage zmap1 = VCreateImageLike(src1[0]);
  GLM2(src1,X,contrast,permflag,nopermtable,nosigntable,signswitch,zmap1);


  if (numperm == 0) {
    double z = VImageVar(zmap1);
    stddev = (float)(sqrt(z)); /* update stddev */
  }
  float mode=0;
  if (centering) mode = VGetMode(zmap1); 
  if (numperm > 0) VZScale(zmap1,mode,stddev);
  VBilateralFilter(zmap1,dst1,(int)radius,(double)(rvar),(double)svar,(int)numiter);


  /* ini histograms */
  VGetHistRange(dst1,&hmin,&hmax);
  size_t nbins = 20000;
  gsl_histogram *hist0 = gsl_histogram_alloc (nbins);
  gsl_histogram_set_ranges_uniform (hist0,hmin,hmax);
  gsl_histogram *histz = gsl_histogram_alloc (nbins);
  gsl_histogram_set_ranges_uniform (histz,hmin,hmax);
  HistoUpdate(dst1,histz);



  /* random permutations */
#pragma omp parallel for shared(src1,permtable) schedule(dynamic)
  for (nperm = 0; nperm < numperm; nperm++) {
    if (nperm%20 == 0) fprintf(stderr," perm  %4d  of  %d\r",nperm,(int)numperm);

    VImage zmap = VCreateImageLike(src1[0]);
    VImage dst  = VCreateImageLike (zmap);      
    GLM2(src1,X,contrast,permflag,permtable[nperm],signtable[nperm],signswitch,zmap);
    float mode=0;
    if (centering) mode = VGetMode(zmap);
    VZScale(zmap,mode,stddev);
    VBilateralFilter(zmap,dst,(int)radius,(double)(rvar),(double)svar,(int)numiter);

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



  /* write output to disk */
  out_list = VCreateAttrList ();
  VHistory(VNumber(options),options,prg_name,&list1,&out_list);
  VSetGeoInfo(geolist,out_list);
  VAppendAttr (out_list,"image",NULL,VImageRepn,fdrimage);
  fp = VOpenOutputFile (out_filename, TRUE);
  if (! VWriteFile (fp, out_list)) exit (1);
  fclose(fp);
  fprintf (stderr, "\n%s: done.\n", argv[0]);
  exit(0);
}
