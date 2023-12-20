/*
** Semi-blind machine learning with ensemble learning (SML-EL)
**
** G.Lohmann, MPI-KYB, Aug 2023
*/
#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_fit.h>


#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/

extern void VCheckImage(VImage src);
extern void RankNorm(gsl_vector *x);
extern void RNorm(double *x,size_t n);
extern void ZNorm(gsl_vector *);
extern void YNorm(double *x,size_t n);
extern void MatNorm0(gsl_matrix *X);
extern void SelectTrainTest(int *train,int *test,size_t ntrain,size_t ntest,size_t numsubjects,gsl_rng *rx);

extern gsl_vector *ReadRegressor(VString filename);
extern VImage ReadSelectionMap(VAttrList sel_list);
extern long *ReadSelectionI(VAttrList sel_list);
extern long *ReadSelectionJ(VAttrList sel_list);
extern VImage ReadSelectionImage(VAttrList sel_list);
extern int VReadCoremaps(VArgVector in_files,size_t numsubjects,int ftype,VImage *coremap);
extern size_t GetNumEdges(VArgVector in_files);
extern size_t GetNumComponents(VArgVector in_files);
extern void VReadBundles(VArgVector in_files,size_t nedges,size_t numsubjects,
			 gsl_matrix_float *X0,gsl_matrix_float *X1,size_t *n0,size_t *n1);
extern void RowNormalizeBundle(gsl_matrix_float *X);
extern void MatColumnCentering(gsl_matrix *X);

extern gsl_matrix *ReadConfounds(VArgVector confound_files,int numsubjects,FILE *,int);
extern void Confounds(gsl_matrix *C,gsl_vector *y,gsl_vector *bx,double,int *table,int flag);

extern VImage ReadMask(VString mask_filename,int verbose);
extern size_t ApplyMask(VImage mask,VImage map,long *I,long *J,char *table,size_t,int);

extern void Corrmap(gsl_matrix_float *X,gsl_vector *ytrain,int *train,double *W);

extern void XResidualFit(gsl_vector *,gsl_vector *,gsl_vector *,double *corr,double *r2,double *,double *);
extern gsl_vector *Multifit(gsl_matrix *X,gsl_vector *regressor);
extern double XFit(gsl_matrix *X,gsl_vector *ytrain,gsl_vector *ytest,
		   gsl_vector *xtrain,gsl_vector *xtest,int *train,int *test);
extern void Fit(gsl_vector *xdata,gsl_vector *ydata,gsl_vector *zdata);


extern size_t Sel0(char *etable,size_t nedges,gsl_rng *rx,gsl_matrix_long *S);
extern double SelProb(gsl_matrix_float *Xcorr,gsl_vector *,gsl_vector *,gsl_vector *,gsl_vector *,
		      int *train,int *test,char *etable,gsl_matrix_long *S,gsl_rng *rx,double,int);

extern double Loop(gsl_matrix_float *Xcorr,
		   gsl_vector *ytrain,gsl_vector *ytest,
		   int *train,int *test,gsl_matrix_long *S,size_t npls,gsl_vector *);

extern double Check(gsl_vector *ytest,gsl_vector *result,int type);
extern void Corrmap(gsl_matrix_float *X,gsl_vector *ytrain,int *train,double *W);

extern int Split(gsl_vector *ytrain,gsl_vector *xtrain,gsl_vector *ytest,gsl_vector *xtest,
		 gsl_vector *sytrain,gsl_vector *sxtrain,gsl_vector *sytest,gsl_vector *sxtest,
		 int *train,int *test,
		 int *strain,int *stest,size_t ktrain,size_t ktest,gsl_rng *rx);


double PearsonCorr(gsl_vector *ytest,gsl_vector *fitted)
{
  if (ytest==NULL || fitted==NULL) return 0;
  return gsl_stats_correlation(ytest->data,1,fitted->data,1,ytest->size);
}



/* univariate fit */
void PredictFit(gsl_vector *xtrain,gsl_vector *ytrain,gsl_vector *xtest)
{
  size_t i;
  double c0,c1,cov00,cov01,cov11,sumsq,y,y_err;

  gsl_fit_linear(xtrain->data,1,ytrain->data,1,xtrain->size,&c0,&c1,&cov00,&cov01,&cov11,&sumsq);
  for (i=0; i<xtrain->size; i++) {
    gsl_fit_linear_est(xtrain->data[i],c0,c1,cov00,cov01,cov11,&y,&y_err);
    xtrain->data[i] = y;
  }
  for (i=0; i<xtest->size; i++) {
    gsl_fit_linear_est(xtest->data[i],c0,c1,cov00,cov01,cov11,&y,&y_err);
    xtest->data[i] = y;
  }
}




/* coefficient of determination */
double CoD(gsl_vector *observed,gsl_vector *fitted)
{
  size_t i,n=observed->size;   
  double *y = observed->data;
  double *f = fitted->data;
  
  double mean = gsl_stats_mean(y,1,n);
  double ss_total=0;
  for (i=0; i<n; i++) ss_total += (y[i]-mean)*(y[i]-mean);

  double ss_res=0;
  for (i=0; i<n; i++) ss_res += (y[i]-f[i])*(y[i]-f[i]);

  double r2 = 0;
  if (ss_total > 0) r2 = 1.0 - ss_res/ss_total;
  return r2;
}


int main (int argc, char *argv[])
{
  static VArgVector train_files;
  static VArgVector test_files;
  static VString out_filename="";
  static VString ytrain_filename="";
  static VString ytest_filename="";
  static VString xtrain_filename="";
  static VString xtest_filename="";
  static VLong  xdim       = 500;
  static VLong  nensembles = 1000;
  static VShort xpls    = 10;
  static VLong  seed    = 601931;
  static VShort nproc = 0;
  static VOptionDescRec options[] = {
    { "train", VStringRepn, 0, & train_files, VRequiredOpt, NULL,"Training files" },
    { "test", VStringRepn, 0, & test_files, VRequiredOpt, NULL,"Test files" },
    { "ytrain", VStringRepn, 1, & ytrain_filename, VRequiredOpt, NULL,"Training Regressor" },
    { "ytest", VStringRepn, 1, & ytest_filename, VOptionalOpt, NULL,"Test Regressor" },
    { "xtrain", VStringRepn, 1, & xtrain_filename, VRequiredOpt, NULL,"Training supplementary info" },
    { "xtest", VStringRepn, 1, & xtest_filename, VRequiredOpt, NULL,"Test supplementary info" },
    { "out", VStringRepn, 1, & out_filename, VRequiredOpt, NULL,"Output txt-file" },
    { "dimX", VLongRepn, 1, & xdim, VOptionalOpt, NULL,"Dimension of matrix" },
    { "npls", VShortRepn, 1, & xpls, VOptionalOpt, NULL,"Number of components in PLS" },
    { "nensembles", VLongRepn, 1, & nensembles, VOptionalOpt, NULL,"Number of ensembles" },
    { "seed", VLongRepn, 1, & seed, VOptionalOpt, NULL,"Seed for random number generator" },
    { "j",VShortRepn,1,(VPointer) &nproc,VOptionalOpt,NULL,"Number of processors to use, '0' to use all"}
  };
  size_t i,j;

  /* parse command line */
  if (! VParseCommand (VNumber (options), options, & argc, argv)) {
    VReportUsage (argv[0], VNumber (options), options, NULL);
    exit (EXIT_FAILURE);
  }
  if (argc > 1) {
    VReportBadArgs (argc, argv);
    exit (EXIT_FAILURE);
  }
  gsl_set_error_handler_off();

#ifdef _OPENMP
  int num_procs=omp_get_num_procs();
  if (nproc > 0 && nproc < num_procs) num_procs = nproc;
  omp_set_num_threads(num_procs);
#endif /* _OPENMP */

 
  /* hyperparameters */
  size_t ntrain = train_files.number;
  size_t ntest = test_files.number;
  size_t dimX  = (size_t)xdim;
  size_t npls  = (size_t)xpls; 
  size_t niter = (size_t)nensembles;
  size_t numsubjects = ntrain+ntest;
  fprintf(stderr,"# ntrain,ntest: %lu %lu\n",ntrain,ntest);


  /* random generator */
  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *rx = gsl_rng_alloc(T);
  gsl_rng_set(rx,(unsigned long int)seed);


  /* read regressor files */
  gsl_vector *ytrain = ReadRegressor(ytrain_filename);
  gsl_vector *ytest = NULL;
  int testflag=0;
  if (strlen(ytest_filename) > 1) { ytest = ReadRegressor(ytest_filename); testflag = 1; }
  else ytest = gsl_vector_calloc(ntest);
  if (ytrain->size != ntrain) VError(" inconsistent training set sizes");
  if (ytest->size != ntest) VError(" inconsistent test set sizes");

  
  /* read supplementary non-imaging info */
  gsl_vector *xtrain = ReadRegressor(xtrain_filename);
  gsl_vector *xtest = ReadRegressor(xtest_filename);
  if (xtrain->size != ntrain) VError(" inconsistent training set sizes");
  if (xtest->size != ntest) VError(" inconsistent test set sizes");
  gsl_vector *xregressor = gsl_vector_calloc(numsubjects);
  for (i=0; i<ntrain; i++) xregressor->data[i] = xtrain->data[i];
  for (i=0; i<ntest; i++) xregressor->data[i+ntrain] = xtest->data[i];
  ZNorm(xregressor);
  for (i=0; i<ntrain; i++) xtrain->data[i] = xregressor->data[i];
  for (i=0; i<ntest; i++) xtest->data[i] = xregressor->data[i+ntrain];
  
  
 
  /* read data */
  fprintf(stderr,"# reading data...\n");
  size_t nedges = GetNumEdges(train_files);
  size_t xedges = GetNumEdges(test_files);
  if (nedges != xedges) VError(" number of edges in training/test data do not match (%lu %lu)",nedges,xedges);
  gsl_matrix_float *DataTrain = gsl_matrix_float_calloc(ntrain,nedges);
  if (!DataTrain) VError(" err allocating DataTrain");
  gsl_matrix_float *DataTest = gsl_matrix_float_calloc(ntest,nedges);
  if (!DataTest) VError(" err allocating DataTest");
  
  size_t k0=0,k1=0;
  VReadBundles(train_files,nedges,(size_t)ntrain,DataTrain,NULL,&k0,&k1);
  if (k0 == 0 && k1 == 0) VError(" err reading input bundles");
  RowNormalizeBundle(DataTrain);

  VReadBundles(test_files,nedges,(size_t)ntest,DataTest,NULL,&k0,&k1);
  if (k0 == 0 && k1 == 0) VError(" err reading input bundles");
  RowNormalizeBundle(DataTest);

  
  /* concatenate */
  float *p0,*p1;
  gsl_matrix_float *Xcorr = gsl_matrix_float_calloc(numsubjects,nedges);
  for (i=0; i<ntrain; i++) {
    p0 = gsl_matrix_float_ptr(DataTrain,i,0);
    p1 = gsl_matrix_float_ptr(Xcorr,i,0);
    for (j=0; j<nedges; j++) p1[j] = p0[j];
  }
  for (i=0; i<ntest; i++) {
    p0 = gsl_matrix_float_ptr(DataTest,i,0);
    p1 = gsl_matrix_float_ptr(Xcorr,ntrain+i,0);
    for (j=0; j<nedges; j++) p1[j] = p0[j];
  }
  
  int *ptrain = (int *)VCalloc(ntrain,sizeof(int));
  int *ptest  = (int *)VCalloc(ntest,sizeof(int));
  for (i=0; i<ntrain; i++) ptrain[i] = i;
  for (i=0; i<ntest; i++) ptest[i] = ntrain+i;
  
  char *etable = (char *) VCalloc(nedges,sizeof(char));
  for (i=0; i<nedges; i++) etable[i] = 1;
  
  gsl_vector *aresult  = gsl_vector_calloc(ntest);
  gsl_vector *cresult  = gsl_vector_calloc(ntest);
  gsl_vector *xresult  = gsl_vector_calloc(ntest);


  /* mult */
  double xmult = 0.5;
  size_t mult = 3;
  size_t n = ntrain+ntest;
  if (ntrain == ntest) { xmult = -1; mult = 1; }
  
  size_t ktrain = ntrain;
  size_t ktest = 0;
  if (xmult < 0) {
    ktrain = ntrain;
    ktest = ntest;
    mult=1;
  }
  else {
    ktrain = (size_t)(xmult * (float)n);
    ktest = n-ktrain;
  }
  int *strain = (int *)VCalloc(ktrain,sizeof(int));
  int *stest = (int *)VCalloc(ktest,sizeof(int));
  gsl_vector *sytrain  = gsl_vector_calloc(ktrain);
  gsl_vector *sxtrain  = gsl_vector_calloc(ktrain);
  gsl_vector *sytest  = gsl_vector_calloc(ktest);
  gsl_vector *sxtest  = gsl_vector_calloc(ktest);
  
  
  /* ensembles */
  gsl_matrix_long *S0 = gsl_matrix_long_calloc(niter,dimX);
  size_t niter2 = niter;
  if (mult > 1) niter2 = (size_t)((float)niter/(float)mult);
  gsl_matrix_long *S2 = gsl_matrix_long_calloc(niter2,dimX);
  if (ktrain+ktest != n) VError("n: %lu,  %lu",n,ktrain+ktest);
  double ya1=0,ya3=0,yx1=0,yx3=0,yc1=0,yc3=0;
  
  /* ini result structs */
  double ymean = gsl_stats_mean(ytrain->data,1,ytrain->size);
  double ysd = gsl_stats_sd(ytrain->data,1,ytrain->size);
  gsl_vector_add_constant(ytrain,-ymean);
  gsl_vector_scale(ytrain,1.0/ysd);
  if (testflag > 0) {
    gsl_vector_add_constant(ytest,-ymean);
    gsl_vector_scale(ytest,1.0/ysd);
  }

  
  /* A:  prediction without using supplementary info */
  Sel0(etable,nedges,rx,S0);
  Loop(Xcorr,ytrain,NULL,ptrain,ptest,S0,npls,aresult);
  if (testflag > 0) {
    ya1 = PearsonCorr(ytest,aresult);
    ya3 = CoD(ytest,aresult);
    fprintf(stderr,"# Prediction results (without supplementary info, BML),  corr: %.4f,  R2: %.4f\n",ya1,ya3);
  }

  /* X: prediction based solely on supplementary info */
  PredictFit(xtrain,ytrain,xtest);
  if (testflag > 0) {
    yx1 = PearsonCorr(ytest,xtest); 
    yx3 = CoD(ytest,xtest);
    fprintf(stderr,"# Prediction results (supplementary info only, SI),      corr: %.4f,  R2: %.4f\n",yx1,yx3);
  }

  /* C: prediction with bias-control, alpha is adjusted */
  double zcorr=0,xdiff=100,up=1.2,down=0.8,xalpha=1;
  int k=0,kiter=30;
  
  double xcorr = PearsonCorr(ytrain,xtrain);
  if (xcorr < 0) VError(" xcorr < 0");
  while (fabs(xdiff) > 0.02 && xalpha > 0.01 && k < kiter) {
    gsl_vector_set_zero(cresult);
    gsl_vector_set_zero(xresult);

    for (j=0; j<mult; j++) {
      int rtcode = Split(ytrain,xtrain,xtest,xtest,sytrain,sxtrain,sytest,sxtest,ptrain,ptest,strain,stest,ktrain,ktest,rx);
      if (rtcode < 0) VError("split: %d",rtcode);
      SelProb(Xcorr,sytrain,sxtrain,sytest,sxtest,strain,stest,etable,S2,rx,(double)xalpha,0);
      Loop(Xcorr,ytrain,NULL,ptrain,ptest,S2,npls,xresult);
      gsl_vector_add(cresult,xresult);
    }
    gsl_vector_scale(cresult,1.0/(double)mult);
    zcorr = PearsonCorr(cresult,xtest);
    xdiff = zcorr-xcorr;
    if (k > 1) { up=1.1; down=0.9; }
    if (zcorr > xcorr) xalpha *= down;
    if (zcorr < xcorr) xalpha *= up;
    k++;
  }
  if (testflag > 0) {
    yc1 = PearsonCorr(ytest,cresult);
    yc3 = CoD(ytest,cresult);
    fprintf(stderr,"# Prediction results (with supplementary info, SML),     corr: %.4f,  R2: %.4f\n",yc1,yc3);
  }
  
  FILE *fp = fopen(out_filename,"w");
  fprintf(stderr,"#\n         Observed      SI        BML        SML\n\n");
  fprintf(fp,"#\n         Observed      SI        BML        SML\n\n");
  for (i=0; i<ntest; i++) {
    fprintf(stderr," %5lu  %9.6f  %9.6f  %9.6f  %9.6f\n",i,ytest->data[i],xtest->data[i],aresult->data[i],cresult->data[i]);
    fprintf(fp," %5lu  %9.6f  %9.6f  %9.6f  %9.6f\n",i,ytest->data[i],xtest->data[i],aresult->data[i],cresult->data[i]);
  }
  fclose(fp);

  fprintf (stderr, "%s: done.\n", argv[0]);
  exit(0);
}
