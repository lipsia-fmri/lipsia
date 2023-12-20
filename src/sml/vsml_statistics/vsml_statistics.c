/*
** Semi-blind machine learning with ensemble learning (SML-EL)
** Statistics over many different training/test splits 
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


#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/

extern void ZNorm(gsl_vector *);
extern void SelectTrainTest(int *train,int *test,size_t ntrain,size_t ntest,size_t numsubjects,gsl_rng *rx);

extern size_t GetNumEdges(VArgVector in_files);
extern size_t GetNumComponents(VArgVector in_files);
extern void VReadBundles(VArgVector in_files,size_t nedges,size_t numsubjects,
			 gsl_matrix_float *X0,gsl_matrix_float *X1,size_t *n0,size_t *n1);
extern void RowNormalizeBundle(gsl_matrix_float *);
extern gsl_matrix *ReadConfounds(VArgVector confound_files,int numsubjects,FILE *,int);
extern void MatColumnCentering(gsl_matrix *);
extern gsl_vector *ReadRegressor(VString filename);

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

extern double PearsonCorr(gsl_vector *ytest,gsl_vector *result);
extern double SpearmanCorr(gsl_vector *ytest,gsl_vector *result);
extern double CoD(gsl_vector *observed,gsl_vector *fitted);
extern double Accuracy(gsl_vector *ytest,gsl_vector *fitted);

extern int Split(gsl_vector *ytrain,gsl_vector *xtrain,gsl_vector *ytest,gsl_vector *xtest,
		 gsl_vector *sytrain,gsl_vector *sxtrain,gsl_vector *sytest,gsl_vector *sxtest,
		 int *train,int *test,
		 int *strain,int *stest,size_t ktrain,size_t ktest,gsl_rng *rx);


void PrintResult(gsl_vector *ytest,gsl_vector *result1,gsl_vector *result2,gsl_vector *result3,char *filename)
{
  size_t i; 
  FILE *fx = fopen(filename,"w");
  for (i=0; i<ytest->size; i++) {
    fprintf(fx," %f  %f  %f  %f\n",
	    ytest->data[i],result1->data[i],result2->data[i],result3->data[i]);
  }
  fclose(fx);    
}


int main (int argc, char *argv[])
{
  static VArgVector in_files;
  static VArgVector x_files;
  static VString out_filename="";
  static VString regressor_filename="";
  static VLong  xdim    = 500;
  static VLong  nensembles = 1000;
  static VLong  nsamples = 50;
  static VShort numtrain = 200;
  static VShort numtest  = 100;
  static VShort xpls    = 10;
  static VLong  seed    = 601931;
  static VShort nproc = 0;
  static VOptionDescRec options[] = {
    { "in", VStringRepn, 0, & in_files, VRequiredOpt, NULL,"Input files" },
    { "out", VStringRepn, 1, & out_filename, VOptionalOpt, NULL,"Output txt-file" },
    { "ntrain", VShortRepn, 1, & numtrain, VOptionalOpt, NULL,"Number of training samples" },
    { "ntest", VShortRepn, 1, & numtest, VOptionalOpt, NULL,"Number of test samples" },
    { "regressor", VStringRepn, 1, & regressor_filename, VRequiredOpt, NULL,"Regressor" },
    { "xxx", VStringRepn, 0, & x_files, VOptionalOpt, NULL,"Supplementary info file(s)" },
    { "dimX", VLongRepn, 1, & xdim, VOptionalOpt, NULL,"Dimension of matrix X" },
    { "npls", VShortRepn, 1, & xpls, VOptionalOpt, NULL,"Number of components in PLS" },
    { "nensembles", VLongRepn, 1, & nensembles, VOptionalOpt, NULL,"Number of ensembles" },
    { "numsamples", VLongRepn, 1, & nsamples, VOptionalOpt, NULL,"Number of samples" },
    { "seed", VLongRepn, 1, & seed, VOptionalOpt, NULL,"Seed for random number generator" },
    { "j",VShortRepn,1,(VPointer) &nproc,VOptionalOpt,NULL,"Number of processors to use, '0' to use all"}
  };
  size_t i,j,sample;

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

   
  /* Various parameters */
  size_t numsubjects = in_files.number;
  size_t ntrain = (size_t)numtrain;
  size_t ntest = (size_t)numtest;
  size_t dimX  = (size_t)xdim;
  size_t npls  = (size_t)xpls; 
  size_t niter = (size_t)nensembles;
  size_t numsamples = (size_t)nsamples;
  
  int *qtrain = (int *)VCalloc(ntrain,sizeof(int));
  int *qtest  = (int *)VCalloc(ntest,sizeof(int));
  if (ntrain < ntest) VError(" ntrain < ntest (%d < %d)",ntrain,ntest);
 
  /* prep output file */
  FILE *fs=NULL;
  if (out_filename != NULL && strcmp(out_filename,"xxx")!=0)  fs = fopen(out_filename,"w");
  VString in_filename = ((VString *) in_files.vector)[0];
  fprintf(stderr,"# %s\n",in_filename);
  fprintf(stderr,"# %s\n",regressor_filename);
  fprintf(stderr,"# ntrain: %lu, ntest: %lu,  dimX: %lu, npls: %lu, nensembles: %ld, numsamples: %lu, seed: %ld\n",
	  ntrain,ntest,dimX,npls,niter,numsamples,seed);
  if (fs != NULL) {
    fprintf(fs,"# %s\n",in_filename);
    fprintf(fs,"# ntrain: %lu, ntest: %lu,  dimX: %lu, npls: %lu,  nensembles: %ld, numsamples: %ld, seed: %ld\n",
	    ntrain,ntest,dimX,npls,niter,numsamples,seed);
  }


  /* read regressor file */
  gsl_vector *regressor = ReadRegressor(regressor_filename);
  if (regressor->size != numsubjects)
    VError(" inconsistent number of inout connectomes versus entries in regressor file, %lu vs %lu",numsubjects,regressor->size);


  /* read X-file, supplementary info */
  VBoolean verbose=FALSE;
  gsl_matrix *X = ReadConfounds(x_files,(int)numsubjects,NULL,verbose);
  if (X != NULL) {
    MatColumnCentering(X);
    Multifit(X,regressor);
  }

  
  /* ini random generator */
  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *rx = gsl_rng_alloc(T);
  gsl_rng_set(rx,(unsigned long int)seed);

   
  /* Select random folds */
  gsl_matrix_int *trainfolds = gsl_matrix_int_calloc(numsamples,ntrain);
  gsl_matrix_int *testfolds = gsl_matrix_int_calloc(numsamples,ntest);
  for (sample=0; sample<numsamples; sample++) { 
    SelectTrainTest(qtrain,qtest,ntrain,ntest,numsubjects,rx);
    for (i=0; i<ntrain; i++) gsl_matrix_int_set(trainfolds,sample,i,qtrain[i]);
    for (i=0; i<ntest; i++) gsl_matrix_int_set(testfolds,sample,i,qtest[i]);
  }

  
  /* read data */
  size_t nedges = GetNumEdges(in_files);
  gsl_matrix_float *Xcorr = gsl_matrix_float_calloc(numsubjects,nedges);
  if (!Xcorr) VError(" err allocating Xcorr");
  
  size_t k0=0,k1=0;
  VReadBundles(in_files,nedges,(size_t)numsubjects,Xcorr,NULL,&k0,&k1);
  if (k0 == 0 && k1 == 0) VError(" err reading input bundles");
  RowNormalizeBundle(Xcorr);
  char *etable = (char *) VCalloc(nedges,sizeof(char));
  for (i=0; i<nedges; i++) etable[i] = 1;

  
  /* main loop */
  gsl_vector *ytrain = gsl_vector_calloc(ntrain);
  gsl_vector *ytest  = gsl_vector_calloc(ntest);  
  gsl_vector *xtrain = gsl_vector_calloc(ntrain);
  gsl_vector *xtest  = gsl_vector_calloc(ntest);
  
  gsl_vector *aresult  = gsl_vector_calloc(ntest);
  gsl_vector *bresult  = gsl_vector_calloc(ntest);
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


  
  /* ini result structs */
  double ya1,ya3,yb1,yb3,yc1,yc3,yx1=0,yx3=0;
  double sa1=0,sa3=0,sb1=0,sb3=0,sc1=0,sc3=0,sx1=0,sx3=0;
  double ssa1=0,ssa6=0,ssb1=0,ssb6=0,ssc1=0,ssc6=0;
  double cx1=0,cx6=0,bx1=0,bx6=0,ax1=0,ax6=0,resa1=0,resa6=0,resb1=0,resb6=0,resc1=0,resc6=0;
  double ymean=0,ysd=0,alpha=1.0;


  /* samples... */
  fprintf(stderr,"# A: Blind machine learning (no supplementary info)\n");
  fprintf(stderr,"# B: Semi-blind machine learning without bias control\n");
  fprintf(stderr,"# C: Semi-blind machine learning with bias control\n");
  fprintf(stderr,"# AX: bias in A\n");
  fprintf(stderr,"# BX: bias in B\n");
  fprintf(stderr,"# CX: bias in C\n");
  fprintf(stderr,"# alpha: parameter used for bias control\n\n");

  if (fs != NULL) {
    fprintf(fs,"# A: Blind machine learning (no supplementary info)\n");
    fprintf(fs,"# B: Semi-blind machine learning without bias control\n");
    fprintf(fs,"# C: Semi-blind machine learning with bias control\n");
    fprintf(fs,"# X: Predictions based solely on supplementaru info, no fMRI\n");
    fprintf(fs,"# AX: bias in A\n");
    fprintf(fs,"# BX: bias in B\n");
    fprintf(fs,"# CX: bias in C\n");
    fprintf(fs,"# alpha: parameter used for bias control\n\n");
  }

  fprintf(stderr,"# %12s %15s  %15s  %15s  %15s  %15s  %15s %15s\n","A","B","C","X","AX","BX","CX","alpha");
  if (fs != NULL) 
    fprintf(fs,"# %12s %15s  %15s  %15s  %15s  %15s  %15s %15s\n","A","B","C","X","AX","BX","CX","alpha");
  
  for (sample=0; sample<numsamples; sample++) {
   
    /* train/test split */
    int *ptrain = gsl_matrix_int_ptr(trainfolds,sample,0);
    int *ptest = gsl_matrix_int_ptr(testfolds,sample,0);
    for (j=0; j<ntrain; j++) ytrain->data[j] = regressor->data[ptrain[j]];
    for (j=0; j<ntest; j++) ytest->data[j] = regressor->data[ptest[j]];

    ymean = gsl_stats_mean(ytrain->data,1,ytrain->size);
    ysd = gsl_stats_sd(ytrain->data,1,ytrain->size);
    gsl_vector_add_constant(ytrain,-ymean);
    gsl_vector_add_constant(ytest,-ymean);
    gsl_vector_scale(ytrain,1.0/ysd);
    gsl_vector_scale(ytest,1.0/ysd);

    
    /* A:  prediction without using supplementary info */
    Sel0(etable,nedges,rx,S0);
    Loop(Xcorr,ytrain,NULL,ptrain,ptest,S0,npls,aresult);
    ya1 = PearsonCorr(ytest,aresult);
    ya3 = CoD(ytest,aresult);
    ya3 = Accuracy(ytest,aresult);
    sa1 += ya1;
    sa3 += ya3;

    
    /* X: prediction based solely on supplementary info */
    if (X==NULL) goto skip;
    XFit(X,ytrain,ytest,xtrain,xtest,ptrain,ptest);
    yx1 = PearsonCorr(ytest,xtest); 
    yx3 = CoD(ytest,xtest);
    yx3 = Accuracy(ytest,xtest);
    sx1 += yx1;
    sx3 += yx3;


    
    /* B: prediction without bias-control */
    alpha=1.0;
    gsl_vector_set_zero(bresult);
    gsl_vector_set_zero(xresult);
    for (j=0; j<mult; j++) {
      int rtcode = Split(ytrain,xtrain,ytest,xtest,sytrain,sxtrain,sytest,sxtest,ptrain,ptest,strain,stest,ktrain,ktest,rx);
      if (rtcode < 0) VError("split: %d",rtcode);
      SelProb(Xcorr,sytrain,sxtrain,sytest,sxtest,strain,stest,etable,S2,rx,(double)alpha,0);
      Loop(Xcorr,ytrain,NULL,ptrain,ptest,S2,npls,xresult);
      gsl_vector_add(bresult,xresult);
    }
    gsl_vector_scale(bresult,1.0/(double)mult);
    yb1 = PearsonCorr(ytest,bresult);
    yb3 = CoD(ytest,bresult);
    yb3 = Accuracy(ytest,bresult);
    sb1 += yb1;
    sb3 += yb3;


    
    /* C: prediction with bias-control, alpha is adjusted */
    double xcorr = PearsonCorr(ytrain,xtrain);
    if (xcorr < 0) VError(" xcorr < 0");
    alpha=1;
    double xdiff=100,up=1.2,down=0.8,xalpha=1;
    int k=0,kiter=30;

    while (fabs(xdiff) > 0.02 && xalpha > 0.01 && k < kiter) {
      gsl_vector_set_zero(cresult);
      gsl_vector_set_zero(xresult);

      for (j=0; j<mult; j++) {
	int rtcode = Split(ytrain,xtrain,ytest,xtest,sytrain,sxtrain,sytest,sxtest,ptrain,ptest,strain,stest,ktrain,ktest,rx);
	if (rtcode < 0) VError("split: %d",rtcode);
	SelProb(Xcorr,sytrain,sxtrain,sytest,sxtest,strain,stest,etable,S2,rx,(double)xalpha,0);
	Loop(Xcorr,ytrain,NULL,ptrain,ptest,S2,npls,xresult);
	gsl_vector_add(cresult,xresult);
      }
      gsl_vector_scale(cresult,1.0/(double)mult);
      double zcorr = PearsonCorr(cresult,xtest);
      xdiff = zcorr-xcorr;
      alpha = xalpha;
      if (k > 1) { up=1.1; down=0.9; }
      if (zcorr > xcorr) xalpha *= down;
      if (zcorr < xcorr) xalpha *= up;
      k++;
    }
    yc1 = PearsonCorr(ytest,cresult);
    yc3 = CoD(ytest,cresult);
    yc3 = Accuracy(ytest,cresult);
    sc1 += yc1;
    sc3 += yc3;

    /* if (sample == 0) PrintResult(ytest,xtest,aresult,cresult,"test.txt"); */
 

    /* influence of supplementary info on final result */
    XResidualFit(xtest,ytest,aresult,&ax1,&ax6,&resa1,&resa6);
    XResidualFit(xtest,ytest,bresult,&bx1,&bx6,&resb1,&resb6);
    XResidualFit(xtest,ytest,cresult,&cx1,&cx6,&resc1,&resc6);
    ssa1 += ax1;
    ssa6 += ax6;
    ssb1 += bx1;
    ssb6 += bx6;
    ssc1 += cx1;
    ssc6 += cx6;

    
 
  skip: ;
    fprintf(stderr," %3lu  %6.3f %7.4f   %6.3f %7.4f   %6.3f %7.4f   %6.3f %7.4f   %6.3f %7.4f   %6.3f %7.4f   %6.3f %7.4f   %6.3f\n",
	    sample,ya1,ya3,yb1,yb3,yc1,yc3,yx1,yx3,ax1,ax6,bx1,bx6,cx1,cx6,alpha);

    if (fs != NULL) {
      fprintf(fs," %3lu  %6.3f %7.4f   %6.3f %7.4f   %6.3f %7.4f   %6.3f %7.4f   %6.3f %7.4f   %6.3f %7.4f   %6.3f %7.4f   %6.3f\n",
	      sample,ya1,ya3,yb1,yb3,yc1,yc3,yx1,yx3,ax1,ax6,bx1,bx6,cx1,cx6,alpha);
    }
  }
  fprintf(stderr,"\n");


  double nx = (double)numsamples;
  sa1 /= nx;
  sa3 /= nx;
  sb1 /= nx;
  sb3 /= nx;
  sc1 /= nx;
  sc3 /= nx;
  sx1 /= nx;
  sx3 /= nx;
  ssb1 /= nx;
  ssb6 /= nx;
  ssc1 /= nx;
  ssc6 /= nx;
  ssa1 /= nx;
  ssa6 /= nx;

  
  fprintf(stderr,"#ave: %6.3f %7.4f   %6.3f %7.4f   %6.3f %7.4f   %6.3f %7.4f   %6.3f %7.4f   %6.3f %7.4f   %6.3f %7.4f\n",
	  sa1,sa3,sb1,sb3,sc1,sc3,sx1,sx3,ssa1,ssa6,ssb1,ssb6,ssc1,ssc6);

  if (fs != NULL) {
    fprintf(fs,"\n");
    fprintf(fs,"#ave: %6.3f %7.4f   %6.3f %7.4f   %6.3f %7.4f   %6.3f %7.4f   %6.3f %7.4f   %6.3f %7.4f   %6.3f %7.4f\n",
	    sa1,sa3,sb1,sb3,sc1,sc3,sx1,sx3,ssa1,ssa6,ssb1,ssb6,ssc1,ssc6);
    fclose(fs);
  }

  fprintf (stderr, "%s: done.\n", argv[0]);
  exit(0);
}
