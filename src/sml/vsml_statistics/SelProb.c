#include <viaio/Vlib.h>

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>


/* random selection selection */
size_t Sel0(char *table,size_t nedges,gsl_rng *rx,gsl_matrix_long *S)
{
  size_t iter,i,j;
  size_t niter = S->size1;
  size_t dimX = S->size2;

  long *ibase  = (long *) VCalloc(nedges,sizeof(long));
  long *select = (long *) VCalloc(dimX,sizeof(long));
  size_t n=0;
  for (i=0; i<nedges; i++) {
    if (table[i] > 0) {
      ibase[n] = i;
      n++;
    }
  }

  for (iter=0; iter<niter; iter++) {
    gsl_ran_sample(rx,select,dimX,ibase,n,sizeof(long));
    for (j=0; j<dimX; j++) {
      gsl_matrix_long_set(S,iter,j,select[j]);
    }
  }
  VFree(select);
  VFree(ibase);
  return n;
}


/* univariate correlation with regressor */
void Corrmap(gsl_matrix_float *X,gsl_vector *ytrain,int *train,double *W)
{
  size_t k,ntrain=ytrain->size;
  size_t nedges=X->size2;
 
#pragma omp parallel for schedule(guided)
  for (k=0; k<nedges; k++) {
    size_t i;
    double z=0;
    gsl_vector *xdata = gsl_vector_calloc(ntrain);
    for (i=0; i<ntrain; i++) {
      xdata->data[i] = (double)gsl_matrix_float_get(X,train[i],k);
    }
    z = gsl_stats_correlation(xdata->data,1,ytrain->data,1,ntrain);
    W[k] = z;
    gsl_vector_free(xdata);
  }
}



/* sigmoid */
void Sigmoid(double *W,size_t nedges,double alpha)
{
  size_t i;
  double x;

  double mean = gsl_stats_mean(W,1,nedges);
  double sd = gsl_stats_sd_m(W,1,nedges,mean);
  
  for (i=0; i<nedges; i++) {
    x = (W[i]-mean)/sd;
    W[i] = 1.0/(1.0+exp(-alpha*x));
  }
}



/* Probabilistic sampling */
double SelProb(gsl_matrix_float *Xcorr,
	       gsl_vector *ytrain,gsl_vector *xtrain,gsl_vector *ytest,gsl_vector *xtest,int *train,int *test,
	       char *etable,gsl_matrix_long *S,gsl_rng *rx,double alpha,int seltype)
{
  size_t iter,i,j;
  size_t nedges = Xcorr->size2;
  size_t niter = S->size1;
  size_t dimX = S->size2;

  /* corr maps */
  double *X1 = (double *) VCalloc(nedges,sizeof(double *));
  double *Y0 = (double *) VCalloc(nedges,sizeof(double *));
  double *W = (double *) VCalloc(nedges,sizeof(double *));
  
  Corrmap(Xcorr,xtest,test,X1);
  Corrmap(Xcorr,ytrain,train,Y0);

  
  /* product */
  for (i=0; i<nedges; i++) {
    W[i] = X1[i]*Y0[i];
  }
  VFree(X1);
  VFree(Y0);


  /* scaling */
  Sigmoid(W,nedges,alpha);

  
  /* probabilistic sampling */
  for (i=0; i<nedges; i++) if (W[i] < 0) W[i] = 0;
  for (i=0; i<nedges; i++) if (etable[i] == 0) W[i] = 0;
  gsl_ran_discrete_t *prob = gsl_ran_discrete_preproc(nedges,W);
  VFree(W);

  for (iter=0; iter<niter; iter++) {
    for (i=0; i<dimX; i++) {
      j = (long)gsl_ran_discrete(rx,prob);
      gsl_matrix_long_set(S,iter,i,j);
    }
  }
  gsl_ran_discrete_free(prob);

  return 0;
}



int Split(gsl_vector *ytrain,gsl_vector *xtrain,gsl_vector *ytest,gsl_vector *xtest,
	  gsl_vector *sytrain,gsl_vector *sxtrain,gsl_vector *sytest,gsl_vector *sxtest,
	  int *train,int *test,
	  int *strain,int *stest,size_t ktrain,size_t ktest,gsl_rng *rx)
{
  size_t i,j,k;
  size_t ntrain = ytrain->size;
  size_t ntest = xtest->size;
  if (ktest < ntest) {
    fprintf(stderr," ktest: %lu,  ntest: %lu\n",ktest,ntest);
    return -1;
  }

  k = ktest-ntest;
  if (k > ntrain) return -2;
  int *sel = (int *)VCalloc(k,sizeof(int));
  
  int *ibase = (int *)VCalloc(ntrain,sizeof(int));
  for (i=0; i<ntrain; i++) ibase[i]=i;

  int *tab = (int *)VCalloc(ntrain,sizeof(int));
  for (i=0; i<ntrain; i++) tab[i] = 0;
  
  gsl_ran_choose(rx,sel,k,ibase,ntrain,sizeof(int));
  for (i=0; i<k; i++) tab[sel[i]] = 1;

  j = 0;
  for (i=0; i<ntrain; i++) {
    if (tab[i] == 0) {
      strain[j] = train[i];
      sytrain->data[j] = ytrain->data[i];
      sxtrain->data[j] = xtrain->data[i];
      j++;
    }
  }
  if (j != ktrain) return -3;

  for (i=0; i<ntest; i++)  {
    stest[i] = test[i];
    sxtest->data[i] = xtest->data[i];
    sytest->data[i] = ytest->data[i];
  }
  for (i=0; i<k; i++) {
    stest[i+ntest] = train[sel[i]];
    sxtest->data[i+ntest] = xtrain->data[sel[i]];    
    sytest->data[i+ntest] = ytrain->data[sel[i]];
  }
  VFree(ibase);
  VFree(sel);
  VFree(tab);
  return 0;
}


