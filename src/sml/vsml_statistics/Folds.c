#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_permutation.h>

extern double CoD(gsl_vector *regressor,gsl_vector *result);
extern double Check(gsl_vector *ytest,gsl_vector *result,int type);

/* Ini folds, cross-validation */
void IniFolds(int **train,int **test,int ntrain,int ntest,int numfolds,gsl_rng *rx)
{
  int i,j,k,fold;
  int numsubjects = ntrain+ntest;

  int *base = (int *)VCalloc(numsubjects,sizeof(int));
  for (i=0; i<numsubjects; i++) base[i] = i;
  gsl_ran_shuffle(rx,base,(size_t)numsubjects,sizeof(int));

  j=0;
  for (fold=0; fold<numfolds; fold++) {
    for (i=0; i<ntest; i++) {
      test[fold][i] = base[j];
      j++;
    }
  }

  for (fold=0; fold<numfolds; fold++) {
    j=0;
    for (k=0; k<numfolds; k++) {
      if (fold == k) continue;
      for (i=0; i<ntest; i++) {
	train[fold][j] = test[k][i];
	j++;
      }
    }
  }
  VFree(base);  
}


/* Ini folds, random selection of non-overlapping train-test sets */ 
void SelectTrainTest(int *train,int *test,size_t ntrain,size_t ntest,size_t numsubjects,gsl_rng *rx)
{
  size_t i,j;

  int *ibase = (int *)VCalloc(numsubjects,sizeof(int));
  for (i=0; i<numsubjects; i++) ibase[i] = i;

  int *table = (int *)VCalloc(numsubjects,sizeof(int));
  for (i=0; i<numsubjects; i++) table[i] = 0;
  
  gsl_ran_choose(rx,train,ntrain,ibase,numsubjects,sizeof(int));
  for (i=0; i<ntrain; i++) table[train[i]] = 1;

  j=0;
  for (i=0; i<numsubjects; i++) {
    if (table[i] == 1) continue;
    ibase[j] = i;
    j++;
  }

  if (j < ntest) VError(" SelectTrainTest, j < ntest, j: %lu, ntest: %lu",j,ntest);
  gsl_ran_choose(rx,test,ntest,ibase,j,sizeof(int));

  VFree(ibase);
  VFree(table);
}



/*
** Ini folds, random selection of non-overlapping train-test sets
** train, test sets from two different segments, train segment=2/3.
*/
void TrainTestSets(int *trainset,int *testset,size_t ntrainset,size_t ntestset,size_t numsubjects,
		   int *train,int *test,size_t ntrain,size_t ntest,gsl_rng *rx)
{
  size_t i,j;

  int *ibase = (int *)VCalloc(numsubjects,sizeof(int));
  for (i=0; i<numsubjects; i++) ibase[i] = i;
  int *table = (int *)VCalloc(numsubjects,sizeof(int));
  for (i=0; i<numsubjects; i++) table[i] = 0;

  gsl_ran_choose(rx,trainset,ntrainset,ibase,numsubjects,sizeof(int));
  for (i=0; i<ntrainset; i++) table[trainset[i]] = 1;

  j=0;
  for (i=0; i<numsubjects; i++) {
    if (table[i] == 0) {
      testset[j] = i;
      j++;
    }
  }
  if (j != ntestset) VError(" fold testset, %lu %lu",j,ntestset);
  VFree(table);
  VFree(ibase);
}



void PrintFolds(gsl_matrix_int *trainfolds, gsl_matrix_int *testfolds)
{
  int i,sample;
  size_t len=512;
  size_t numsamples = trainfolds->size1;
  int ntrain = trainfolds->size2;
  int ntest = testfolds->size2;
  
  fprintf(stderr," ntrain: %d, ntest: %d\n",ntrain,ntest);

  char *filename = (char *)VCalloc(len,sizeof(char));
  memset(filename,0,len);
  sprintf(filename,"folds_%d_%d.txt",ntrain,ntest);
  fprintf(stderr," %s\n",filename);
  				  
  FILE *fx = fopen(filename,"w");
  
  for (sample=0; sample<numsamples; sample++) {
    fprintf(fx,"#----------------------------------------\n");
    fprintf(fx,"# train %d\n",sample);
    
    for (i=0; i<ntrain; i++) {
      fprintf(fx," %d\n",gsl_matrix_int_get(trainfolds,sample,i));
    }
    fprintf(fx,"# test %d\n",sample);
    for (i=0; i<ntest; i++) {
      fprintf(fx," %d\n",gsl_matrix_int_get(testfolds,sample,i));
    }
  }
  fclose(fx);

  exit(0);
}
