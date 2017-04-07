
#include <viaio/Vlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_combination.h>

/* binomial coefficient, n choose k */
double binom (int n,int k)
{
  double prod,ix,nx,kx;
  int i;
  nx = (double)n;
  kx = (double)k;

  prod = ix = 1;
  for (i=1; i<=k; i++) {
    prod *= (nx-kx+ix)/ix;
    ix++;
  }
  return prod;
}

char **randomCombinations(long maxiter,int n,unsigned long int seed)
{
  size_t iter,i,j,k;
  size_t znum,maxnum=90000;

  /* ini random */
  fprintf(stderr," random shuffle\n");
  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *rx = gsl_rng_alloc(T);
  gsl_rng_set(rx,(unsigned long int)seed);


  /* ini data structs */
  int  *src = (int *) calloc(n,sizeof(int));
  char *tmp = (char *)calloc(n,sizeof(char));
  for (i=0; i<n; i++) src[i] = i;

  char  **table = (char **)calloc(maxiter,sizeof(char *));
  for (i=0; i<maxiter; i++)
    table[i] = (char *)calloc(n,sizeof(char));
 

  /* generate permutation table */
  for (iter=0; iter<maxiter; iter++) {
    znum = 0;

  skip:
    if (znum >= maxnum) VError(" znum %d",znum);
    k = n/2;
    if (n%2 == 1 && iter%2 == 0) k++;  /* for odd number of images */
    gsl_ran_shuffle (rx,src,n,sizeof(int));
    
    memset(tmp,'0',n*sizeof(char));
    for (j=0; j<k; j++) tmp[src[j]] = '1';

    /* check if already there */
    for (j=0; j<iter; j++) { 
      if (memcmp(tmp,table[j],(size_t)n) == 0 && znum < maxnum) {
	znum++;
	goto skip;
      }
    }
    memcpy(table[iter],tmp,(size_t)n);
  }
  return table;
}

char **XXrandomCombinations(long maxiter,int n,unsigned long int seed)
{
  size_t iter,i,j;
  double p=0.5;

  /* ini random */
  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *rx = gsl_rng_alloc(T);
  gsl_rng_set(rx,(unsigned long int)seed);

  /* ini data structs */
  char  *tmp = (char *)calloc(n,sizeof(char));
  char  **table = (char **)calloc(maxiter,sizeof(char *));
  for (i=0; i<maxiter; i++)
    table[i] = (char *)calloc(n,sizeof(char));

  /* generate permutation table */
  iter = 0;
  while (iter < maxiter) {
    memset(tmp,'0',n*sizeof(char));

    int ic=0;
    for (j=0; j<n; j++) {
      if (gsl_ran_bernoulli(rx,p) > 0) {
	tmp[j] = '1';
	ic++;
      }
    }
    for (j=0; j<iter; j++) {  /* check if already there */
      if (memcmp(tmp,table[j],n) == 0) goto skip;
    }
    memcpy(table[iter],tmp,(size_t)n);
    iter++;
  skip: ;
  }
  return table;
}


char **allCombinations(long maxiter,int n,int iskip)
{
  size_t iter,i,j,num,b;

  /* ini data structs */
  char  *tmp = (char *)calloc(n,sizeof(char));  
  char  **table = (char **)calloc(maxiter,sizeof(char *));
  for (i=0; i<maxiter; i++)
    table[i] = (char *)calloc(n,sizeof(char));

  size_t xmax = (size_t) gsl_pow_int(2.0,(int)n);

  iter = 0;
  for (i=0; i<xmax; i+=iskip) {
    /* decimal to binary */
    memset(tmp,'0',n*sizeof(char));
    j=0;
    num = i;
    while(num > 0) {
      b=num%2;
      if (b == 1) tmp[j] = '1';
      num=num/2;
      j++;
    }
    /* fprintf(stderr," iter: %5d  %s\n",iter,tmp); */
    memcpy(table[iter],tmp,(size_t)n);
    iter++;
    if (iter >= maxiter) return table;
  }
  return table;
}
