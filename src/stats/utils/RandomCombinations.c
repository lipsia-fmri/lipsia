/*
** random swap functions
**
** M. Kuhlmann, MPI-KYB, Sept 2015
**/

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

double factorial(int n) {
  int i;
  double ix = 1;
  double prod = 1;
  for (i=1; i<=n; i++) {
    prod *= ix;
    ix++;
  }
  return prod;
}

char **enumerate_swaps(int n)
{
  size_t i,j,num,b;

  long maxiter = (size_t) gsl_pow_int(2.0,(int)n);
  
  char **swaps = (char **)VCalloc(maxiter,sizeof(char *));
  for(i = 0; i < maxiter; i++)
    swaps[i] = (char *)VCalloc(n,sizeof(char));

  /* generate swap table */
  for (i=0; i<maxiter; i++) {
    /* decimal to binary */
    j = 0;
    num = i;
    while(num > 0)
    {
      b=num%2;
      swaps[i][j] = b;
      num=num/2;
      j++;
    }
  }
  
  return swaps;
}

char **randomSwaps(long maxiter,int n,unsigned long int seed)
{
  size_t i,j;
  
  char **swaps = (char **)VCalloc(maxiter,sizeof(char *));
  for(i = 0; i < maxiter; i++)
    swaps[i] = (char *)VCalloc(n,sizeof(char));
  
  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *rx = gsl_rng_alloc(T);
  gsl_rng_set(rx,(unsigned long int)seed);
  
  for(i = 0; i < maxiter; i++)
    for(j = 0; j < n; j++)
      swaps[i][j] = gsl_ran_bernoulli(rx, 0.5);
  
  return swaps;
}

char **randomSwaps_noduplicates(long maxiter,int n,unsigned long int seed)
{
  size_t iter,i,j;
  
  /* if number of maxiter exceeds number of all possible combinations */
  if(maxiter > (size_t) gsl_pow_int(2.0,(int)n))
    VError("Error: Number of iterations larger than possible permutations!");
  
  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *rx = gsl_rng_alloc(T);
  gsl_rng_set(rx,(unsigned long int)seed);
    
  char **swaps = (char **)VCalloc(maxiter,sizeof(char *));
  for(i = 0; i < maxiter; i++)
    swaps[i] = (char *)VCalloc(n,sizeof(char));
    
  /* generate swap table */
  iter = 0;
  while (iter < maxiter) {
    for (j=0; j<n; j++)
      swaps[iter][j] = gsl_ran_bernoulli(rx,0.5);
    /* check if already there */
    /* decrementing iter will result in overwriting the swaps */
    /* since there are no duplicates in can only be decremented once */
    for (j=0; j<iter; j++)
      if (memcmp(swaps[iter],swaps[j],n*sizeof(char)) == 0)
	iter--;
    iter++;
  }
  
  return swaps;
}

/*
 * shuffle blocks; needed for permutation_test
 */
int **randomBlocks(long int maxiter, int nimages, int block_size, unsigned long int seed)
{
  int nblocks = nimages/block_size;
  int i,j,iter;
  
  /* ini with real order */
  int **table = (int **)VCalloc(maxiter,sizeof(int*));
  for(i = 0; i < maxiter; i++)
  {
    table[i] = (int *)VCalloc(2*nimages,sizeof(int));
    for(j = 0; j < nimages; j++)
    {
      table[i][2*j] = j;
      table[i][2*j+1] = j+nimages;
    }
  }
  
  /* each row of table has order: 0, nimages, 1, nimages+1, ... */
  /* shuffle blocks of size 2*block_size */
  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *rx = gsl_rng_alloc(T);
  gsl_rng_set(rx,(unsigned long int)seed);
  
  for (iter=0; iter<maxiter; iter++)
    for (j=0; j<nblocks; j++)
      gsl_ran_shuffle(rx, (void *)&table[iter][j*2*block_size], 2*block_size, sizeof(int));
  return table;
}
/*
 * shuffle blocks; needed for permutation_test
 * maxiter is the number of random permutations
 * nimages is the number of images in one condition (2*nimages in total)
 * block_size is the size of the permutation blocks
 */
int **randomBlocks_noduplicates(long int maxiter, int nimages, int block_size, unsigned long int seed)
{
  int nblocks = nimages/block_size;
  int i,j,iter;
  
  /* ini with real order */
  int **table = (int **)VCalloc(maxiter,sizeof(int*));
  for(i = 0; i < maxiter; i++)
  {
    table[i] = (int *)VCalloc(2*nimages,sizeof(int));
    for(j = 0; j < nimages; j++)
    {
      table[i][2*j] = j;
      table[i][2*j+1] = j+nimages;
    }
  }
  
  /* each row of table has order: 0, nimages, 1, nimages+1, ... */
  /* shuffle blocks of size 2*block_size */
  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *rx = gsl_rng_alloc(T);
  gsl_rng_set(rx,(unsigned long int)seed);
  
  iter = 0;
  while (iter < maxiter) {
    /* shuffle */
    for (j = 0; j<nblocks; j++)
      gsl_ran_shuffle(rx, (void *)&table[iter][j*2*block_size], 2*block_size, sizeof(int));
    
    /* sort for one side */
    int k, l, m;
    int index = 0;
    int temp = 0;
    for (m = 0; m < nblocks; m++) {
      for (k = 0; k < block_size-1; k++) {
	index = k;
	for (l = k+1; l < block_size; l++) {
	  if (table[iter][m*2*block_size+2*index] < table[iter][m*2*block_size+2*l])
	    index = l;
	}
	/* swap pairs 'k' and 'index' */
	temp = table[iter][m*2*block_size+2*k];
	table[iter][m*2*block_size+2*k] = table[iter][m*2*block_size+2*index];
	table[iter][m*2*block_size+2*index] = temp;
	temp = table[iter][m*2*block_size+2*k+1];
	table[iter][m*2*block_size+2*k+1] = table[iter][m*2*block_size+2*index+1];
	table[iter][m*2*block_size+2*index+1] = temp;
      }
    }
    
    /*
     * decrementing iter will result in overwriting the last swaps row
     * since there are no duplicates it can only be decremented once
     */
    for(j=0; j<iter; j++)
      if (memcmp(table[iter],table[j],2*nimages*sizeof(int)) == 0)
	iter--;
    iter++;
  }
  return table;
}