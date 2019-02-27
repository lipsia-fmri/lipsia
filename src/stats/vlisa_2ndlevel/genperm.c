/* 
** generate permutation table using exchangeability information
**
** G.Lohmann, July 2018
*/
#include <viaio/Vlib.h>
#include <viaio/file.h>
#include <viaio/option.h>
#include <viaio/mu.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

extern gsl_matrix *XRead2ndLevel(VString);
extern void XReadExchange(VString filename,int *exchange,int n);
extern gsl_matrix *PseudoInv(gsl_matrix *A,gsl_matrix *B);



void genperm(long seed,int *exchange,int signswitch,
	     int **permtable,int **signtable,int nimages,int numperm,gsl_vector *contrast)
{
  int i,j,k,m,kmax;


  /* ini permutation and sign-switching tables */
  for (i = 0; i < numperm; i++) {  
    permtable[i] = (int *) VCalloc(nimages,sizeof(int));
    signtable[i] = (int *) VCalloc(nimages,sizeof(int));
    for (j=0; j<nimages; j++) permtable[i][j] = j;
    for (j=0; j<nimages; j++) signtable[i][j] = 1;
  }


  /* ini random number generator */
  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *rx = gsl_rng_alloc(T);
  gsl_rng_set(rx,(unsigned long int)seed);


  /* sign switching only */
  if (signswitch >= 0) {
    for (i=0; i<numperm; i++) {  
      for (j=0; j<nimages; j++) {	
	if (gsl_ran_bernoulli (rx,(double)0.5) == 1) signtable[i][j] = -1;
      }
    }
    return;  /* stop here, no permutations */
  }


  /* get max number of exchange groups */
  kmax=0;
  for (j=0; j<nimages; j++) {
    if (exchange[j] > kmax) kmax = exchange[j]; 
  }


  /* permutations of rows of the design matrix */
  int *base = (int *) VCalloc(nimages,sizeof(int));
  for (i = 0; i < numperm; i++) {  

    /* for each group */
    for (k=0; k<=kmax; k++) {
      m = 0;
      for (j=0; j<nimages; j++) {
	if (exchange[j] == k) {
	  base[m] = j;
	  m++;
	}
      }
      if (m < 1) continue;

      /* random shuffle within group */
      gsl_ran_shuffle (rx, base, m, sizeof (int));

      m = 0;
      for (j=0; j<nimages; j++) {
	if (exchange[j] == k) {
	  permtable[i][j] = base[m];
	  m++;
	}
      }
    }
  }
  VFree(base);
}



/* for zero variance replace permutations with sign switching (only for non-nuisance columns) */
int SignSwitch (gsl_matrix *X,gsl_vector *contrast,int *permflag)
{
  int i,j;
  double s1=0,s2=0,nx=0,u=0,mean=0,var=0,tiny=1.0e-6;

  int signswitch = -1;
  for (j=0; j<X->size2; j++) {
    if (permflag[j] == 0) continue;
    s1 = s2 = nx = 0;
    for (i=0; i<X->size1; i++) { 
      u = gsl_matrix_get(X,i,j);
      s1 += u;
      s2 += u*u;
      nx++;
    }
    mean = s1/nx;
    var = (s2 - nx * mean * mean) / (nx - 1.0);
    if (var < tiny && signswitch >= 0)
      VError(" Implausible design, two columns have zero variance: %d %d",signswitch+1,j+1);
    if (var < tiny) signswitch = j;
  }
  return signswitch;
}
