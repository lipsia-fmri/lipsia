/* 
** pseudoinverse
**
** G.Lohmann, July 2004
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "viaio/Vlib.h"

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

extern void printmat(gsl_matrix *R,char *str);
extern void printvec(gsl_vector *x,char *str);


/*
**    B = A^-1 = A^+
*/
gsl_matrix *PseudoInv(gsl_matrix *A,gsl_matrix *B)
{
  int j,k,l,n,m;
  double u,x;
  gsl_matrix *U=NULL,*V=NULL;
  gsl_vector *w=NULL;

  m = A->size1;
  n = A->size2;

  if (B == NULL) {
    B = gsl_matrix_calloc (n,m);
  }

  U = gsl_matrix_calloc (m, n);
  V = gsl_matrix_calloc (n, n);
  w = gsl_vector_calloc (n);


  /* singular value decomposition */
  gsl_matrix_memcpy (U,A);
  gsl_linalg_SV_decomp_jacobi(U,V,w);
  gsl_matrix_set_zero(B);

  
  k=0;
  for (j=0; j<n; j++) {
    if (fabs(w->data[j]) > 1.0e-6) k++;
  }
  if (k < 1) VError(" Design matrix is singular\n");
  

  for (k=0; k<n; k++) {
    for (l=0; l<m; l++) {
      for (j=0; j<n; j++) {
	u = gsl_vector_get(w,j);
        if (fabs(u) > 1.0e-6) {
          x = gsl_matrix_get(B,k,l);
	  x += gsl_matrix_get(V,k,j)*gsl_matrix_get(U,l,j)/u;
          gsl_matrix_set(B,k,l,x);
        }
      }
    }
  }
  gsl_matrix_free(U);
  gsl_matrix_free(V);
  gsl_vector_free(w);
  return B;
}
