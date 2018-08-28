#include <viaio/Vlib.h>
#include <viaio/VImage.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include "gsl_utils.h"

#define ROUND(x) (int)(((x) >= 0) ? ((x) + 0.5) : ((x) - 0.5))


void dfs(VFloat** Dfs,gsl_matrix_float* pinvX,gsl_matrix_float* X,gsl_matrix_float* con,int numlags) 
{
  int numcon = 1;
  int n = X->size2;
  /* these vals represent the matlab interval 2:n */
  int k1 = 1;
  int k2 = n-1;
  int i;
    
    
  /* X^T */
  gsl_matrix_float* transX = gsl_matrix_float_alloc(X->size2, X->size1);
  gsl_matrix_float_transpose_memcpy(transX,X);
        
  gsl_vector_float* CorX2 = gsl_vector_float_alloc(numcon);
  gsl_vector_float_set_zero(CorX2);
    
  /* cpinvX=contrast^T*pinvX */
  gsl_matrix_float* cpinvX = fmat_x_mat(con, pinvX, NULL);
    
  /* CovX0=cpinvX*cpinvX' */
  gsl_matrix_float* CovX0 = fmat_x_matT(cpinvX,cpinvX,NULL);
    
  if(numlags == 1) {
    /* CovX1=cpinvX(:,k1)*cpinvX(:,k1-1)' */
    gsl_matrix_float_view sub1 = 
      gsl_matrix_float_submatrix(cpinvX,0,k1,cpinvX->size1,k2);
    gsl_matrix_float_view sub2 = 
      gsl_matrix_float_submatrix(cpinvX,0,k1-1,cpinvX->size1,k2); 
    gsl_matrix_float* CovX1 = fmat_x_matT(&sub1.matrix, &sub2.matrix, NULL);
   
    gsl_vector_float_view CovX1diag = gsl_matrix_float_diagonal(CovX1);
    gsl_vector_float_view CovX0diag = gsl_matrix_float_diagonal(CovX0);
  
    gsl_vector_float_div (&CovX1diag.vector, &CovX0diag.vector);
    gsl_vector_float_mul (&CovX1diag.vector, &CovX1diag.vector);

    gsl_vector_float_memcpy (CorX2, &CovX1diag.vector);
  }
  else {
    int lag;
    gsl_matrix_float_view sub1, sub2;
    gsl_matrix_float* CovX1 = gsl_matrix_float_alloc(cpinvX->size1,cpinvX->size1);
    for(lag=0;lag<numlags;lag++){
      sub1 = gsl_matrix_float_submatrix(cpinvX,0,0,cpinvX->size1,k2-lag);
      sub2 = gsl_matrix_float_submatrix(cpinvX,0,lag+1,cpinvX->size1,k2-lag);
      fmat_x_matT(&sub1.matrix, &sub2.matrix,CovX1);
 
      gsl_vector_float_view CovX1diag = gsl_matrix_float_diagonal(CovX1);
      gsl_vector_float_view CovX0diag = gsl_matrix_float_diagonal(CovX0);

      gsl_vector_float_div (&CovX1diag.vector, &CovX0diag.vector);
      gsl_vector_float_mul (&CovX1diag.vector, &CovX1diag.vector);

      gsl_vector_float_add(CorX2,&CovX1diag.vector);
    }
    gsl_matrix_float_free(CovX1);
  }

  float dfresid  = n-rank(transX);
  float dfmin = dfresid;

  for (i=0;i<numcon;i++) {
    (*Dfs)[i] = (VFloat)ROUND(dfresid / (1 + 2*gsl_vector_float_get(CorX2,i)));
    if ((*Dfs)[i]<dfmin) dfmin = (*Dfs)[i];
  }
  (*Dfs)[numcon] = dfresid;
    
  /*
    if(dfmin < 100)
    fprintf(stderr, "  warning, dfmin= %f\n",dfmin);
  */
	
  gsl_matrix_float_free(cpinvX);
  gsl_matrix_float_free(transX);
  gsl_matrix_float_free(CovX0);
  gsl_vector_float_free(CorX2);
}
