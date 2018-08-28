#include <viaio/VImage.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <viaio/Vlib.h>
#include "gsl_utils.h"


void prewhite(gsl_matrix_float** pinvX,gsl_matrix_float** invM,gsl_matrix_float* X,int numlags) 
{
  int n = X->size2;
     
  /* preparation for whitening */

  /* X^T */
  gsl_matrix_float* transX = gsl_matrix_float_alloc(X->size2, X->size1);
  gsl_matrix_float_transpose_memcpy(transX,X);
     
  /* identity */
  gsl_matrix_float *R = gsl_matrix_float_alloc(n,n);
  gsl_matrix_float_set_identity(R);
     
  /* pseudo inverse */
  *pinvX = fmat_PseudoInv(transX,NULL);

  /* temporary matrix buffer for binary matrix operations */
  gsl_matrix_float *result;
  gsl_matrix_float *result2;

  result = fmat_x_mat(transX,*pinvX,NULL);
  gsl_matrix_float_sub(R,result);
  gsl_matrix_float_free(result);
     
  /* result buffer */
  gsl_matrix_float *M;

  if(numlags == 1){
    /* create a n-x-n matrix Diag1 where the k=1 sub- and suberdiagonals 
     * are 1 and the remaining entries are 0. 
     * Equivalent matlab code:
     *
     * keep = 1:n;
     * indk1=((keep(2:n)-keep(1:n-1))==1);
     * Diag1=diag(indk1,1)+diag(indk1,-1);
     */
    gsl_matrix_float *Diag1 = gsl_matrix_float_alloc(n,n);
    gsl_matrix_float_set_zero(Diag1);
    /* set superdiagonal to 1 */
    gsl_vector_float_view superview = gsl_matrix_float_superdiagonal(Diag1,1);
    gsl_vector_float_set_all(&superview.vector,1);
    /* set subdiagonal to 1 */
    gsl_vector_float_view subview = gsl_matrix_float_subdiagonal(Diag1,1);
    gsl_vector_float_set_all(&subview.vector,1);

    float M11 = trace(R);
    result = fmat_x_mat(R, Diag1, NULL);
    float M12 = trace(result);
    float M21 = M12/2.0;
    fmat_x_mat(R,Diag1,result);
    result2 = fmat_x_mat(result,R,NULL);
    fmat_x_mat(result2,Diag1,result);
    float M22 = trace(result)/2;
    gsl_matrix_float_free(result);
    gsl_matrix_float_free(result2);

    M = gsl_matrix_float_alloc(2,2);
   

    /* M=[M11 M12; M21 M22] */
    gsl_matrix_float_set(M,0,0,M11);
    gsl_matrix_float_set(M,0,1,M12);
    gsl_matrix_float_set(M,1,0,M21);
    gsl_matrix_float_set(M,1,1,M22);

    gsl_matrix_float_free(Diag1);

  }
  else {

    int i,j;
    M = gsl_matrix_float_alloc(numlags+1,numlags+1);
    gsl_matrix_float_set_zero(M);
    gsl_matrix_float* Di = gsl_matrix_float_alloc(n,n);
    gsl_matrix_float* Dj = gsl_matrix_float_alloc(n,n);

    /* This double loop can be supported by a look-up-table for 
     * the diagonal matrixes. */
    for(i=0;i<numlags+1;i++){
      gsl_matrix_float_set_zero(Di);
      /* set maindiagonal to 1 */
      if(i==0){
	gsl_vector_float_view mdiag = gsl_matrix_float_diagonal(Di);
	gsl_vector_float_set_all(&mdiag.vector,1);
      }
      /* set i. super- and subdiagonal to 1 */ 
      else {
	/* set superdiagonal to 1 */
	gsl_vector_float_view superview = gsl_matrix_float_superdiagonal(Di,i);
	gsl_vector_float_set_all(&superview.vector,1);
	/* set subdiagonal to 1 */
	gsl_vector_float_view subview = gsl_matrix_float_subdiagonal(Di,i);
	gsl_vector_float_set_all(&subview.vector,1);
      }
      for(j=0;j<numlags+1;j++){
	gsl_matrix_float_set_zero(Dj);
	/* set maindiagonal to 1 */
	if(j==0){
	  gsl_vector_float_view mdiag = gsl_matrix_float_diagonal(Dj);
	  gsl_vector_float_set_all(&mdiag.vector,1);
	}
	/* set j. super- and subdiagonal to 1 */ 
	else {
	  /* set superdiagonal to 1 */
	  gsl_vector_float_view superview = gsl_matrix_float_superdiagonal(Dj,j);
	  gsl_vector_float_set_all(&superview.vector,1);
	  /* set subdiagonal to 1 */
	  gsl_vector_float_view subview = gsl_matrix_float_subdiagonal(Dj,j);
	  gsl_vector_float_set_all(&subview.vector,1);
	}   

	result = fmat_x_mat(R,Di,NULL);
	result2 = fmat_x_mat(result,R,NULL);
	fmat_x_mat(result2,Dj,result);
	gsl_matrix_float_free(result2);

	gsl_matrix_float_set(M,i,j,trace(result)/(1+(i>0)));
	gsl_matrix_float_free(result);

      }
    }

    gsl_matrix_float_free(Di);
    gsl_matrix_float_free(Dj);
  }    

  /* calculation of the return values */
  *invM = fInv(M,NULL);
  /* pinvX has already been initialized */     

  /* free memory */
  gsl_matrix_float_free(transX);
  gsl_matrix_float_free(R);
  gsl_matrix_float_free(M);
}

