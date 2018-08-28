#include <viaio/VImage.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <viaio/Vlib.h>
#include "gsl_utils.h"

void whitecov(VImage rho_vol, gsl_matrix_float *Y,
	      gsl_matrix_float* invM, gsl_matrix_float* pinvX,
	      gsl_matrix_float* X, int numlags, int slice)
{
    
  /* set up variables */
  int n = Y->size1;
  /* these values represent the matlab interval 2:n */
  int k1 = 1;
  int k2 = n-1;
  /* loop counter */
  int i,j;
    
  /* matrix buffers */
  gsl_matrix_float *result,*buffer;
    
  /* : X^T */
  gsl_matrix_float* transX = gsl_matrix_float_alloc(X->size2,X->size1);
  gsl_matrix_float_transpose_memcpy(transX,X);

  /* : betahat_ls = pinvX * Y; */
  gsl_matrix_float* betahat_ls = fmat_x_mat(pinvX, Y, NULL);

  /* gsl_blas_sgemm(CblasNoTrans,CblasNoTrans,1.0,pinvX,Y,0.0,betahat_ls); */


  /* : resid = Y - X' * betahat_ls; */
  result = fmat_x_mat(transX, betahat_ls, NULL);
  gsl_matrix_float_sub(Y,result);
  gsl_matrix_float_free(result);

  /* rename */
  gsl_matrix_float* resid = Y;          
    
  if (numlags == 1) {

    /* : Cov0=sum(resid.*resid,1); */
    result = gsl_matrix_float_alloc(resid->size1, resid->size2);
    gsl_matrix_float_memcpy (result, resid);
    gsl_matrix_float_mul_elements (result, result);
    gsl_vector_float* Cov0 = fsum(result, 1,NULL);
    gsl_matrix_float_free(result);
        
    /* : Cov1=sum(resid(k1,:).*resid(k1-1,:),1); */
    gsl_matrix_float_view sub1 =
      gsl_matrix_float_submatrix(resid, k1,0,k2,resid->size2);
    gsl_matrix_float_view sub2 = 
      gsl_matrix_float_submatrix(resid, k1-1, 0,k2,resid->size2);
    result = gsl_matrix_float_alloc(k2, resid->size2);
    gsl_matrix_float_memcpy (result, &sub1.matrix);
    gsl_matrix_float_mul_elements (result, &sub2.matrix);               
    gsl_vector_float* Cov1 = fsum(result,1,NULL);
    gsl_matrix_float_free(result);
        
    /* : Covadj=invM*[Cov0; Cov1]; */
    buffer = gsl_matrix_float_alloc(2,Cov0->size);
    float* pBuffer = buffer->data;
    float* pCov0 = Cov0->data;
    float* pCov1 = Cov1->data;
    for(i=0;i<Cov0->size;i++) {
      *pBuffer++ = *pCov0++;            
    }
    for(i=0;i<Cov1->size;i++) {
      *pBuffer++ = *pCov1++;
    }
    gsl_matrix_float* Covadj = fmat_x_mat(invM, buffer, NULL);
    gsl_matrix_float_free(buffer);        

    /* :  rho_vol(:,slice,1)=(Covadj(2,:)./ ...
     * (Covadj(1,:)+(Covadj(1,:)<=0)).*(Covadj(1,:)>0))'; */
    int cols = Covadj->size2;
    float* p = Covadj->data;
    float val;        
    for(i=0;i<cols;i++) {
      val = *(p+cols) / ((*p <= 0) ? *p +1 : *p) * (*p > 0);
      VPixel(rho_vol, i, slice, 0, VFloat) = val;
      p++;
    }

    gsl_matrix_float_free(Covadj);
    gsl_vector_float_free(Cov0);
    gsl_vector_float_free(Cov1);
  }

  else {
    /* :
     * for lag=0:numlags
     *    Cov(lag+1,:)=sum(resid(1:(n-lag),:).*resid((lag+1):n,:));
     * end
     */
    int lag;
    gsl_matrix_float* Cov = gsl_matrix_float_alloc(numlags+1,resid->size2);        
    gsl_vector_float* v;
    for(lag=0;lag<=numlags;lag++) {
      gsl_matrix_float_view sub1 = 
	gsl_matrix_float_submatrix(resid,0,0,n-lag,resid->size2);
      gsl_matrix_float_view sub2 = 
	gsl_matrix_float_submatrix(resid,lag,0,n-lag,resid->size2 );
      result = gsl_matrix_float_alloc(n-lag, resid->size2);
      gsl_matrix_float_memcpy (result, &sub1.matrix);
      gsl_matrix_float_mul_elements (result, &sub2.matrix);
      v = fsum(result, 1,NULL);
      gsl_matrix_float_free(result);
                        
      memcpy((float*)(Cov->data+(lag*Cov->size2)),
	     (float*)v->data,
	     sizeof(float)*v->size);
            
      gsl_vector_float_free(v);
    }
                
    /* :  Covadj=invM*Cov; */
    gsl_matrix_float* Covadj = fmat_x_mat(invM, Cov, NULL);
        
    /* : rho_vol(:,slice,:)= ( Covadj(0:(numlags+1),:) ...
     *  .*( ones(numlags,1)*((Covadj(1,:)>0)./ ...
     * (Covadj(1,:)+(Covadj(1,:)<=0)))) )'; */
    gsl_matrix_float* line = gsl_matrix_float_alloc(1,Covadj->size2);
    float* pC = Covadj->data;
    float* pB = line->data; 
    for(i=0;i<line->size2;i++) {
      *pB = (float)(*pC > 0) / (float)(*pC+(*pC<=0));            
      pB++;
      pC++;
    }
        
    buffer = gsl_matrix_float_alloc(numlags,1);
    gsl_matrix_float_set_all(buffer,1);
    result = fmat_x_mat(buffer, line,NULL);
        
    gsl_matrix_float_view sub = 
      gsl_matrix_float_submatrix(Covadj,1,0,numlags,Covadj->size2);        
    gsl_matrix_float_free(buffer);
    buffer = gsl_matrix_float_alloc(result->size1, result->size2);
    gsl_matrix_float_memcpy (buffer, &sub.matrix);
    gsl_matrix_float_mul_elements (buffer, result);

    float* p = buffer->data;
    for(i=0;i<buffer->size1;i++) {
      for(j=0;j<buffer->size2;j++) {
	VPixel(rho_vol,j,slice,i, VFloat) = (float)*p++;
      }
    }
        
    gsl_matrix_float_free(buffer);
    gsl_matrix_float_free(result);
    gsl_matrix_float_free(line);
    gsl_matrix_float_free(Covadj);
    gsl_matrix_float_free(Cov);
  }
    
  gsl_matrix_float_free(betahat_ls);
  gsl_matrix_float_free(transX);
}
