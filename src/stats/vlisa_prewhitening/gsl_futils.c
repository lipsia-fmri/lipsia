/* 
** gsl routines for float data type
**
** G.Lohmann, July 2004
*/
#include <viaio/Vlib.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sort_vector.h>
#include "gsl_utils.h" 


/*
** y = Ax
*/
gsl_vector_float *fmat_x_vector(gsl_matrix_float *A, gsl_vector_float *x, gsl_vector_float *y)
{
  if (y == NULL) {
    y = gsl_vector_float_alloc (A->size1);
  }

  gsl_blas_sgemv(CblasNoTrans, 1.0, A, x, 0.0, y);

  return y;
}


/*
**    C = A x B^T
*/
gsl_matrix_float *fmat_x_matT(gsl_matrix_float *A,gsl_matrix_float *B,gsl_matrix_float *C)
{
  if (C == NULL) {
    C = gsl_matrix_float_alloc( A->size1, B->size1 );
  }

  gsl_blas_sgemm( CblasNoTrans, CblasTrans, 1.0, A, B, 0.0, C );

  return C;
}


/*
**    C = A^T x B
*/
gsl_matrix_float *fmatT_x_mat(gsl_matrix_float *A,gsl_matrix_float *B,gsl_matrix_float *C)
{
  if (C == NULL) {
    C = gsl_matrix_float_alloc( A->size2, B->size2 );
  }

  gsl_blas_sgemm( CblasTrans, CblasNoTrans, 1.0, A, B, 0.0, C );

  return C;
}


/*
**    C = A x B
*/
gsl_matrix_float *fmat_x_mat(gsl_matrix_float *A,gsl_matrix_float *B,gsl_matrix_float *C)
{
  if (C == NULL) {
    C = gsl_matrix_float_alloc( A->size1, B->size2 );
  }

  gsl_blas_sgemm( CblasNoTrans, CblasNoTrans, 1.0, A, B, 0.0, C );

  return C;
}


/*
**  y = Ax,  A double, x float
*/
gsl_vector_float *dmat_x_fvector(gsl_matrix *A, gsl_vector_float *x, gsl_vector_float *y)
{
  int i,j,nrows,ncols;
  float *ptr2,*ptr3,sum;
  double *ptr1;

  nrows = A->size1;
  ncols = A->size2;

  if (y == NULL) {
     y = gsl_vector_float_alloc (nrows);
  }

  if (x->size != ncols || y->size != nrows) {
    fprintf(stderr," fmat_x_vect: incongruent dimensions\n");
    exit(0);
  }

  ptr1 = A->data;
  ptr3 = y->data;

  for (i=0; i<nrows; i++) {
    sum = 0;
    ptr2 = x->data;
    for (j=0; j<ncols; j++) {
      sum += (*ptr1++) * (*ptr2++);
    }
    *ptr3++ = sum;
  }
  return y;
}


/*
**  z = x^T y
*/
float fskalarproduct(gsl_vector_float *x,gsl_vector_float *y)
{
  int i,n;
  float *ptr1,*ptr2,sum;

  n = x->size;
  if (y->size != n) {
    fprintf(stderr," fskalarproduct: incongruent vector sizes: %d %ld",n,y->size);
    exit(0);
  }

  ptr1 = x->data;
  ptr2 = y->data;
  sum = 0;
  for (i=0; i<n; i++) {
    sum += (*ptr1) * (*ptr2);
    ptr1++;
    ptr2++;
  }
  return sum;
}

/*
**    printout
*/
void fmatprint(FILE *fp,gsl_matrix_float *A,const char *format)
{
  int i,j;
  for (i=0; i<A->size1; i++) {
    for (j=0; j<A->size2; j++) {
      fprintf(fp,format,fmget(A,i,j));
    }
    fprintf(fp,"\n");
  }
  fprintf(fp,"\n");

}

/*
**    B = A^T
*/
gsl_matrix_float *ftranspose(gsl_matrix_float *A,gsl_matrix_float *B)
{
  int i,j,n,m;

  n = A->size1;
  m = A->size2;
  
  if (B == NULL) {
    B = gsl_matrix_float_alloc(m,n);
  }
  else if (B->size1 != m || B->size2 != n) {
    gsl_matrix_float_free(B);
    B = gsl_matrix_float_alloc(m,n);
  }

  for (i=0; i<n; i++) {
    for (j=0; j<m; j++) {
      fmset(B,j,i,fmget(A,i,j));
    }
  }
  return B;
}


/*
**    B = A^-1 = A^+
*/
gsl_matrix_float *fmat_PseudoInv(gsl_matrix_float *A,gsl_matrix_float *B)
{
  int i,j,k,l,n,m;
  double u,x;
  double *dbl_pp;
  float *flt_pp;

  m = A->size1;
  n = A->size2;

  if (B == NULL) {
    B = gsl_matrix_float_alloc (n,m);
  }
  else if (B->size1 != n || B->size2 != m) {
    gsl_matrix_float_free(B);
    B = gsl_matrix_float_alloc (n, m);
  }

  gsl_matrix *U = gsl_matrix_alloc (m, n);
  gsl_matrix *V = gsl_matrix_alloc (n, n);
  gsl_matrix *X = gsl_matrix_alloc (n, n);
  gsl_vector *w = gsl_vector_alloc (n);
  gsl_vector *work = gsl_vector_alloc (n);

  /* singular value decomposition */
  flt_pp = A->data;
  dbl_pp = U->data;
  for (i=0; i<A->size1 * A->size2; i++) 
    *dbl_pp++ = *flt_pp++;

  /* gsl_linalg_SV_decomp_jacobi(U,V,w); */
  gsl_linalg_SV_decomp_mod (U,X,V,w,work);


  /* check if singular */
  k=0;
  for (j=0; j<n; j++) {
    if (fabs(w->data[j]) > 1.0e-6) k++;
  }
  if (k < 2) VError(" Design matrix is singular");


  /* exception from Gaby */
  int j0 = 0;
  double xmin, xmax, tiny=10.0e-6;
  xmax = gsl_vector_get(w,0);
  xmin = tiny;
  for (j=n-1; j >= 0; j--) {
    u = gsl_vector_get(w,j);
    if (u > 0 && u/xmax > tiny) {
      j0 = j;
      goto skip;
    }
  }
 skip: ;
  if (j0 < n-1) {
    fprintf(stderr," Warning: Matrix almost singular\n");
    xmin = gsl_vector_get(w,j0) - tiny;
    if (xmin < 0) xmin = 0;
  }


  /* Fill the result matrix */
  gsl_matrix_float_set_zero(B);
  for (k=0; k<n; k++) {
    for (l=0; l<m; l++) {
      for (j=0; j<n; j++) {
	u = gsl_vector_get(w,j);
        if (fabs(u) > xmin) {
          x = gsl_matrix_float_get(B,k,l);
	  x += gsl_matrix_get(V,k,j)*gsl_matrix_get(U,l,j)/u;
          gsl_matrix_float_set(B,k,l,x);
        }
      }
    }
  }
  gsl_matrix_free(U);
  gsl_matrix_free(V);
  gsl_matrix_free(X);
  gsl_vector_free(w);
  gsl_vector_free(work);
  return B;
}


/*
** returns the trace of a matrix M 
*/
float trace(gsl_matrix_float *M)
{
  gsl_vector_float_view diag = gsl_matrix_float_diagonal(M);
  float sum = 0;
  int i = 0;
  for(i = 0; i<diag.vector.size; i++) {
    sum += gsl_vector_float_get(&diag.vector,i);
  }
  return sum;
}

/*
**  rank(A)
*/
int rank(gsl_matrix_float* mat) 
{
  int m = mat->size1;
  int n = mat->size2;
  int i;
    
  /* Create matrix buffer from source matrix. This prevents the overwriting of
   * matrix mat by gsl_linalg_SV_decomp().*/
  gsl_matrix* U = gsl_matrix_alloc(m,n);
  double* dpU =  U->data;
  float* fpMat = mat->data;
  for(i=0;i<mat->size1* mat->size2;i++) {
    *dpU++ = *fpMat++;
  }
    
  /* create additional matrix buffers */
  gsl_matrix *V;
  gsl_matrix *X;
  gsl_vector *S;
  gsl_vector *work;
  S = gsl_vector_alloc(n);
  V = gsl_matrix_alloc(n,n);
  X = gsl_matrix_alloc (n, n);
  work = gsl_vector_alloc (n);
   
  /* SVD */
  /* gsl_linalg_SV_decomp_jacobi(U,V,S); */
  gsl_linalg_SV_decomp_mod (U,X,V,S,work);

  int rank = 0;
  /* counting the nonzero singular values with a tolerance of EPSILON*/
  for(i=0;i<S->size;i++){
    if(S->data[i] > EPSILON) {
      rank++;
    }
  }
    
  gsl_matrix_free(U);
  gsl_matrix_free(V);
  gsl_matrix_free(X);
  gsl_vector_free(S);
  gsl_vector_free(work);
    
  return rank;
}


/*
**  fsum
*/
gsl_vector_float *fsum(gsl_matrix_float *matrix, int dim,gsl_vector_float *v) 
{
        
  int row,col;
  float sum=0;
    
  /* build sum over all columns */
  if (dim == 1) {
        
    if (v == NULL) {
      v = gsl_vector_float_alloc(matrix->size2);
    }
    if(v->size != matrix->size2) {
      fprintf(stderr, "Warning in fsum: vector size doesn't match related matrix dimension. Resizing ..");
      gsl_vector_float_free(v);
      v = gsl_vector_float_alloc(matrix->size2);
    }        
        
    for (col=0;col<matrix->size2;col++){            
      sum=0;
      for(row=0;row<matrix->size1;row++) {
	sum+=matrix->data[col+row*matrix->size2];
      }
      v->data[col] = sum;
    }
  }

  /* build sum over all rows */
  else {

    if (v == NULL) {
      v = gsl_vector_float_alloc(matrix->size1);
    }
        
    if (v->size != matrix->size1) {
      fprintf(stderr, "Warning in fsum: vector size doesn't match related matrix dimension. Resizing ..");
      gsl_vector_float_free(v);
      v = gsl_vector_float_alloc(matrix->size1);
    }

    for (row = 0;row<matrix->size1;row++) {
      sum = 0;
      for(col=0;col<matrix->size2;col++) {
	sum += matrix->data[col+row*matrix->size2];
      }
      v->data[row] = sum;
    }        
  }
  return v;
}


/*
**  funique
*/
gsl_vector_float *funique(gsl_vector_float* V) 
{
    
  /* working copy */
  gsl_vector_float* v = gsl_vector_float_alloc(V->size);
  /* pointer to result */
  gsl_vector_float* res;                      
  /* holds value of last saved element */
  float max=0;
  /* counter */
  int i;
    
  gsl_vector_float_memcpy(v,V);
        
  /* sort elements from V. */
  gsl_sort_vector_float(v);

  /* count number of different elements */
  int nelements = 0;
    
  /* first run: count different values */
  float* p = v->data;
  for(i=0;i<v->size;i++) {
    /* the first element in the list is special */
    if(i==0) {
      max = *p;
      nelements++;            
    }
    else {
      if(*p > max) {
	max = *p;
	nelements++;
      }                
    }
    p++;
  }
    
  /* second run: allocate mem & copy every value exactly once. */
  res = gsl_vector_float_alloc(nelements);
  p = v->data;
  float *s = res->data;
  for(i=0;i<v->size;i++) {
    if(i==0) {
      max = *p;
      *s++ = *p;
    }
    else {
      if(*p > max) {
	max = *p;
	*s++ = *p;
      }
    }
    p++;
  }

  gsl_vector_float_free(v);
    
  return res;
}

/*
**  subcols
*/
gsl_matrix_float *fmat_subcols(gsl_matrix_float* mat,gsl_vector_float* cols) 
{
    
  /* some tests */    
    
  if((cols->size < 0) || (cols->size > mat->size2)) {
    fprintf(stderr,"column vector: invalid dimensions");
    exit(-1);
  }
    
  float min, max;
  gsl_vector_float_minmax(cols, &min, &max);
    
  if((min < 0) || (max > mat->size2)) { 
    fprintf(stderr,"column vector values exceed matrix dimensions!");
    exit(-1);
  }
    
  /* TODO reduce column vector to unique values with funique. We will omit 
   * this due to performance reasons. */
    
  /* the return value */
  gsl_matrix_float* ret = gsl_matrix_float_alloc(mat->size1, cols->size);
  /* the copy buffer */
  gsl_vector_float* buff = gsl_vector_float_alloc(mat->size1);        
  /* counter */
  int i;
    
  /* copy columns from matrix to return buffer */
  for (i = 0; i < cols->size; ++i) {
    gsl_matrix_float_get_col(buff, mat, (int)cols->data[i]);
    gsl_matrix_float_set_col(ret, i, buff);
  }
    
  gsl_vector_float_free(buff);
    
  return ret;
}

/*
**  toeplitz
*/
gsl_matrix_float *fmat_toeplitz(gsl_vector_float* v, gsl_matrix_float* A)
{
  int i,j;

  if(A == NULL)
    A = gsl_matrix_float_alloc(v->size, v->size);
  else {
    if((A->size1 != v->size)  || (A->size2 != v->size)) {
      fprintf(stderr, "Warning fmat_toeplitz: incongruent matrix dimensions. Trying to\
                    correct it.");
      gsl_matrix_float_free(A);
      A = gsl_matrix_float_alloc(v->size, v->size);
    }
  }

  for(i=0;i<A->size1;i++){
    for(j=0;j<A->size2;j++) {
      gsl_matrix_float_set(A,i,j,gsl_vector_float_get(v, fabs(i-j)));
    }
  }
   
  return A;
}

/*
**  A^-1
*/
gsl_matrix_float *fInv(gsl_matrix_float *M, gsl_matrix_float *result) 
{
  int s; /*permutation signum for LU-Decomposition*/
  int i;

  int m = M->size1;
  int n = M->size2;

  if(m != n) {
    fprintf(stderr, "dInv: not a square matrix\n");
    exit(0);
  }

  if(result == NULL) { 
    result = gsl_matrix_float_alloc(m,m);
  }
  if((result->size1 != m) || (result->size2 != n)) {
    fprintf(stderr,"dInv: incongruent matrix dimensions.\n");
    exit(0);
  }
   
  gsl_matrix *ludecomp = gsl_matrix_alloc(m,m);
  gsl_permutation *perm = gsl_permutation_alloc(m); 
  gsl_matrix *res = gsl_matrix_alloc(m,m);

  /* convert from float to double */
  double *d_ptr = ludecomp->data;
  float *f_ptr = M->data;
  for(i=0;i<M->size1*M->size2;i++) {
    *d_ptr = (double)*f_ptr;
    d_ptr++;f_ptr++;
  }

  /* LU decomposition */
  gsl_linalg_LU_decomp(ludecomp,perm,&s);

  /* inverting matrix */
  gsl_linalg_LU_invert(ludecomp,perm,res);
    
  /* convert from double to float */
  d_ptr = res->data;
  f_ptr = result->data;
  for(i=0;i<res->size1*res->size2;i++) {
    *f_ptr = (float)*d_ptr;
    d_ptr++;f_ptr++;
  }   


  gsl_matrix_free(ludecomp);
  gsl_matrix_free(res);
  gsl_permutation_free(perm);
  return result;
}
