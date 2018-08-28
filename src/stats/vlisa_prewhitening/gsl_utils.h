
/* float */
extern gsl_vector_float * fmat_x_vector(gsl_matrix_float *,gsl_vector_float *,gsl_vector_float *);
extern gsl_matrix_float * fmat_x_mat(gsl_matrix_float *,gsl_matrix_float *,gsl_matrix_float *);
extern gsl_matrix_float *fmat_x_matT(gsl_matrix_float *,gsl_matrix_float *,gsl_matrix_float *);
extern gsl_matrix_float *fmatT_x_mat(gsl_matrix_float *,gsl_matrix_float *,gsl_matrix_float *);
extern float fskalarproduct(gsl_vector_float *,gsl_vector_float *);
extern void fmatprint(FILE *,gsl_matrix_float *,const char *);
extern gsl_matrix_float * fmat_PseudoInv(gsl_matrix_float *,gsl_matrix_float *);
extern gsl_matrix_float *ftranspose(gsl_matrix_float *,gsl_matrix_float *);

/* new functions used by vwhiteglm */
/* returns the trace of a matrix M */
extern float trace(gsl_matrix_float *M);

/* calculates the inverse with the use of LU decomposition */
extern gsl_matrix_float* fInv(gsl_matrix_float *M, gsl_matrix_float *result); 

/* rank of a matrix */
extern int rank(gsl_matrix_float* mat);

/* Sums all elements of along the given dimension */
extern gsl_vector_float* fsum(gsl_matrix_float* mat, int dim, gsl_vector_float* v);

extern gsl_vector_float* funique(gsl_vector_float* V);

extern gsl_matrix_float* fmat_subcols(gsl_matrix_float* mat, gsl_vector_float* cols);

extern gsl_matrix_float* fmat_toeplitz(gsl_vector_float* v, gsl_matrix_float* A);



/* double */
extern gsl_vector * dmat_x_vector(gsl_matrix *,gsl_vector *,gsl_vector *);
extern gsl_matrix * dmat_x_mat(gsl_matrix *,gsl_matrix *,gsl_matrix *);
extern gsl_matrix * dmatT_x_mat(gsl_matrix *,gsl_matrix *,gsl_matrix *);
extern gsl_matrix * dmat_x_matT(gsl_matrix *,gsl_matrix *,gsl_matrix *);
extern double dskalarproduct(gsl_vector *,gsl_vector *);
extern gsl_matrix * dmat_PseudoInv(gsl_matrix *,gsl_matrix *);
extern gsl_matrix *dtranspose(gsl_matrix *,gsl_matrix *);



#define dvset gsl_vector_set 
#define dvget gsl_vector_get 
#define dmset gsl_matrix_set 
#define dmget gsl_matrix_get

#define fvset gsl_vector_float_set 
#define fvget gsl_vector_float_get 
#define fmset gsl_matrix_float_set 
#define fmget gsl_matrix_float_get

/* epsilon interval */
#define EPSILON 0.000001
