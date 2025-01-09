typedef struct cstruct
{
  gsl_matrix_int *xmap;
  gsl_vector_int **addr;
  size_t numcylinders;
  double radius;
} Cylinders;
