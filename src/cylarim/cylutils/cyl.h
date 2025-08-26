typedef struct cstruct
{
  gsl_matrix_int *xmap;
  gsl_vector_int **addr;
  size_t numcylinders;
  double radius;
} Cylinders;


#define n3bins 3
#define p3bins 3
#define npeak 2
#define nshape 3
#define nR2 6

/* minimal number of data points per cylinder */
#define NPTS 20


/* Chebyshev polynomials */
#define T0(x) (1)
#define T1(x) (x)
#define T2(x) (2.0*(x)*(x)-1.0)
#define T3(x) (4.0*gsl_pow_3((x)) - 3.0*(x))
#define T4(x) (8.0*gsl_pow_4((x))-8.0*(x)*(x)+1.0)
#define T5(x) (16.0*gsl_pow_5((x)) - 20.0*gsl_pow_3((x)) + 5.0*(x))

/* Legendre polynomials */
#define P0(x) (1)
#define P1(x) (x)
#define P2(x) (0.5*(3.0*(x)*(x)-1.0))
#define P3(x) (0.5*(5.0*(x)*(x)*(x)-3.0*(x)))
#define P4(x) ((35.0*gsl_pow_4((x)) - 30.0*(x)*(x) + 3.0)/8.0)
#define P5(x) ((63.0*gsl_pow_5((x)) - 70.0*gsl_pow_3((x)) + 15.0*(x))/8.0)

/* standard polynomials */
#define S0(x) (1)
#define S1(x) (x)
#define S2(x) ((x)*(x))
#define S3(x) ((x)*(x)*(x))
#define S4(x) ((x)*(x)*(x)*(x))


