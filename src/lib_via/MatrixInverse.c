/*! \file
  matrix inverse of 2x2 and 3x3 matrices


\par Author:
Gabriele Lohmann, MPI-CBS, Jan 2005
*/

#include <math.h>
#include <float.h>

#define ABS(x) ((x) > 0 ? (x) : -(x))

/*
** compute matrix inverse of a 3x3 matrix
*/
/*!
\fn float VMatrixInverse_3x3(float a[3][3],float ainv[3][3])
\brief matrix inverse of a 3x3 matrix, returns determinant of input matrix.
\param a     3x3 matrix
\param ainv  resulting 3x3 inverse
*/
float
VMatrixInverse_3x3(float a[3][3],float ainv[3][3])
{
  int i,j;
  float detA;

  /* get its inverse : */
  ainv[0][0] =  a[1][1]*a[2][2] - a[1][2]*a[2][1];
  ainv[1][0] = -a[1][0]*a[2][2] + a[1][2]*a[2][0];
  ainv[2][0] =  a[1][0]*a[2][1] - a[1][1]*a[2][0];

  ainv[0][1] = -a[0][1]*a[2][2] + a[0][2]*a[2][1];
  ainv[1][1] =  a[0][0]*a[2][2] - a[0][2]*a[2][0];
  ainv[2][1] = -a[0][0]*a[2][1] + a[0][1]*a[2][0];

  ainv[0][2] =  a[0][1]*a[1][2] - a[0][2]*a[1][1];
  ainv[1][2] = -a[0][0]*a[1][2] + a[0][2]*a[1][0];
  ainv[2][2] =  a[0][0]*a[1][1] - a[0][1]*a[1][0];

  /* determinant */
  detA = a[0][0]*ainv[0][0] + a[0][1]*ainv[1][0] + a[0][2]*ainv[2][0];
  if (ABS(detA) < 1.0e-6) return detA;  /* singular matrix */

  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      ainv[i][j] /= detA;
    }
  }
  return detA;
}


/*
** compute matrix inverse of a 2x2 matrix
*/
/*!
\fn float VMatrixInverse_2x2(float a[2][2],float ainv[2][2])
\brief matrix inverse of a 2x2 matrix, returns determinant of input matrix.
\param a     2x2 matrix
\param ainv  resulting 2x2 inverse
*/
float
VMatrixInverse_2x2(float a[2][2],float ainv[2][2])
{
  int i,j;
  float detA;

  for (i=0; i<2; i++)
    for (j=0; j<2; j++) ainv[i][j] = 0;

  detA = (a[0][0]*a[1][1] - a[0][1]*a[1][0]);
  if (ABS(detA) < 1.0e-6) return detA;  /* singular matrix */

  ainv[0][0] =  a[1][1]/detA;
  ainv[0][1] = -a[0][1]/detA;
  ainv[1][0] = -a[1][0]/detA;
  ainv[1][1] =  a[0][0]/detA;

  return detA;
}
