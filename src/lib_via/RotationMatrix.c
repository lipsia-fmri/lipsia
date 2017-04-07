/*! \file
  3D rotation matrices.


\par Author:
Gabriele Lohmann, MPI-CBS
*/



#include <math.h>
#include <float.h>


/*
** compute rotation matrix from quaternions
*/

/*!
\fn void VQuaternionsRotation(float *q,float **rot)
\brief Compute a 3x3 rotation matrix from quaternions.
\param *q     1x4 vector of quaternions
\param **rot  the resulting 3x3 rotation matrix
*/
void
VQuaternionsRotation(float *q,float **rot)
{
  rot[0][0] = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
  rot[0][1] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
  rot[0][2] = 2.0 * (q[1] * q[3] + q[0] * q[2]);

  rot[1][0] = 2.0 * (q[1] * q[2] + q[0] * q[3]);
  rot[1][1] = q[0] * q[0] + q[2] * q[2] - q[1] * q[1] - q[3] * q[3];
  rot[1][2] = 2.0 * (q[2] * q[3] - q[0] * q[1]);

  rot[2][0] = 2.0 * (q[1] * q[3] - q[0] * q[2]);
  rot[2][1] = 2.0 * (q[2] * q[3] + q[0] * q[1]);
  rot[2][2] = q[0] * q[0] + q[3] * q[3] - q[1] * q[1] - q[2] * q[2];
}




/*
** compute rotation matrix from roll,pitch and yaw parameters
*/
/*!
\fn void VRotationMatrix(double roll,double pitch,double yaw,double rot[3][3])
\brief Compute a 3x3 rotation matrix from roll,pitch and yaw.
\param roll  roll parameter
\param pitch pitch parameter
\param yaw   yaw parameter
\param rot   the resulting 3x3 rotation matrix
*/
void
VRotationMatrix(double roll,double pitch,double yaw,double rot[3][3])
{
  double cr,cy,cp,sr,sy,sp;

  cr = cos(roll);
  cp = cos(pitch);
  cy = cos(yaw);

  sr = sin(roll);
  sp = sin(pitch);
  sy = sin(yaw);

  rot[0][0] = cr * cy + sr * sp * sy;
  rot[0][1] = sr * cp;
  rot[0][2] = sr * sp * cy - sy * cr;
  rot[1][0] = cr * sp * sy - sr * cy;
  rot[1][1] = cr * cp;
  rot[1][2] = sr * sy + cy * cr * sp;
  rot[2][0] = cp * sy;
  rot[2][1] = - sp;
  rot[2][2] = cp * cy;
}
