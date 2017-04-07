/*! \file
  Compute 3D shape features and moments.

The following functions compute several shape features, including
central moments, center of gravity, and volume size.

\par Author:
Gabriele Lohmann, MPI-CBS
*/


/* From the Vista library: */
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <via/via.h>

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>



/* From the standard C libaray: */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SQR(x) ((x)*(x))
#define ABS(x) ((x) > 0 ? (x) : -(x))


static double power(double a, int k)
{
  register double res;
  register int i;

  res = 1.0;
  for (i=0; i<k; i++)
    res *= a;
  
  return res;
}



/*!
\fn  void VolumeCentroid(Volume v,double mean[3])
\brief Compute center of gravity of a volume.
\param v  input volume.
\param *mean output array.
The center of gravity is returned in mean[0],mean[1],mean[2].
Sufficient memory must have been allocated for <mean> prior to the call.
*/
void VolumeCentroid(Volume v, double mean[3])
{
  double npixels;
  long c,i,ce;
  VTrack t;

  mean[0] = 0;
  mean[1] = 0;
  mean[2] = 0;
  npixels = 0;

  for (i=0; i<v->nbuckets; i++) { 
    for (t = v->bucket[i].first; t != NULL; t = t->next)  {
      mean[0] += (double) (t->band * t->length);
      mean[1] += (double) (t->row * t->length);

      ce = t->col + t->length;
      for (c = t->col; c<ce; c++)
	mean[2] += (double) c;

      npixels += t->length;
    }
  }
  if (npixels > 0) {
    mean[0] /= npixels;
    mean[1] /= npixels;
    mean[2] /= npixels;
  }
}





/*!
\fn double VolumeMoment(Volume v,double mean[3],int m0,int m1,int m2)
\brief Compute central moments of a volume.
\param v  input volume
\param *mean input array containing 
the center of gravity as mean[0],mean[1],mean[2]  (slice,row,column).
If mean is NULL, the center of gravity is taken to be (0,0,0).
\param m0  first index of moment
\param m1  second index of moment 
\param m2  third index of moment
*/
double VolumeMoment(Volume v,double mean[3],long m0,long m1,long m2)
{
  double res=0;
  double a = 1.0,sum = 0;
  long i,ca,ce,ci;
  double b,r,c;
  VTrack t;
  double g0,g1,g2;

  g0 = mean[0]; 
  g1 = mean[1]; 
  g2 = mean[2];
  
  res = 0;
  for (i=0; i<v->nbuckets; i++) { 
    for (t = v->bucket[i].first; t != NULL; t = t->next)  {
      b  = (double) t->band - g0;
      r  = (double) t->row  - g1;
      a  = power(b,m0) * power(r,m1);
      ca = t->col;
      ce = ca + t->length;
      sum = 0;
      for (ci=ca; ci<ce; ci++) {
	c = (double)(ci - g2);
	sum += power(c,m2);
      }
      res += sum * a;
    }
  }
  return res;
}



/*!
\fn  void VBinCentroid(VImage src,double mean[3])
\brief Compute central moments of a binary raster image.
\param src  input image (bit repn)
\param *mean output array.
The center of gravity is returned in mean[0],mean[1],mean[2].
Sufficient memory must have been allocated for <mean> prior to the call.
*/
void VBinCentroid(VImage src,double mean[3])
{
  double npixels;
  long b,r,c;

  mean[0] = 0;
  mean[1] = 0;
  mean[2] = 0;
  npixels = 0;
  for (b=0; b<VImageNBands(src); b++) {
    for (r=0; r<VImageNRows(src); r++) {
      for (c=0; c<VImageNColumns(src); c++) {
	if (VPixel(src,b,r,c,VBit) > 0) {
	  mean[0] += b;
	  mean[1] += r;
	  mean[2] += c;
	  npixels++;
	}
      }
    }
  }
  if (npixels > 0) {
    mean[0] /= npixels;
    mean[1] /= npixels;
    mean[2] /= npixels;
  }
}




/*!
\fn  double VBinMoment(VImage src,double mean[3],int m0, int m1, int m2)
\brief Compute central moments of a binary raster image.
\param src  input image (bit repn)
\param *mean input array containing 
the center of gravity as mean[0],mean[1],mean[2]  (slice,row,column).
If mean is NULL, the center of gravity is taken to be (0,0,0);
\param m0  first index of moment
\param m1  second index of moment 
\param m2  third index of moment
*/
double VBinMoment(VImage src, double mean[3], long m0, long m1, long m2)
{
  double res;
  long b,r,c;
  double g0,g1,g2;

  g0 = mean[0]; 
  g1 = mean[1]; 
  g2 = mean[2];

  res = 0;
  for (b=0; b<VImageNBands(src); b++) {
    for (r=0; r<VImageNRows(src); r++) {
      for (c=0; c<VImageNColumns(src); c++) {
	if (VPixel(src,b,r,c,VBit) > 0) 
	  res += power((double)(b - g0),m0) 
	    * power((double)(r - g1),m1) 
	    * power((double)(c - g2),m2);
      }
    }
  }
  return res;
}


/*!
  \fn long VolumeSize(Volume v)
  \brief compute volume size, output the number of voxels of a single volume
  \param v   a single volume
*/
long VolumeSize(Volume v) 
{
  long i,isize;
  VTrack t;

  isize = 0;
  for (i=0; i<VolumeNBuckets(v); i++) {
    for (t = VFirstTrack(v,i); VTrackExists(t); t = VNextTrack(t))
      isize += VTrackLength(t);
  }
  return isize;
}



/*!
  \fn long VBinSize(VImage src)
  \brief count number of foreground voxels.
  \param src  input image (bit repn)
*/
long VBinSize(VImage src)
{
  long n;
  long i;
  VBit *bin_pp;

  if (VPixelRepn(src) != VBitRepn) VError(" input image must be bit repn");

  n = 0;
  bin_pp = (VBit *) VPixelPtr(src,0,0,0);
  for (i=0; i<VImageNPixels(src); i++) {
    if (*bin_pp > 0) n++;
    bin_pp++;
  }
  return n;
}



/*!
  \fn float VolumeDir(Volume vol,float *e,float x[3])
  \brief computer principal direction of a volume from its interia matrix
  \param vol   input volume
  \param e     output largest eigenvalue
  \param x     output first eigenvector
*/
float VolumeDir(Volume vol,float *e,float x[3])
{
  gsl_matrix *a=NULL;
  static gsl_matrix *evec=NULL;
  static gsl_vector *eval=NULL;
  static gsl_eigen_symmv_workspace *workspace=NULL;
  double m020,m002,m200,m110,m101,m011;
  double norm,angle;
  double mean[3];
  float tiny=1.0e-5;

  mean[0] = mean[1] = mean[2] = 0;
  VolumeCentroid(vol,mean);

  m020 = VolumeMoment(vol,mean,0,2,0);
  m002 = VolumeMoment(vol,mean,0,0,2);
  m200 = VolumeMoment(vol,mean,2,0,0);

  m110 = VolumeMoment(vol,mean,1,1,0);
  m101 = VolumeMoment(vol,mean,1,0,1);
  m011 = VolumeMoment(vol,mean,0,1,1);


  /* inertia matrix */
  if (a == NULL) {
    a = gsl_matrix_calloc(3,3);
    workspace = gsl_eigen_symmv_alloc(3);
  }
  gsl_matrix_set(a,0,0,m020 + m002);
  gsl_matrix_set(a,1,1,m200 + m002);
  gsl_matrix_set(a,2,2,m200 + m020);

  gsl_matrix_set(a,0,1,-m110);
  gsl_matrix_set(a,1,0,-m110);

  gsl_matrix_set(a,0,2,-m101);
  gsl_matrix_set(a,2,0,-m101);

  gsl_matrix_set(a,1,2,-m011);
  gsl_matrix_set(a,2,1,-m011);

  gsl_eigen_symmv(a,eval,evec,workspace);
  gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_VAL_DESC);

  x[0] = gsl_matrix_get(evec,0,0);
  x[1] = gsl_matrix_get(evec,1,0);
  x[2] = gsl_matrix_get(evec,2,0);

  *e = gsl_vector_get(eval,0);

  angle = 0;
  norm = sqrt((double)(SQR(x[0]) + SQR(x[1]) + SQR(x[2])));
  if (norm > tiny) {
    angle = ABS(x[1])/ norm;
  }

  return angle;
}




/*!
  \fn void VolumeEigen(Volume vol,gsl_vector *eval,gsl_matrix *evec)
  \brief computer principal directions of a volume from its interia matrix
  \param vol    input volume
  \param eval   output eigenvalues
  \param evec   output matrix of eigenvectors (columns)
*/
void VolumeEigen(Volume vol,gsl_vector *eval,gsl_matrix *evec)
{  
  static gsl_matrix *a=NULL;
  static gsl_eigen_symmv_workspace *workspace=NULL;
  double m020,m002,m200,m110,m101,m011;
  double mean[3];

  mean[0] = mean[1] = mean[2] = 0;
  VolumeCentroid(vol,mean);

  m020 = VolumeMoment(vol,mean,0,2,0);
  m002 = VolumeMoment(vol,mean,0,0,2);
  m200 = VolumeMoment(vol,mean,2,0,0);

  m110 = VolumeMoment(vol,mean,1,1,0);
  m101 = VolumeMoment(vol,mean,1,0,1);
  m011 = VolumeMoment(vol,mean,0,1,1);


  /* inertia matrix */
  if (a == NULL) {
    a = gsl_matrix_calloc(3,3);
    workspace = gsl_eigen_symmv_alloc(3);
  }

  gsl_matrix_set(a,0,0,m020 + m002);
  gsl_matrix_set(a,1,1,m200 + m002);
  gsl_matrix_set(a,2,2,m200 + m020);

  gsl_matrix_set(a,0,1,-m110);
  gsl_matrix_set(a,1,0,-m110);

  gsl_matrix_set(a,0,2,-m101);
  gsl_matrix_set(a,2,0,-m101);

  gsl_matrix_set(a,1,2,-m011);
  gsl_matrix_set(a,2,1,-m011);

  gsl_eigen_symmv(a,eval,evec,workspace);
  gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_VAL_DESC);

  return;
}


