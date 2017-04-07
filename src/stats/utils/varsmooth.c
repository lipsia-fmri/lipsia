/*
** smooth variance computations
**
** M. Kuhlmann, MPI-KYB, Sept 2015
*/

#include <viaio/Vlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_combination.h>

VImage variance_smoothing(VImage *src[], float width)
{
  size_t nslices, nrows, ncols;
  
  /* image dims */
  nslices = VImageNBands(src[0]);
  nrows  = VImageNRows(src[0]);
  ncols  = VImageNColumns(src[0]);
  
  VImage var = VCreateImage(nslices,nrows,ncols,VFloatRepn);
  VFillImage(var,VAllBands,0);
  
  /*double *data,int n,double *a,double *v)*/

  int j;
  double ave,var,nx,s,u;

  nx = (double)n;
  ave = 0;
  for (j=0; j<n; j++) ave += data[j];
  ave /= nx;

  var = u = 0;
  for (j=0; j<n; j++) {
    s = data[j]-ave;
    u   += s;
    var += s*s;
  }
  var=(var-u*u/nx)/(nx-1);
  
  /* smooth variance */
  
  
  
  return var;
}