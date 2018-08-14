/*
** t-tests
**
** G.Lohmann, Jan 2017
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>



#define TINY 1.0e-10
#define ABS(x) ((x) > 0 ? (x) : -(x))
#define SQR(x) ((x)*(x))


/* 
** convert t to z values 
*/
double t2z(double t,double df)
{
  double p=0,z=0;
  double a,b,x;
  extern double gsl_sf_beta_inc(double,double,double);


  /* t to p */
  x = df/(df+(t*t));
  if (x <= 0 || x > 1) return 0;

  a = 0.5*df;
  b = 0.5;
  p = 0.5 * gsl_sf_beta_inc(a,b,x);


  /* p to z */
  z = gsl_cdf_ugaussian_Qinv(p);
  return z;
}


void avevar(double *data,int n,double *a,double *v)
{
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
  
  *v = var;
  *a = ave;
}


/* welch test, unequal variance */
double welchtest(double *data1,double *data2,int n1,int n2)
{
  double tiny=TINY;
  double ave1,ave2,var1,var2;
  double z=0,t=0,df=0;
  double nx1 = (double)n1;
  double nx2 = (double)n2;

  avevar(data1,n1,&ave1,&var1);
  if (var1 < tiny) return 0;
  avevar(data2,n2,&ave2,&var2);
  if (var2 < tiny) return 0;
  t = (ave1 - ave2)/sqrt(var1/nx1 + var2/nx2);
  df = SQR(var1/nx1+var2/nx2)/(SQR(var1/nx1)/(nx1-1)+SQR(var2/nx2)/(nx2-1));
  z  = t2z((double) t,(double) df);
  if (t < 0) z = -z;
  return z;
}

/* paired twosample test */
double paired_ttest(double *data1,double *data2,int n)
{
  int j;
  double ave,var,tiny=TINY;
  for (j=0; j<n; j++) {
    data1[j] -= data2[j];
  }
  avevar(data1,n,&ave,&var);
  if (var < tiny) return 0;
  double nx = (double)n;
  double df = nx-1.0;
  double t = sqrt(nx) * ave/sqrt(var);
  double z = t2z(t,df);
  if (t < 0) z = -z;
  return z;
}


/* paired twosample test */
double xtest2(double *data1,double *data2,int n)
{
  int j;
  double nx,ave1,ave2,var1,var2,sd,u,df,cov;
  double tiny=TINY;
  double t=0,z=0;

  avevar(data1,n,&ave1,&var1);
  avevar(data2,n,&ave2,&var2);
  if (var1 < tiny || var2 < tiny) return 0;

  nx = (double)n;
  df = nx-1;
  cov = 0;
  for (j=0; j<n; j++)
    cov += (data1[j]-ave1)*(data2[j]-ave2);
  cov /= df;

  t = 0;
  u = (var1+var2-2.0*cov);
  if (u < tiny) return 0;
  sd = sqrt(u/nx);
  if (sd < tiny) return 0;
  t = (double)(ave1-ave2)/sd;
  z = t2z(t,df);
  if (t < 0) z = -z;
  return z;
}


/* twosample t-test, pooled variance  */
double ttest2(double *data1,double *data2,int n1,int n2)
{
  int j;
  double nx1,nx2,ave1,ave2,sd,s1,s2,df;
  double tiny=TINY;
  double t=0,z=0;

  nx1 = (double)n1;
  nx2 = (double)n2;

  ave1 = ave2 = 0;
  for (j=0; j<n1; j++) ave1 += data1[j];
  for (j=0; j<n2; j++) ave2 += data2[j];
  if (ABS(ave1) < tiny) return 0;
  if (ABS(ave2) < tiny) return 0;
  ave1 /= nx1;
  ave2 /= nx2;

  s1 = 0;
  for (j=0; j<n1; j++) s1 += SQR(data1[j]-ave1);

  s2 = 0;
  for (j=0; j<n2; j++) s2 += SQR(data2[j]-ave2);

  df = nx1 + nx2 - 2.0;
  sd = sqrt((1.0/nx1 + 1.0/nx2) * (s1+s2)/df);
  if (sd < tiny) return 0;

  t = (double)(ave1-ave2)/sd;
  z = t2z(t,df);
  if (t < 0) z = -z;
  return z;
}


/* onesample t-test */
double ttest1(double *data,int n)
{
  double ave=0,var=0,tiny=TINY;
  avevar(data,n,&ave,&var);
  if (var < tiny) return 0;
  double nx = (double)n;
  double df = nx-1.0;
  double t  = sqrt(nx) * ave/sqrt(var);
  double z  = t2z(t,df);
  if (t < 0) z = -z;
  return z;
}
