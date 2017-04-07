/*
** conversions: p->z, z->p, t->z
**
** G.Lohmann, April 2007
*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>

#define ABS(x) ((x) > 0 ? (x) : -(x))



/* 
** convert t to p values 
*/
double t2p(double t,double df)
{
  double a,b,x;
  extern double gsl_sf_beta_inc(double,double,double);

  x = df/(df+(t*t));
  if (x < 0 || x > 1) return 1;
  a = 0.5*df;
  b = 0.5;

  return 0.5 * gsl_sf_beta_inc(a,b,x);
}


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


/* 
** approximation: convert t to p values 
*/
float
t2z_approx(float t,float df)
{
  float z=0,u;

  u = df*log(1.0+t*t/df)*(1.0-0.5/df);
  if (u <= 0) return 0;
  z = sqrt(u);
  return z;
}


/* 
** convert p to t values 
*/
double 
p2t(double px,double df)
{
  double p,t0,t1,t,step=0.00001;

  t  = 0;
  t0 = 0;
  t1 = 20;
  while (ABS(t0-t1) > step) {
    t = (t0 + t1)*0.5;
    p = t2p(t,df);
    if (p < px)  t1 = t;
    else t0 = t;
  }
  p = t2p(t,df);

  return t;
}


/* 
** convert z to p values 
*/
double z2p(double z)
{
  if (z < 0) 
    return gsl_cdf_ugaussian_Q(-z);
  else 
    return gsl_cdf_ugaussian_Q(z);
}


/* 
** convert p to z values 
*/
double p2z(double p)
{
  return gsl_cdf_ugaussian_Qinv(p);
}


