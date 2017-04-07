/*
** concordance correlation coefficient
**  Ref: Lin 1989, Biometrics, 45(1):255-268
**
** G.Lohmann, Jul. 2009
*/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define SQR(x) ((x) * (x))
#define ABS(x) ((x) > 0 ? (x) : -(x))

/*
** concordance corr coeff
*/
double CCC(double *data1,double *data2,int n)
{
  int i;
  double u,sum1,sum2,mean1,mean2,var1,var2,sxy,nx;
  double tiny=1.0e-6;
  
  nx = (double)n;
  sum1 = sum2 = 0;
  for (i=0; i<n; i++) {
    sum1 += data1[i];
    sum2 += data2[i];
  }
  mean1 = sum1/nx;
  mean2 = sum2/nx;

  sum1 = sum2 = 0;
  for (i=0; i<n; i++) {
    sum1 += SQR(data1[i] - mean1);
    sum2 += SQR(data2[i] - mean2);
  }
  var1 = sum1/nx;
  var2 = sum2/nx;

  sum1 = 0;
  for (i=0; i<n; i++) {
    sum1 += (data1[i] - mean1)*(data2[i] - mean2);
  }
  sxy = sum1/nx;
  
  u = var1 + var2 + SQR(mean1 - mean2);
  if (u < tiny) return 0;
  return (2.0*sxy) / u;
}

