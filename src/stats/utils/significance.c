/* 
** significances and alpha-quantiles
**
** M. Kuhlmann, MPI-KYB, Sept 2015
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

int compare_function(const void *a,const void *b) {
float *x = (float *) a;
float *y = (float *) b;
if (*x < *y) return -1;
else if (*x > *y) return 1; return 0;
}

/* returns P (X >= T) given empirical distribution */
float significance(float * emp_dist, int count, float T)
{
  int i;
  double tiny = 1e-8;
  int c_smaller = 0, c_leq = 0;
  for(i = 0; i < count; i++) {
    if(emp_dist[i] < T) {
      c_smaller++;
    } else {
      c_leq++;
    }
  }
  
  float p = (float)(c_leq)/count;
  
  if(p < tiny)
    return 1.0/count;
  else
    return p;
}

/* returns P (X >= T) given sorted empirical distribution */
float significance_sorted(float * emp_dist, int count, float T)
{
  if(T > emp_dist[count-1])
    return 1.0/count;
  if(T <= emp_dist[0])
    return 1;
  
  int lower = 0, upper = count;
  int center = floor(upper+lower)/2;
  while(upper-lower > 1) {
    if(T < emp_dist[center])
      upper = center;
    else
      lower = center;
    center = floor(upper+lower)/2;
  }
  float p = (float)(count-(lower+1))/count;
  
  return p;
}

/* returns t from emp_dist, s.t. P(X >= t) <= alpha */
float get_alpha(float * emp_dist, int count, float alpha)
{
  int get_n = floor(count*alpha);
  printf("%d\n", get_n);
  float threshold = 0;
  
  /* duplicate emp_dist */
  float * temp_dist = calloc(count,sizeof(float));
  memcpy(temp_dist,emp_dist,count*sizeof(float));

  /* sort temp_dist ascending
   * still missing */
  qsort(temp_dist,count,sizeof(float),compare_function);

  if(count-get_n < count)
    threshold = temp_dist[count - get_n];
  else
    threshold = 0./0.;

  free(temp_dist);

  return threshold;
}
