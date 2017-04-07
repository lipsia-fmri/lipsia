/*
** sorting.
**
** G.Lohmann, MPI-CBS
*/

#include <stdio.h>
#include <stdlib.h>

typedef struct {
  short b;
  short r;
  short c;
  float val;
} VPoint;


void VPoint_hpsort(unsigned long n, VPoint ra[]) 
{
  unsigned long i,ir,j,l;
  VPoint rra;

  if (n < 2) return;
  l=(n >> 1)+1;
  ir=n;
  for (;;) {
    if (l > 1) {
      rra.b=ra[--l].b;
      rra.r=ra[l].r;
      rra.c=ra[l].c;
      rra.val=ra[l].val;
    } 
    else {
      rra.b=ra[ir].b;
      rra.r=ra[ir].r;
      rra.c=ra[ir].c;
      rra.val=ra[ir].val;

      ra[ir].b=ra[1].b;
      ra[ir].r=ra[1].r;
      ra[ir].c=ra[1].c;
      ra[ir].val=ra[1].val;
 
      if (--ir == 1) {
	ra[1].b=ra[ir].b;
	ra[1].r=ra[ir].r;
	ra[1].c=ra[ir].c;
	ra[1].val=ra[ir].val;
	break;
      }
    }
    i=l;
    j=l+l;
    while (j <= ir) {
      if (j < ir && ra[j].val < ra[j+1].val) j++;
      if (rra.val < ra[j].val) {

	ra[i].b=ra[j].b;
	ra[i].r=ra[j].r;
	ra[i].c=ra[j].c;
	ra[i].val=ra[j].val;

	i=j;
	j <<= 1;
      } else j=ir+1;
    }
    ra[i].b=rra.b;
    ra[i].r=rra.r;
    ra[i].c=rra.c;
    ra[i].val=rra.val;
  }
}
