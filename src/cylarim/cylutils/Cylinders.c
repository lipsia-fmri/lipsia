/*
**  Cylarim:  Cylinder-based Algorithm for Laminar fMRI Data Analysis
**  Create cylinder data structs
**
** G.Lohmann, MPI-KYB, Tuebingen, Oct 2024
*/


/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>

/* From the standard C library: */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../cylutils/cyl.h"

#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/

/*
void XWriteOutput(VImage image,VAttrList geolist,char *filename)
{
  VAttrList out_list = VCreateAttrList();
  VSetGeoInfo(geolist,out_list);
  VAppendAttr(out_list,"image",NULL,VImageRepn,image);
  FILE *fp = fopen(filename,"w");
  if (! VWriteFile (fp, out_list)) exit (1);
  fclose(fp);
}
*/


int verbose=1;

/* path from P1 to P2 must be fully inside grey matter */
int CheckConnected(VImage rim,int b1,int r1,int c1,int b2,int r2,int c2)
{
  int b,r,c;
  double t,x,y,z,step=0.02;
  double u[3],v[3];

  if (VPixel(rim,b1,r1,c1,VUByte) == 0) return -1;
  if (VPixel(rim,b2,r2,c2,VUByte) == 0) return -1;
  
  int nslices = VImageNBands(rim);
  int nrows = VImageNRows(rim);
  int ncols = VImageNColumns(rim);

  u[0] = (double)c1;
  u[1] = (double)r1;
  u[2] = (double)b1;
  
  v[0] = (double)c2;
  v[1] = (double)r2;
  v[2] = (double)b2;

  for (t=0; t<=1; t+=step) {
    x = u[0] + (v[0]-u[0])*t;
    y = u[1] + (v[1]-u[1])*t;
    z = u[2] + (v[2]-u[2])*t;
    b = (int)(z+0.5);
    r = (int)(y+0.5);
    c = (int)(x+0.5);
    if (b < 0 || b >= nslices) continue;
    if (r < 0 || r >= nrows) continue;
    if (c < 0 || c >= ncols) continue;
    if (b==b1 && r==r1 && c==c1) continue;
    if (b==b2 && r==r2 && c==c2) continue;
    if (VPixel(rim,b,r,c,VUByte) != 3) return -1;  /* not inside grey matter */
  }
  return 1;
}

/* must be 6-adjacent to some voxel inside this cylinder */
int CheckAdj(VImage visited,VImage rim,int cid,int bx,int rx,int cx)
{
  int b,r,c;
  int nslices = VImageNBands(visited);
  int nrows = VImageNRows(visited);
  int ncols = VImageNColumns(visited);
  int i,j;
  for (b=bx-1; b<=bx+1; b++) {
    if (b < 0 || b >= nslices) continue;
    for (r=rx-1; r<=rx+1; r++) {
      if (r < 0 || r >= nrows) continue;
      for (c=cx-1; c<=cx+1; c++) {
	if (c < 0 || c >= ncols) continue;
	if (fabs(b-bx) + fabs(r-rx) + fabs(c-cx) > 1) continue;  /* 6-adj */
	i = VPixel(rim,b,r,c,VUByte);
	j = VPixel(visited,b,r,c,VInteger);
	if (i==3 && j==cid) return 1;
      }
    }
  }
  return -1;
}


void CrossProduct(double *a,double *b,double *c)
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

/* Euclidean distance from z[3] to the nearest point of the center line
** which passes through x[3] and y[3]
*/
double CylDist(double x[],double y[],double z[],double dnorm)
{
  int i;
  double d[3],az[3],dz[3];
  for (i=0; i<3; i++) {
    d[i]  = y[i] - x[i];
    az[i] = z[i] - x[i];
  }
  CrossProduct(az,d,dz);
  double s=0;
  for (i=0; i<3; i++) {
    s += dz[i]*dz[i];
  }
  s /= dnorm;
  return s;
}

/* check cylinder length */
double CylLength(double x[],double y[],double z[],double dnorm)
{
  int i;
  double s=0;
  for (i=0; i<3; i++) {
    s += (z[i]-x[i])*(y[i]-x[i]);
  }
  s /= dnorm;
  return s;
}	  


size_t Connect(VImage rim,VImage xmapimage,VImage marked,VImage visited,
	       gsl_matrix_int *xmap,gsl_vector *reso,
	       int b1,int r1,int c1,int b2,int r2,int c2,double radius,int *tab,int kmax,int cid)
{
  int i,k=0,m=0,mlen=0;
  double dnorm=0,t,step=0.03;
  double x[3],y[3],u[3];

  double xreso = gsl_vector_min(reso);
  double rad2 = radius*radius;
  int wn = (int)(radius/xreso + 0.5) + 1;

  int nslices = VImageNBands(rim);
  int nrows = VImageNRows(rim);
  int ncols = VImageNColumns(rim);

  
  /* x,y are the end points of the center line of the cylinder connecting the two rims */
  x[0] = reso->data[0]*(double)c1;
  x[1] = reso->data[1]*(double)r1;
  x[2] = reso->data[2]*(double)b1;
  
  y[0] = reso->data[0]*(double)c2;
  y[1] = reso->data[1]*(double)r2;
  y[2] = reso->data[2]*(double)b2;

  dnorm=0;
  for (i=0; i<3; i++) dnorm += (y[i]-x[i])*(y[i]-x[i]);
  if (dnorm < TINY) return 0;

  
  /* make sure centerline is included even if radius < 0 */
  for (k=0; k<kmax; k++) tab[k] = -1;
  k=0;
  for (t=0; t<=1.05; t+=step) {
    if (k >= kmax-1) continue;
    
    for (i=0; i<3; i++) u[i] = (1.0-t)*x[i] + t*y[i];
    int c = (int)(u[0]/reso->data[0]+0.5);
    int r = (int)(u[1]/reso->data[1]+0.5); 
    int b = (int)(u[2]/reso->data[2]+0.5);
    if (b < 0 || b >= nslices) continue;
    if (r < 0 || r >= nrows) continue;
    if (c < 0 || c >= ncols) continue;

    if (VPixel(rim,b,r,c,VUByte) == 0) continue;
    if (VPixel(visited,b,r,c,VInteger) == cid) continue;
    VPixel(visited,b,r,c,VInteger) = cid;
    VPixel(marked,b,r,c,VInteger) = cid;
    if (k >= kmax) { VWarning(" kmax"); continue; }
    int ii = VPixel(xmapimage,b,r,c,VInteger);
    tab[k++] = ii;
  }
  int kn=k;
  if (radius < 0 || wn < 1) return k;
  

  int *xtab = (int *)VCalloc(kmax,sizeof(int));
  double *dist = (double *)VCalloc(kmax,sizeof(double));


  /* walk along the center line from x to y,
  ** mark voxels surrounding the center line as potential cylinder points
  **/
  mlen = (int)(1.0/step+0.5);
  int kk=0;

  
#pragma omp parallel for
  for (m=0; m<=mlen; m++) {
    int i,ii,b,r,c,bb,rr,cc;
    double t,s,dx,u[3],z[3];

    t = step*(double)m;
    for (i=0; i<3; i++) u[i] = (1.0-t)*x[i] + t*y[i];

    c = (int)(u[0]/reso->data[0]+0.5);
    r = (int)(u[1]/reso->data[1]+0.5); 
    b = (int)(u[2]/reso->data[2]+0.5);
    if (b < 0 || b >= nslices) continue;
    if (r < 0 || r >= nrows) continue;
    if (c < 0 || c >= ncols) continue;

    for (bb=b-wn; bb<=b+wn; bb++) {
      if (bb < 0 || bb >= nslices) continue;
      z[2] = reso->data[2]*(double)bb;
      for (rr=r-wn; rr<=r+wn; rr++) {
	if (rr < 0 || rr >= nrows) continue;
	z[1] = reso->data[1]*(double)rr;
	for (cc=c-wn; cc<=c+wn; cc++) {
	  if (cc < 0 || cc >= ncols) continue;
	  z[0] = reso->data[0]*(double)cc;

	  if (VPixel(rim,bb,rr,cc,VUByte) == 0) continue;
	  if (VPixel(marked,bb,rr,cc,VInteger) == cid) continue;
	  ii = VPixel(xmapimage,bb,rr,cc,VInteger);
	  if (ii < 0) continue;

	  /* check if inside cylinder */
	  dx = CylDist(x,y,z,dnorm);
	  if (dx > rad2) continue;	  
	  s = CylLength(x,y,z,dnorm);
	  if (s < -step || s > 1.0+step) continue;
	  
#pragma omp critical
	  {
	    if (kk >= kmax) VError(" cylinder is too big: %lu,  max allowed: %d",kk,kmax);
	    VPixel(marked,bb,rr,cc,VInteger) = cid;
	    xtab[kk] = ii;
	    dist[kk] = sqrt(dx);	  
	    kk++;
	  }
	}
      }
    }
  }
  size_t nn = kk;
 
  /* sort by distance from center line */
  size_t *perm = (size_t *)VCalloc(nn,sizeof(size_t));
  gsl_sort_index(perm,dist,1,nn);

  
  /* add voxels to cylinder ordered by distance from center line */
  int iter=0;
  k=kn;
  for (iter=0; iter<2; iter++) {
    for (i=0; i<nn; i++) {
      int j = xtab[perm[i]];
      int b = gsl_matrix_int_get(xmap,j,0);
      int r = gsl_matrix_int_get(xmap,j,1);
      int c = gsl_matrix_int_get(xmap,j,2);
    
      if (VPixel(rim,b,r,c,VUByte) == 0) continue;
      if (VPixel(visited,b,r,c,VInteger) == cid) continue;  /* already booked */
      if (CheckAdj(visited,rim,cid,b,r,c) < 0) continue; /* must be 6-adjacent */
      int ii = VPixel(xmapimage,b,r,c,VInteger);
      if (ii < 0) continue;
    
      VPixel(visited,b,r,c,VInteger) = cid;
      tab[k] = ii;
      k++; 
    }
  }

  VFree(dist);
  VFree(perm);
  VFree(xtab);
  return (size_t)k;
}

/* fill values to cylinder struct */
void FillCylinder(int *tab,gsl_vector_int *avec)
{
  int i,j;
  for (i=0; i<avec->size; i++) {
    j = tab[i];
    if (j < 0) VError(" FillCylinder, j<0");
    avec->data[i] = j;
  }
}



/* get spatial resolution, can be non-isotropic */
void GetResolution(VImage src,gsl_vector *reso)
{
  VString str;
  float vx=0,vy=0,vz=0;
  reso->data[0] = reso->data[1] = reso->data[2] = 1;
  if (VGetAttr (VImageAttrList (src), "voxel", NULL,VStringRepn, (VPointer) & str) == VAttrFound) {
    sscanf(str,"%f %f %f",&vx,&vy,&vz);
    reso->data[0] = vx;
    reso->data[1] = vy;
    reso->data[2] = vz;
  }
  else VError(" voxel resolution unknown");
}


Cylinders *VCylinder(VImage rim,VImage metric,double radius)
{
  size_t i,j,k;
  int b,r,c;
  
  int nslices = VImageNBands(rim);
  int nrows = VImageNRows(rim);
  int ncols = VImageNColumns(rim);
  
  /* search radius in voxels, max cortical thickness assumed to be 6mm */
  gsl_vector *reso = gsl_vector_calloc(3);
  GetResolution(rim,reso);
  if (verbose) fprintf(stderr," resolution: %.4f x %.4f x %.4f mm\n",reso->data[0],reso->data[1],reso->data[2]);
  double xreso = gsl_vector_min(reso);
  int wn = (int)(6.0/xreso+0.5);  /* max cortical thickness assumed to be 6mm, wn is the search radius */
  if (wn < 1) VError(" implausible voxel resolution,  %f",xreso);
  
  /* rim maps */
  VUByte *pu = VImageData(rim);
  size_t n1=0,n2=0,nx=0;
  for (i=0; i<VImageNPixels(rim); i++) {
    if (pu[i] == 1) n1++;
    if (pu[i] == 2) n2++;
    if (pu[i] > 0) nx++;
  }
  size_t maxint = (size_t)gsl_pow_int((double)2.0,31)-1;
  if (n1+n2 > maxint) VError(" too many rim points (%lu), max is: %lu",n1+n2,maxint);
  if (nx > maxint) VError(" too many cortex points (%lu), max is: %lu",nx,maxint);
  
  /* fill voxel address maps */
  gsl_matrix_int *map1 = gsl_matrix_int_calloc(n1,3);
  gsl_matrix_int *map2 = gsl_matrix_int_calloc(n2,3);
  gsl_matrix_int *xmap = gsl_matrix_int_calloc(nx,3);
  
  i=j=k=0;
  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	if (VPixel(rim,b,r,c,VUByte) == 1) {
	  gsl_matrix_int_set(map1,i,0,b);
	  gsl_matrix_int_set(map1,i,1,r);
	  gsl_matrix_int_set(map1,i,2,c);
	  i++;
	}
	if (VPixel(rim,b,r,c,VUByte) == 2) {
	  gsl_matrix_int_set(map2,j,0,b);
	  gsl_matrix_int_set(map2,j,1,r);
	  gsl_matrix_int_set(map2,j,2,c);
	  j++;
	}
	if (VPixel(rim,b,r,c,VUByte) > 0) {
	  gsl_matrix_int_set(xmap,k,0,b);
	  gsl_matrix_int_set(xmap,k,1,r);
	  gsl_matrix_int_set(xmap,k,2,c);
	  k++;
	}
      }
    }
  }

  /* xmapimage */
  VImage xmapimage = VCreateImage(nslices,nrows,ncols,VIntegerRepn);
  VFillImage(xmapimage,VAllBands,-1);
  k=0;
  for (b=0; b<nslices; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
	if (VPixel(rim,b,r,c,VUByte) == 0) continue;
	VPixel(xmapimage,b,r,c,VInteger) = k;
	k++;
      }
    }
  }


  /* mapimage, rim only */
  VImage mapimage = VCreateImage(nslices,nrows,ncols,VIntegerRepn);
  VFillImage(mapimage,VAllBands,0);
  for (i=0; i<n1; i++) {
    int b1 = gsl_matrix_int_get(map1,i,0);
    int r1 = gsl_matrix_int_get(map1,i,1);
    int c1 = gsl_matrix_int_get(map1,i,2);
    VPixel(mapimage,b1,r1,c1,VInteger) = i;
  }
  
  for (i=0; i<n2; i++) {
    int b2 = gsl_matrix_int_get(map2,i,0);
    int r2 = gsl_matrix_int_get(map2,i,1);
    int c2 = gsl_matrix_int_get(map2,i,2);
    VPixel(mapimage,b2,r2,c2,VInteger) = i;
  }

  /*
  VSetAttr(VImageAttrList(xmapimage),"voxel",NULL,VStringRepn,"0.1 0.1 0.1");
  XWriteOutput(xmapimage,NULL,"xmap.v");
  exit(0);
  */

  
  /* bipartite graph, assignment problem */
  int *match1 = (int *)VCalloc(n1,sizeof(int));
  int *match2 = (int *)VCalloc(n2,sizeof(int));
  for (i=0; i<n1; i++) match1[i] = -1;
  for (i=0; i<n2; i++) match2[i] = -1;
  double maxthickness=10.0;  /* maximal cortical thickness (10mm) */
  size_t progress=0;
 
  /* for each point on the CSF-rim find closest point on the WM-rim */
#pragma omp parallel for shared(progress)
  for (i=0; i<n1; i++) {
    int b1,r1,c1,b2,r2,c2,ii,j0;
    double d,dmin,x[3],y[3],db=0,dr=0,dc=0;
    progress++;
    
    b1 = gsl_matrix_int_get(map1,i,0);
    r1 = gsl_matrix_int_get(map1,i,1);
    c1 = gsl_matrix_int_get(map1,i,2);

    x[0] = reso->data[0]*(double)c1;
    x[1] = reso->data[1]*(double)r1;
    x[2] = reso->data[2]*(double)b1;
    
    j0 = -1;
    dmin = maxthickness;
    for (b2=b1-wn; b2<=b1+wn; b2++) {
      if (b2 < 0 || b2 >= nslices) continue;
      y[2] = reso->data[2]*(double)b2;
      db = (x[2]-y[2]);
     
      for (r2=r1-wn; r2<=r1+wn; r2++) {
	if (r2 < 0 || r2 >= nrows) continue;
	y[1] = reso->data[1]*(double)r2;
	dr = (x[1]-y[1]);
	
	for (c2=c1-wn; c2<=c1+wn; c2++) {
	  if (c2 < 0 || c2 >= ncols) continue;
	  y[0] = reso->data[0]*(double)c2;
	  dc = (x[0]-y[0]);
	  
	  ii = VPixel(rim,b2,r2,c2,VUByte);
	  if (ii != 2) continue;

	  d = sqrt(db*db + dr*dr + dc*dc);
	  if (d < dmin) {	    
	    int jj = CheckConnected(rim,b1,r1,c1,b2,r2,c2);
	    if (jj > 0) {
	      dmin = d;
	      j0 = VPixel(mapimage,b2,r2,c2,VInteger);
	    }
	  }
	}
      }       
    }
    if (j0 < 0) continue;

    
#pragma omp critical
    {
      match1[i] = j0;
    }
  }


  /* for each point on the WM-rim find closest point on the CSF-rim */
#pragma omp parallel for shared(progress)
  for (i=0; i<n2; i++) {
    int b1,r1,c1,b2,r2,c2,ii,j0;
    double d,dmin,x[3],y[3],db=0,dr=0,dc=0;
    progress++;
    
    b1 = gsl_matrix_int_get(map2,i,0);
    r1 = gsl_matrix_int_get(map2,i,1);
    c1 = gsl_matrix_int_get(map2,i,2);

    x[0] = reso->data[0]*(double)c1;
    x[1] = reso->data[1]*(double)r1;
    x[2] = reso->data[2]*(double)b1;
  
    j0 = -1;
    dmin = maxthickness;
    for (b2=b1-wn; b2<=b1+wn; b2++) {
      if (b2 < 0 || b2 >= nslices) continue;
      y[2] = reso->data[2]*(double)b2;
      db = (x[2]-y[2]);
      
      for (r2=r1-wn; r2<=r1+wn; r2++) {
	if (r2 < 0 || r2 >= nrows) continue;
	y[1] = reso->data[1]*(double)r2;
	dr = (x[1]-y[1]);
	
	for (c2=c1-wn; c2<=c1+wn; c2++) {
	  if (c2 < 0 || c2 >= ncols) continue;
	  y[0] = reso->data[0]*(double)c2;
	  dc = (x[0]-y[0]);
	  
	  ii = VPixel(rim,b2,r2,c2,VUByte);
	  if (ii != 1) continue;
	  
	  d = sqrt(db*db + dr*dr + dc*dc);
	  if (d < dmin) {

	    int jj = CheckConnected(rim,b1,r1,c1,b2,r2,c2);
	    if (jj > 0) {
	      dmin = d;
	      j0 = VPixel(mapimage,b2,r2,c2,VInteger);
	    }
	  }
	}
      }       
    }
    if (j0 < 0) continue;

    
#pragma omp critical
    {
      match2[i] = j0;
    }
  }
  
  
  /*
  ** assignment between each rim point on one rim to the closest point on the opposite rim
  */
  k=0;
  for (i=0; i<n1; i++) { 
    if (match1[i] < 0) continue;
    k++;
  }
  int m=0;
  for (i=0; i<n2; i++) {
    m = match2[i];
    if (m < 0) continue;
    if (match1[m] == i) continue;  /* exclude duplications */
    k++;
  }
  size_t numcylinders = k;
  gsl_matrix_int *E = gsl_matrix_int_calloc(numcylinders,2);  /* list of coordinates */
 
  k=0;
  for (i=0; i<n1; i++) {
    if (match1[i] < 0) continue;    
    gsl_matrix_int_set(E,k,0,i);
    gsl_matrix_int_set(E,k,1,match1[i]);
    k++;
  }
  for (i=0; i<n2; i++) {
    m = match2[i];
    if (m < 0) continue;
    if (match1[m] == i) continue;
    gsl_matrix_int_set(E,k,0,m);
    gsl_matrix_int_set(E,k,1,i);
    k++;
  }
 
  /* create cylinders */
  VImage visited = VCreateImage(nslices,nrows,ncols,VIntegerRepn);
  VFillImage(visited,VAllBands,0);
  VImage marked = VCreateImage(nslices,nrows,ncols,VIntegerRepn);
  VFillImage(marked,VAllBands,0);
  
  gsl_vector_int **avec = (gsl_vector_int **) VCalloc(numcylinders,sizeof(gsl_vector_int **));
  int kmax = 100000;  /* max number of voxels per cylinder */
  int *tab = (int *)VCalloc((size_t)kmax,sizeof(int));
  
  size_t np = (size_t)((float)numcylinders/100.0);
  double sum=0;
  for (k=0; k<numcylinders; k++) {
    if (k%np==0 && verbose) fprintf(stderr," create cylinders: %7.3f\r",(float)(k)/(float)numcylinders);

    i = gsl_matrix_int_get(E,k,0);
    int b1 = gsl_matrix_int_get(map1,i,0);
    int r1 = gsl_matrix_int_get(map1,i,1);
    int c1 = gsl_matrix_int_get(map1,i,2);

    j = gsl_matrix_int_get(E,k,1);
    int b2 = gsl_matrix_int_get(map2,j,0);
    int r2 = gsl_matrix_int_get(map2,j,1);
    int c2 = gsl_matrix_int_get(map2,j,2);
   
    size_t m = Connect(rim,xmapimage,marked,visited,xmap,reso,
		       b1,r1,c1,b2,r2,c2,radius,tab,kmax,(int)(k+1));
     
    avec[k] = gsl_vector_int_calloc(m);
    FillCylinder(tab,avec[k]);
    sum += (double)m;
  }
  if (verbose) fprintf(stderr," create cylinders: %7.3f\n",1.0);
  fprintf(stderr," num cylinders: %lu             \n",numcylinders);
  fprintf(stderr," cylinder radius: %.3f mm       \n",radius);
  fprintf(stderr," mean cylinder size: %.3f voxels\n",sum/(double)numcylinders);

  
  /* alloc Cylinder struct */
  Cylinders *cyl = (Cylinders *)malloc(sizeof(Cylinders));
  cyl->numcylinders = numcylinders;
  cyl->radius = radius;
  cyl->xmap = xmap;
  cyl->addr = avec;


  /* free memory */
  VDestroyImage(visited);
  VDestroyImage(mapimage);
  VDestroyImage(xmapimage);
  VFree(tab);
  gsl_matrix_int_free(map1);
  gsl_matrix_int_free(map2);
  gsl_matrix_int_free(E);
  gsl_vector_free(reso);

  return cyl;
}
