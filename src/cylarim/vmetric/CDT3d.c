/*! \file
Constrained distance transform.

Perform Chamfer distance transform in 3D while avoiding obstacles, i.e.
get distance from each "inside" point to closest "outside" point
while avoiding "obstacle" points. The distance is measured as the length
of the shortest path to the closest outside point that does not touch an
obstacle point.

\par Reference:
Verwer, B.J., Verbeek, P.W.,  Dekker, S.T. (1989).
An efficient uniform cost algorithm applied to distance transforms.
IEEE Trans. on Pattern Aanalysis and Machine Intelligence,
Vol.11, No.4, pp. 425--429.

\par Author:
Gabriele Lohmann, MPI-CBS
*/

/* From the Vista library: */
#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>

/*
Constrained distance transform.

Perform Chamfer distance transform in 3D while avoiding obstacles, i.e.
get distance from each "inside" point to closest "outside" point
while avoiding "obstacle" points. The distance is measured as the length
of the shortest path to the closest outside point that does not touch an
obstacle point.

Reference:
Verwer, B.J., Verbeek, P.W.,  Dekker, S.T. (1989).
An efficient uniform cost algorithm applied to distance transforms.
IEEE Trans. on Pattern Aanalysis and Machine Intelligence,
Vol.11, No.4, pp. 425--429.
*/



/* get spatial resolution */
void XGetResolution(VImage src,double reso[])
{
  VString str;
  float vx=0,vy=0,vz=0;
  if (VGetAttr (VImageAttrList (src), "voxel", NULL,VStringRepn, (VPointer) & str) == VAttrFound) {
    sscanf(str,"%f %f %f",&vx,&vy,&vz);
  }
  reso[0] = vx;
  reso[1] = vy;
  reso[2] = vz;
}

  

void VCDT3d (VImage src,VImage dest,int inside,int outside,int obstacle)
{
  int b,r,c,k;
  size_t i,iter;
  double d1,d2,d3;
  double u,x;

  /* Check the input image */
  if (VPixelRepn(src) != VUByteRepn) 
    VError("VCDT3d: input image must be of type VUByte");

  int nbands = VImageNBands(src);
  int nrows  = VImageNRows (src);
  int ncols  = VImageNColumns (src);
  size_t npixels = VImageNPixels(src);

  double inf = 99999;
  VFillImage(dest,VAllBands,inf);
  
  VFloat *dest_pp = VImageData(dest);
  VUByte *src_pp = VImageData(src);

  u=0;
  for (i=0; i<npixels; i++) {
    k = (int)(src_pp[i]);
    if (k == inside) u = inf;
    else if (k == outside) u = 0;
    else u = -999.0;
    dest_pp[i] = u;
  }


  double tiny = 0.000001;
  double maxdist = 0;

 
  /* optimal chamfer distances (Borgefors,1984, p. 334) */
  /*
  d1 = 1.0;
  d2 = 1.314;
  d3 = 1.628;
  */

  /* optimal distances (see Beckers, Smeulders, 1992). */
  /*
  d1 = 0.88750;
  d2 = 1.34224;
  d3 = 1.59772;
  */

  /* optimal distances (see Verwer, 1991). */
  /*
  d1 = 0.894;
  d2 = 1.3409;
  d3 = 1.5879;
  */

  /* naive */
  /*
  d1 = 1;
  d2 = sqrt(2.0);
  d3 = sqrt(3.0);
  */


  /* optimal distances (see Kiryati, 1993). */
  d1 = 0.9016;
  d2 = 1.289;
  d3 = 1.615;

  
  /* adjustment for non-isotropic grids */
  double reso[3];
  XGetResolution(src,reso);
  fprintf(stderr," voxel resolution: %.4f  %.4f  %.4f\n",reso[0],reso[1],reso[2]);

  
  double v = reso[0];
  if (reso[1] < v) v = reso[1];
  if (reso[2] < v) v = reso[2];
  reso[0] /= v;
  reso[1] /= v;
  reso[2] /= v;

  double dist1[3],dist2[3],dist3=0;
  
  dist1[0] = d1*reso[0];
  dist1[1] = d1*reso[1];
  dist1[2] = d1*reso[2];

  double s2 = d2/sqrt(2.0);
  dist2[0] = s2*sqrt(reso[0]*reso[0] + reso[1]*reso[1]);
  dist2[1] = s2*sqrt(reso[0]*reso[0] + reso[2]*reso[2]);
  dist2[2] = s2*sqrt(reso[1]*reso[1] + reso[2]*reso[2]);

  double s3 = d3/sqrt(3.0);
  dist3 = s3*sqrt(reso[0]*reso[0] + reso[1]*reso[1] + reso[2]*reso[2]);
  
  /*
  fprintf(stderr,"   d : %f %f %f\n",d1,d2,d3);
  fprintf(stderr,"dist1: %f %f %f\n",dist1[0],dist1[1],dist1[2]);
  fprintf(stderr,"dist2: %f %f %f\n",dist2[0],dist2[1],dist2[2]);
  fprintf(stderr,"dist3: %f\n",dist3);
  exit(0);
  */
  
  /* forward scan */
  for (iter=0; iter<2; iter++) {
    for (b=1; b<nbands; b++) {
      for (r=1; r<nrows-1; r++) {
	for (c=1; c<ncols-1; c++) {

	  u = VPixel(dest,b,r,c,VFloat);
	  if (u < tiny) continue;

	  x = VPixel(dest,b-1,r-1,c-1,VFloat) + dist3;
	  if (x < u && x > tiny) u = x;

	  x = VPixel(dest,b-1,r-1,c,VFloat) + dist2[2];
	  if (x < u && x > tiny) u = x;

	  x = VPixel(dest,b-1,r-1,c+1,VFloat) + dist3;
	  if (x < u && x > tiny) u = x;

	  x = VPixel(dest,b-1,r,c-1,VFloat) + dist2[1];
	  if (x < u && x > tiny) u = x;

	  x = VPixel(dest,b-1,r,c,VFloat) + dist1[2];
	  if (x < u && x > tiny) u = x;

	  x = VPixel(dest,b-1,r,c+1,VFloat) + dist2[1];
	  if (x < u && x > tiny) u = x;

	  x = VPixel(dest,b-1,r+1,c-1,VFloat) + dist3;
	  if (x < u && x > tiny) u = x;

	  x = VPixel(dest,b-1,r+1,c,VFloat) + dist2[2];
	  if (x < u && x > tiny) u = x;

	  x = VPixel(dest,b-1,r+1,c+1,VFloat) + dist3;
	  if (x < u && x > tiny) u = x;


	  x = VPixel(dest,b,r-1,c-1,VFloat) + dist2[0];
	  if (x < u && x > tiny) u = x;

	  x = VPixel(dest,b,r-1,c,VFloat) + dist1[1];
	  if (x < u && x > tiny) u = x;

	  x = VPixel(dest,b,r-1,c+1,VFloat) + dist2[0];
	  if (x < u && x > tiny) u = x;

	  x = VPixel(dest,b,r,c-1,VFloat) + dist1[0];
	  if (x < u && x > tiny) u = x;

	  if (u > tiny)
	    VPixel(dest,b,r,c,VFloat) = u;
	}
      }
    }


    /* backward scan */
    for (b=nbands-2; b>=0; b--) {
      for (r=nrows-2; r>0; r--) {
	for (c=ncols-2; c>0; c--) {

	  u = VPixel(dest,b,r,c,VFloat);
	  if (u < tiny) continue;
	
	  x = VPixel(dest,b,r,c+1,VFloat) + dist1[0];
	  if (x < u && x > tiny) u = x;

	  x = VPixel(dest,b,r+1,c-1,VFloat) + dist2[0];
	  if (x < u && x > tiny) u = x;

	  x = VPixel(dest,b,r+1,c,VFloat) + dist1[1];
	  if (x < u && x > tiny) u = x;

	  x = VPixel(dest,b,r+1,c+1,VFloat) + dist2[0];
	  if (x < u && x > tiny) u = x;

	  x = VPixel(dest,b+1,r-1,c-1,VFloat) + dist3;
	  if (x < u && x > tiny) u = x;

	  x = VPixel(dest,b+1,r-1,c,VFloat) + dist2[2];
	  if (x < u && x > tiny) u = x;
	
	  x = VPixel(dest,b+1,r-1,c+1,VFloat) + dist3;
	  if (x < u && x > tiny) u = x;

	  x = VPixel(dest,b+1,r,c-1,VFloat) + dist2[1];
	  if (x < u && x > tiny) u = x;

	  x = VPixel(dest,b+1,r,c,VFloat) + dist1[2];
	  if (x < u && x > tiny) u = x;

	  x = VPixel(dest,b+1,r,c+1,VFloat) + dist2[1];
	  if (x < u && x > tiny) u = x;


	  x = VPixel(dest,b+1,r+1,c-1,VFloat) + dist3;
	  if (x < u && x > tiny) u = x;

	  x = VPixel(dest,b+1,r+1,c,VFloat) + dist2[2];
	  if (x < u && x > tiny) u = x;

	  x = VPixel(dest,b+1,r+1,c+1,VFloat) + dist3;
	  if (x < u && x > tiny) u = x;

	  if (u > tiny) {
	    VPixel(dest,b,r,c,VFloat) = u;
	    if (u > maxdist) maxdist = u;
	  }
	}
      }
    }
  }

  /* cleanup */
  src_pp = VImageData(src);
  dest_pp = VImageData(dest);
  for (i=0; i<VImageNPixels(dest); i++) {
    u = (double)(dest_pp[i]);
    if (u < 0) dest_pp[i] = 0;
    if (src_pp[i] != inside) dest_pp[i] = 0;
    if (u >= inf - tiny) dest_pp[i] = 0;
    else if (u > maxdist) dest_pp[i] = maxdist;
  }
}

