/*! \file
  3D Skeletonization

Perform a topological thinning.

\par Reference:
 A.Manzanera,T.Bernard,F.Preteux,B.Longuet: 
 "nD skeletonization: a unified mathematical framework".
 Journal of Electronic Engineering, Vol. 11-1, Jan. 2002, pp.25-37.


\par Author:
Gabriele Lohmann, MPI-CBS
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* From the Vista library: */
#include <viaio/Vlib.h>
#include <viaio/file.h>
#include <viaio/mu.h>
#include <via/via.h>


#define ABS(x) ((x) > 0 ? (x) : -(x))


extern VBoolean Beta1_3d(VImage,int,int,int);
extern VBoolean Beta0_3d(VImage,int,int,int);
extern VBoolean Border_3d(VImage,int,int,int);

extern VBoolean Alpha_3d(VImage,VImage,VImage,int,int,int);
extern void GenerateMaps(void);


VImage alpha0_pos=NULL,alpha1_pos=NULL,alpha2_pos=NULL;
VImage alpha0_neg=NULL,alpha1_neg=NULL,alpha2_neg=NULL;



/*!
\fn VImage VSkel3d (VImage src,VImage dest)
\param src   input image (bit repn)
\param dest  output image (bit repn)
*/
VImage
VSkel3d(VImage src,VImage dest)
{
  int b,r,c,nbands,nrows,ncols,npixels;
  int ndel=0;
  VBit *tmp_pp,*dst_pp;
  VImage tmp=NULL,dt=NULL;
  int i,j,n,n0,n1;

  float depth,maxdepth,step=0.5;
  VPoint *array;
  VBit *src_pp;


  nrows  = VImageNRows (src);
  ncols  = VImageNColumns (src);
  nbands = VImageNBands (src);
  npixels = nbands*nrows*ncols;

  dest = VCopyImage(src,dest,VAllBands);
  VFillImage(dest,VAllBands,0);
  tmp = VCreateImage(nbands,nrows,ncols,VBitRepn);
  VFillImage(tmp,VAllBands,0);

  for (b=2; b<nbands-2; b++) {
    for (r=2; r<nrows-2; r++) {
      for (c=2; c<ncols-2; c++) {
	VPixel(dest,b,r,c,VBit) = VPixel(src,b,r,c,VBit);
      }
    }
  }

  /*
  ** generate masks
  */
  GenerateMaps();


  /*
  ** sort foreground voxels by depth
  */
  src_pp  = (VBit *) VPixelPtr(src, 0, 0, 0);
  n = 0;
  for (i=0; i<npixels; i++) {
    if (*src_pp > 0) n++;
    *src_pp = (*src_pp > 0) ? 0 : 1;
    src_pp++;
  }

  dt = VChamferDist3d(src,NULL,VFloatRepn);

  src_pp  = (VBit *) VPixelPtr(src, 0, 0, 0);
  for (i=0; i<npixels; i++) {
    *src_pp = (*src_pp > 0) ? 0 : 1;
    src_pp++;
  }

  array = (VPoint *) VMalloc(sizeof(VPoint) * (n+1));

  array[0].b = 0;
  array[0].r = 0;
  array[0].c = 0;
  array[0].val = 0;

  n = 1;
  for (b=0; b<nbands; b++) {
    for (r=0; r<nrows; r++) {
      for (c=0; c<ncols; c++) {
        if (VPixel(dt, b, r, c,VFloat) > 0.1) {
          array[n].b = b;
          array[n].r = r;
          array[n].c = c;
          array[n].val = VPixel(dt,b,r,c,VFloat);
          n++;
        }
      }
    }
  }
  VDestroyImage(dt);
  VPoint_hpsort((n-1),array);


  /*
  ** now start deleting voxels
  */  
  maxdepth = array[n-1].val;
  fprintf(stderr," maxdepth= %f\n",maxdepth);

  for (depth=0; depth <= maxdepth; depth += step) {

    fprintf(stderr,"   depth= %7.2f\r",depth); fflush(stderr);

    n0 = 1;
    for (j=n-1; j>=1; j--) {
      if (array[j].val < depth-2.0*step) {
        n0 = j;
        break;
      }
    }

    n1 = n;
    for (j=1; j<=n; j++) {
      if (array[j].val > depth) {
        n1 = j-1;
        break;
      }
    }

    /*
    ** delete all deletable voxels within this depth level
    */
    ndel = 1;
    while (ndel > 0) {
      ndel = 0;

      VFillImage(tmp,VAllBands,0);

      for (j=n0; j<=n1; j++) {

	if (array[j].val > depth) continue;

	b = array[j].b;
	r = array[j].r;
	c = array[j].c;
	if (b < 2 || b >= nbands-2) continue;
	if (r < 2 || r >= nrows-2) continue;
	if (c < 2 || c >= ncols-2) continue;


	if (VPixel(dest,b,r,c,VBit) == 0) continue;
	if (Border_3d(dest,b,r,c) == FALSE) continue;

	if (Beta1_3d(dest,b,r,c) == FALSE) continue;
	if (Beta0_3d(dest,b,r,c) == FALSE) continue;

	/* alpha conditions */
	if (Alpha_3d(dest,alpha0_pos,alpha0_neg,b,r,c) == FALSE &&
	    Alpha_3d(dest,alpha1_pos,alpha1_neg,b,r,c) == FALSE &&
	    Alpha_3d(dest,alpha2_pos,alpha2_neg,b,r,c) == FALSE) continue;

	/* mark for deletion */
	VPixel(tmp,b,r,c,VBit) = 1; 
	ndel++;
      }
    
      /* perform parallel deletion */
      tmp_pp = (VBit *) VPixelPtr(tmp,0,0,0);
      dst_pp = (VBit *) VPixelPtr(dest,0,0,0);

      for (i=0; i<npixels; i++) {
	if (*tmp_pp == 1) *dst_pp = 0;
	tmp_pp++;
	dst_pp++;
      }
    }
  }

  VDestroyImage(tmp);
  return dest;
}


VBoolean
Border_3d(VImage image,int band,int row,int col)
{
  int b,r,c;

  for (b=band-1; b<=band+1; b++) {
    for (r=row-1; r<=row+1; r++) {
      for (c=col-1; c<=col+1; c++) {
	if (VPixel(image,b,r,c,VBit) == 0) return TRUE;
      }
    }
  }
  return FALSE;
}



VImage
GenAlpha(int *bmap,int *rmap,int *cmap,int nn)
{
  VImage alpha=NULL;
  int b=0,r=0,c=0;
  int i,j,k,l,n,ii,jj;
  int nangles = 64;
  int x,y,z;
  double rot[3][3];
  int irot[3][3];
  double pi=3.14159265,roll,pitch,yaw;
  double xroll,xpitch,xyaw;

  alpha = VCreateImage(nangles,3,nn,VSByteRepn);
  VFillImage(alpha,VAllBands,0);

  n = 0;
  roll = 0;
  for (i=0; i<4; i++) {

    pitch = 0;
    for (j=0; j<4; j++) {

      yaw = 0;
      for (k=0; k<4; k++) {

	xroll  = roll  * pi / 180.0;
	xpitch = pitch * pi / 180.0;
	xyaw   = yaw   * pi / 180.0;

	VRotationMatrix(xroll,xpitch,xyaw,rot);

	for (ii=0; ii<3; ii++) {
	  for (jj=0; jj<3; jj++) {
	    irot[ii][jj] = (int)VRint(rot[ii][jj]);
	  }
	}

	for (l=0; l<nn; l++) {
	  x = cmap[l];
	  y = rmap[l];
	  z = bmap[l];

	  b = (int) (irot[0][0] * z + irot[0][1] * y + irot[0][2] * x);
	  r = (int) (irot[1][0] * z + irot[1][1] * y + irot[1][2] * x);
	  c = (int) (irot[2][0] * z + irot[2][1] * y + irot[2][2] * x);

	  if (b < -2 || r < -2 || c < -2) continue;
	  if (b > 2 || r > 2 || c > 2) continue;

	  VPixel(alpha,n,0,l,VSByte) = b;
	  VPixel(alpha,n,1,l,VSByte) = r;
	  VPixel(alpha,n,2,l,VSByte) = c;
	}
	n++;
	yaw += 90;
      }
      pitch += 90;
    }
    roll += 90;
  }

  return alpha;
}


void
GenerateMaps()
{
  int b0_on[7] = {1,1,1,1,1,2,0};
  int r0_on[7] = {1,0,1,2,1,1,1};
  int c0_on[7] = {1,1,0,1,2,1,1};

  int b0_off[3] = { 0,-1, 0};
  int r0_off[3] = {-1, 0, 0};
  int c0_off[3] = { 0, 0,-1};

  int b1_on[7] = {0,1,1,1,1, 1,2};
  int r1_on[7] = {1,0,1,2,1, 1,1};
  int c1_on[7] = {0,0,0,0,1,-1,0};

  int b1_off[2] = {-1, 0};
  int r1_off[2] = { 0,-1};
  int c1_off[2] = { 0, 0};

  int b2_on[6] = {1, 1, 1, 1, 1, 2};
  int r2_on[6] = {0,-1, 0, 1, 0, 0};
  int c2_on[6] = {0, 0,-1, 0, 1, 0};

  int b2_off[1] = {-1};
  int r2_off[1] = { 0};
  int c2_off[1] = { 0};

  alpha0_pos = GenAlpha(b0_on,r0_on,c0_on,7);
  alpha1_pos = GenAlpha(b1_on,r1_on,c1_on,7);
  alpha2_pos = GenAlpha(b2_on,r2_on,c2_on,6);

  alpha0_neg = GenAlpha(b0_off,r0_off,c0_off,3);
  alpha1_neg = GenAlpha(b1_off,r1_off,c1_off,2);
  alpha2_neg = GenAlpha(b2_off,r2_off,c2_off,1);
}




VBoolean
Alpha_3d(VImage image,VImage alpha_pos,VImage alpha_neg,int band,int row,int col)
{
  int b,r,c,n1,n2,b0,r0,c0,angle,l;
  int nangles = 64;

  n1 = VImageNColumns(alpha_pos);
  n2 = VImageNColumns(alpha_neg);

  for (angle = 0; angle < nangles; angle++) {

    for (l=0; l<n1; l++) {
      b0 = VPixel(alpha_pos,angle,0,l,VSByte);
      r0 = VPixel(alpha_pos,angle,1,l,VSByte);
      c0 = VPixel(alpha_pos,angle,2,l,VSByte);

      b = b0 + band;
      r = r0 + row;
      c = c0 + col;
      if (VPixel(image,b,r,c,VBit) != 1) goto next;
    }

    for (l=0; l<n2; l++) {
      b0 = VPixel(alpha_neg,angle,0,l,VSByte);
      r0 = VPixel(alpha_neg,angle,1,l,VSByte);
      c0 = VPixel(alpha_neg,angle,2,l,VSByte);

      b = b0 + band;
      r = r0 + row;
      c = c0 + col;
      if (VPixel(image,b,r,c,VBit) != 0) goto next;
    }
    return TRUE;

  next: ;
  }
  return FALSE;
}




VBoolean
Beta0_3d(VImage image,int band,int row,int col)
{
  int b,r,c,u0,u1,u2,u3;

  for (b=band-1; b <= band; b++) {
    for (r=row-1; r <= row; r++) {
      for (c=col-1; c <= col; c++) {

	u0 = VPixel(image,b,r,c,VBit);

	u1 = VPixel(image,b+1,r+1,c,VBit);
	if (u0 == u1) {
	  u2 = VPixel(image,b+1,r,c,VBit);
	  u3 = VPixel(image,b,r+1,c,VBit);
	  if (u2 == u3 && u3 != u0) return FALSE;
	}

	u1 = VPixel(image,b,r+1,c+1,VBit);
	if (u0 == u1) {
	  u2 = VPixel(image,b,r+1,c,VBit);
	  u3 = VPixel(image,b,r,c+1,VBit);
	  if (u2 == u3 && u3 != u0) return FALSE;
	}

	u1 = VPixel(image,b+1,r,c+1,VBit);
	if (u0 == u1) {
	  u2 = VPixel(image,b+1,r,c,VBit);
	  u3 = VPixel(image,b,r,c+1,VBit);
	  if (u2 == u3 && u3 != u0) return FALSE;
	}
      }
    }
  }

  return TRUE;
}


VBoolean
Beta1_3d(VImage image,int band,int row,int col)
{
  int b,r,c,u0,u1,u2,u3;

  for (b=band-1; b <= band; b++) {
    for (r=row-1; r <= row; r++) {
      for (c=col-1; c <= col; c++) {

	u0 = VPixel(image,b,r,c,VBit);
	u1 = VPixel(image,b+1,r+1,c+1,VBit);
	if (u0 != u1) continue;

	u2 = VPixel(image,b+1,r,c,VBit);
	u3 = VPixel(image,b,r+1,c+1,VBit);
	if (u2 == u3 && u0 != u3) return FALSE;

	u2 = VPixel(image,b,r+1,c,VBit);
	u3 = VPixel(image,b+1,r,c+1,VBit);
	if (u2 == u3 && u0 != u3) return FALSE;

	u2 = VPixel(image,b,r,c+1,VBit);
	u3 = VPixel(image,b+1,r+1,c,VBit);
	if (u2 == u3 && u0 != u3) return FALSE;
      }
    }
  }
  return TRUE;
}


/*
VBoolean
Beta1_3d(VImage image,int band,int row,int col)
{
  int b,r,c,u0,u1,u2,u3;

  for (b=band-1; b <= band; b++) {
    for (r=row-1; r <= row; r++) {
      for (c=col-1; c <= col; c++) {

	u0 = VPixel(image,b,r,c,VBit);
	u1 = VPixel(image,b,r,c+1,VBit);
	if (u0 != u1) {
	  u2 = VPixel(image,b+1,r+1,c,VBit);
	  u3 = VPixel(image,b+1,r+1,c+1,VBit);
	  if (u2 != u3 && u0 == u3) return FALSE;
	}

	u0 = VPixel(image,b,r,c,VBit);
	u1 = VPixel(image,b,r+1,c,VBit);
	if (u0 != u1) {
	  u2 = VPixel(image,b+1,r,c+1,VBit);
	  u3 = VPixel(image,b+1,r+1,c+1,VBit);
	  if (u2 != u3 && u0 == u3) return FALSE;
	}

	u0 = VPixel(image,b,r,c,VBit);
	u1 = VPixel(image,b+1,r,c,VBit);
	if (u0 != u1) {
	  u2 = VPixel(image,b,r+1,c+1,VBit);
	  u3 = VPixel(image,b+1,r+1,c+1,VBit);
	  if (u2 != u3 && u0 == u3) return FALSE;
	}
      }
    }
  }

  return TRUE;
}
*/

