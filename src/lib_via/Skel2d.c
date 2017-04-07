/*! \file
  2D Skeletonization

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

/* From the Vista library: */
#include <viaio/Vlib.h>
#include <viaio/file.h>
#include <viaio/mu.h>
#include <viaio/VImage.h>


extern VBoolean Alpha0_2d(VImage,int,int,int);
extern VBoolean Alpha1_2d(VImage,int,int,int);
extern VBoolean Beta_2d(VImage,int,int,int);
extern VBoolean Border_2d(VImage,int,int,int);


/*!
\fn VImage VSkel2d (VImage src,VImage dest)
\param src   input image (bit repn)
\param dest  output image (bit repn)
*/
VImage
VSkel2d(VImage src,VImage dest)
{
  int i,b,r,c,nbands,nrows,ncols,npixels;
  int ndel=0;
  VBit *tmp_pp,*dst_pp;
  VImage tmp=NULL;

  nrows  = VImageNRows (src);
  ncols  = VImageNColumns (src);
  nbands = VImageNBands (src);
  npixels = nrows*ncols;

  dest = VCopyImage(src,dest,VAllBands);
  tmp = VCreateImage(1,nrows,ncols,VBitRepn);

  for (b=0; b<nbands; b++) {

    ndel = 1;
    while (ndel > 0) {
      ndel = 0;

      VFillImage(tmp,VAllBands,0);

      for (r=2; r<nrows-2; r++) {
	for (c=2; c<ncols-2; c++) {
	  if (VPixel(dest,b,r,c,VBit) == 0) continue;

	  if (Border_2d(dest,b,r,c) == FALSE) continue;
	  
	  /* alpha conditions */
	  if (Alpha1_2d(dest,b,r,c) == FALSE && 
	      Alpha0_2d(dest,b,r,c) == FALSE) continue;

	  if (Beta_2d(dest,b,r,c) == FALSE) continue;

	  /* mark for deletion */
	  VPixel(tmp,0,r,c,VBit) = 1; 
	  ndel++;
	}
      }

      /* perform parallel deletion */
      tmp_pp = (VBit *) VPixelPtr(tmp,0,0,0);
      dst_pp = (VBit *) VPixelPtr(dest,b,0,0);

      for (i=0; i<npixels; i++) {
	if (*tmp_pp == 1) *dst_pp = 0;
	tmp_pp++;
	dst_pp++;
      }

    }
  }
  return dest;
}


VBoolean
Border_2d(VImage image,int band,int row,int col)
{
  int r,c;

  for (r=row-1; r<=row+1; r++) {
    for (c=col-1; c<=col+1; c++) {
      if (VPixel(image,band,r,c,VBit) == 0) return TRUE;
    }
  }
  return FALSE;
}


VBoolean
Alpha1_2d(VImage image,int band,int row,int col)
{
  int i,r=0,c=0,n=4;
  int rowon[4] = {-1,0,0,1};
  int colon[4] = { 1,1,2,1};
  int rowoff =  0;
  int coloff = -1;


  /* angle 0 */
  for (i=0; i<n; i++) {
    c = colon[i] + col;
    r = rowon[i] + row;
    if (VPixel(image,band,r,c,VBit) != 1) goto a90;
  }
  c = coloff + col;
  r = rowoff + row;
  if (VPixel(image,band,r,c,VBit) != 0) goto a90;  
  return TRUE;


  /* angle 90 */
 a90:
  for (i=0; i<n; i++) {
    c =  rowon[i] + col;
    r = -colon[i] + row;
    if (VPixel(image,band,r,c,VBit) != 1) goto a180;
  }
  c =  rowoff + col;
  r = -coloff + row;
  if (VPixel(image,band,r,c,VBit) != 0) goto a180;  
  return TRUE;


  /* angle 180 */
 a180:
  for (i=0; i<n; i++) {
    c = -colon[i] + col;
    r = -rowon[i] + row;
    if (VPixel(image,band,r,c,VBit) != 1) goto a270;
  }
  c = -coloff + col;
  r = -rowoff + row;
  if (VPixel(image,band,r,c,VBit) != 0) goto a270;  
  return TRUE;


  /* angle 270 */
 a270:
  for (i=0; i<n; i++) {
    c = -rowon[i] + col;
    r =  colon[i] + row;
    if (VPixel(image,band,r,c,VBit) != 1) goto ende1;
  }
  c = -rowoff + col;
  r =  coloff + row;
  if (VPixel(image,band,r,c,VBit) != 0) goto ende1;
  return TRUE;


 ende1:
  return FALSE;
}



VBoolean
Alpha0_2d(VImage image,int band,int row,int col)
{
  int i,r,c,n=5,k=2;
  int rowon[5] = {0,1,1,1,2};
  int colon[5] = {1,0,1,2,1};
  int rowoff[2] = {-1, 0};
  int coloff[2] = { 0,-1};


  /* angle 0 */
  for (i=0; i<n; i++) {
    c = colon[i] + col;
    r = rowon[i] + row;
    if (VPixel(image,band,r,c,VBit) != 1) goto aa90;
  }
  for (i=0; i<k; i++) {
    c = coloff[i] + col;
    r = rowoff[i] + row;
    if (VPixel(image,band,r,c,VBit) != 0) goto aa90;
  }
  return TRUE;


  /* angle 90 */
 aa90:
  for (i=0; i<n; i++) {
    c =  rowon[i] + col;
    r = -colon[i] + row;
    if (VPixel(image,band,r,c,VBit) != 1) goto aa180;
  }
  for (i=0; i<k; i++) {
    c =  rowoff[i] + col;
    r = -coloff[i] + row;
    if (VPixel(image,band,r,c,VBit) != 0) goto aa180;  
  }
  return TRUE;


  /* angle 180 */
 aa180:
  for (i=0; i<n; i++) {
    c = -colon[i] + col;
    r = -rowon[i] + row;
    if (VPixel(image,band,r,c,VBit) != 1) goto aa270;
  }
  for (i=0; i<k; i++) {
    c = -coloff[i] + col;
    r = -rowoff[i] + row;
    if (VPixel(image,band,r,c,VBit) != 0) goto aa270;  
  }
  return TRUE;


  /* angle 270 */
 aa270:
  for (i=0; i<n; i++) {
    c = -rowon[i] + col;
    r =  colon[i] + row;
    if (VPixel(image,band,r,c,VBit) != 1) goto ende0;
  }
  for (i=0; i<k; i++) {
    c = -rowoff[i] + col;
    r =  coloff[i] + row;
    if (VPixel(image,band,r,c,VBit) != 0) goto ende0;  
  }
  return TRUE;

 ende0:
  return FALSE;
}


VBoolean
Beta_2d(VImage image,int band,int row,int col)
{
  int u0,u1,u2,u3;

  u0 = VPixel(image,band,row,col,VBit);
  u1 = VPixel(image,band,row,col+1,VBit);
  u2 = VPixel(image,band,row+1,col,VBit);
  u3 = VPixel(image,band,row+1,col+1,VBit);
  if (u0 == 0 && u1 == 1 && u2 == 1 && u3 == 0) return FALSE;
  if (u0 == 1 && u1 == 0 && u2 == 0 && u3 == 1) return FALSE;


  u0 = VPixel(image,band,row-1,col,VBit);
  u1 = VPixel(image,band,row-1,col+1,VBit);
  u2 = VPixel(image,band,row,col,VBit);
  u3 = VPixel(image,band,row,col+1,VBit);
  if (u0 == 0 && u1 == 1 && u2 == 1 && u3 == 0) return FALSE;
  if (u0 == 1 && u1 == 0 && u2 == 0 && u3 == 1) return FALSE;

  u0 = VPixel(image,band,row,col-1,VBit);
  u1 = VPixel(image,band,row,col,VBit);
  u2 = VPixel(image,band,row+1,col-1,VBit);
  u3 = VPixel(image,band,row+1,col,VBit);
  if (u0 == 0 && u1 == 1 && u2 == 1 && u3 == 0) return FALSE;
  if (u0 == 1 && u1 == 0 && u2 == 0 && u3 == 1) return FALSE;

  u0 = VPixel(image,band,row-1,col-1,VBit);
  u1 = VPixel(image,band,row-1,col,VBit);
  u2 = VPixel(image,band,row,col-1,VBit);
  u3 = VPixel(image,band,row,col,VBit);
  if (u0 == 0 && u1 == 1 && u2 == 1 && u3 == 0) return FALSE;
  if (u0 == 1 && u1 == 0 && u2 == 0 && u3 == 1) return FALSE;

  return TRUE;
}
