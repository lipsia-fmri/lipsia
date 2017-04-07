/*! \file

  3D topological thinning.


  The algorithm is based on
  Tsao's method. However, the order in which voxels are checked for
  deletion is governed by distance. We first compute the distance
  transform of the image to be thinned such that each foreground voxel recieves
  a label indicating its distance from the nearest background voxel.
  We then place all foreground voxels in an array and sort it by
  distance values (for instance using ``heapsort'' or "shellsort").
  The points that recieve the smallest distance values are the first to be
  considered for deletion, where the same deletion criteria as in Tsao's original
  algorithm is used. If no more points at this distance level can be deleted,
  we move on to the next higher distance value, and so on until all distance
  levels have been processed.

\par Reference:
  Y.F. Tsao, K.S. Fu (1981).
  "A parallel thinning algorithm for 3D pictures",
  Computer Graphics Image Proc., Vol. 17, pp.315--331.

\par Author:
  Gabriele Lohmann, MPI-CBS
*/


#include <stdio.h>
#include <stdlib.h>
#include <viaio/Vlib.h>
#include <via/via.h>


extern int VCheckPoint(VImage,int,int,int,int);


/*!
  \brief 3D topological thinning

  \fn VImage VThin3d(VImage src,VImage dest,int nadj)

  \param src   binary input image
  \param dest  binary output image
  \param nadj  neighbourhood type (6 or 26)
*/
VImage
VThin3d(VImage src,VImage dest,int nadj)
{
  int b,r,c,nbands,nrows,ncols,npixels;
  int i,j,n,n0,n1,ndel,mdel;
  int dir,db,dr,dc;
  VImage dt=NULL,tmp=NULL;
  VBit *src_pp;
  float depth,maxdepth,step=0.5;
  VPoint *array;

  nrows  = VImageNRows (src);
  ncols  = VImageNColumns (src);
  nbands = VImageNBands (src);
  npixels = nbands * nrows * ncols;
  if (dest == NULL)
    dest = VCreateImage(nbands,nrows,ncols,VBitRepn);

  for (r=0; r<nrows; r++) {
    for (c=0; c<ncols; c++) {
      VPixel(src,0,r,c,VBit) = 0;
      VPixel(src,nbands-1,r,c,VBit) = 0;
    }
  }

  /*
  ** distance transform
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

  /*
  ** set up array of black points
  */
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
  if (n < 2) return src; /* no more than 1 point in image */


  /* heap sort array of black pixels by distance */
  VPoint_hpsort((unsigned long)(n-1),array); 

  /* copy src to dest */
  dest = VCopyImage(src,NULL,VAllBands);


  /* matrix for marking voxels for deletion */
  tmp = VCreateImage(nbands,nrows,ncols,VBitRepn);
  VFillImage(tmp,VAllBands,0);

  /*
  ** delete simple non-border points
  */
  maxdepth = array[n-1].val;

  mdel = 0;
  for (depth=0; depth <= maxdepth; depth += step) {

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

    ndel = 1;
    while (ndel > 0) {
      ndel = 0;

      for (dir=0; dir<6; dir++) {

	db = 0; dr = 0; dc = 0;
	switch (dir) {
	case 0: dr = -1;  break;   /* north */
	case 1: dr =  1;  break;   /* south */
	case 2: dc = -1;  break;   /* west  */
	case 3: dc =  1;  break;   /* east  */
	case 4: db = -1;  break;   /* above */
	case 5: db =  1;  break;   /* below */
	}

	for (j=n0; j<=n1; j++) {

	  if (array[j].val > depth) continue;

	  b = array[j].b;
	  r = array[j].r;
	  c = array[j].c;
	  if (b < 1 || b >= nbands-1) continue;
	  if (r < 1 || r >= nrows-1) continue;
	  if (c < 1 || c >= ncols-1) continue;

	  /* check if already deleted */
	  if (VPixel(dest,b,r,c,VBit) == 0) continue;

	  /* check if already marked */
	  if (VPixel(tmp,b,r,c,VBit) == 1) continue;

	  /* check if directed border point */
	  if (VPixel(dest, b+db, r+dr, c+dc, VBit) != 0) continue;

	  /* checking plane condition */
	  if (VBorderPoint(dest,b,r,c) == 0 ||
	      VCheckPoint(dest,b,r,c,dir) == 0)  continue;

	  /* topological correctness */
	  if (VSimplePoint(dest,b,r,c,nadj) == 1) {
	    VPixel(tmp,b,r,c,VBit) = 1;  /* mark for deletion */
	    ndel++;
	  }
	}

	for (j=n0; j<=n1; j++) {
	  if (array[j].val > depth) continue;
	  b = array[j].b;
	  r = array[j].r;
	  c = array[j].c;
	  if (VPixel(tmp,b,r,c,VBit) > 0) {  /* delete all marked voxels */
	    VPixel(dest,b,r,c,VBit) = 0;
	    mdel++;
	  }
	}
      }
    }
  }

  VFree(array);
  VDestroyImage(tmp);


  /*  output */
  VCopyImageAttrs (src, dest);
  return dest;
}



