/*
** TFCE heuristic
**
** M. Kuhlmann, MPI-KYB, Sept 2015
*/

#include <viaio/VImage.h>
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <via/via.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>

#define ABS(x) ((x) > 0 ? (x) : -(x))

VImage tfce(VImage t_image,VImage map,int nstrata)
{
  int b,r,c,i,j,k;
  size_t nslices = VImageNBands(t_image);
  size_t nrows  = VImageNRows(t_image);
  size_t ncols  = VImageNColumns(t_image);
  double tiny = 1e-16;
  
  /* if pixel map is not given, compute it */
  if(!map) {
    VImage mask = VCreateImage(nslices,nrows,ncols,VBitRepn);
    VFillImage(mask,VAllBands,0);
    
    size_t nvox = 0;
    float u;
    for (b=0; b<nslices; b++) {
      for (r=0; r<nrows; r++) {
	for (c=0; c<ncols; c++) {
	  u = VPixel(t_image,b,r,c,VFloat);
	  if (ABS(u) >= tiny) {
	    VPixel(mask,b,r,c,VBit) = 1;
	    nvox++;
	  }
	}
      }
    }

    if(nvox == 0)
      VError("TFCE: Number of voxels!=0 is zero.");
    /* voxel addresses */
    map = VCreateImage(1,3,nvox,VShortRepn);
    if (map == NULL) VError(" error allocating addr map");
    VFillImage(map,VAllBands,0);
    i = 0;
    for (b=0; b<nslices; b++) {
      for (r=0; r<nrows; r++) {
	for (c=0; c<ncols; c++) {
	  if (VGetPixel(mask,b,r,c) < 0.5) continue;
	  VPixel(map,0,0,i,VShort) = b;
	  VPixel(map,0,1,i,VShort) = r;
	  VPixel(map,0,2,i,VShort) = c;
	  i++;
	}
      }
    }
    VDestroyImage(mask);
  }
  
  /* stratification */
  /* set z strata */
  size_t nvox = VImageNColumns(map);
  double z;
  
  /* get an idea of z-value distribution */
  double zmax=0;
  for (i=0; i<nvox; i++) {
    b = VPixel(map,0,0,i,VShort);
    r = VPixel(map,0,1,i,VShort);
    c = VPixel(map,0,2,i,VShort);
    z = VPixel(t_image,b,r,c,VFloat);
    if (z > zmax) zmax = z;
  }

  /* ini strata */
  double step = zmax / (double)nstrata;
  double *stratum = (double *) VCalloc(nstrata+1,sizeof(double));
  for (j=0; j<=nstrata; j++) stratum[j] = step*j;

  float *clustersize = (float *) VCalloc(nvox,sizeof(float));
  if (!clustersize) VError(" err allocating clustersize");
  float *cluster_cont = (float *) VCalloc(nvox,sizeof(float));
  if (!cluster_cont) VError(" err allocating cluster_cont");
  for (i=0; i<nvox; i++) {
    clustersize[i] = 0;
    cluster_cont[i] = 0;
  }

  /* output, create TFCE image */
  VImage tfce  = VCreateImage(nslices,nrows,ncols,VFloatRepn);
  if (!tfce) VError(" err allocating tfce image");
  VFillImage(tfce,VAllBands,0);  
  
  /* ini images */
  VImage bin_image = VCreateImage(nslices,nrows,ncols,VBitRepn);
  if (!bin_image) VError(" err allocating bin_image image");
  VImage label_image = VCreateImage(nslices,nrows,ncols,VShortRepn);
  if (!label_image) VError(" err allocating label_image image");
  
  /* for each stratum */
  double *strat_sq = (double *)VCalloc(nstrata,sizeof(double));
  for (j=0; j<nstrata; j++)
    strat_sq[j] = stratum[j] * stratum[j];
  for (j=0; j<nstrata; j++) {
    VFillImage(bin_image,VAllBands,0);
    for (i=0; i<nvox; i++) {
      b = VPixel(map,0,0,i,VShort);
      r = VPixel(map,0,1,i,VShort);
      c = VPixel(map,0,2,i,VShort);
      z = VPixel(t_image,b,r,c,VFloat);
      if (z > stratum[j]) {
	VPixel(bin_image,b,r,c,VBit) = 1;
      }
    }

    /* get connected components (clusters) */
    long nl=0;
    VLabelImage3d(bin_image,label_image,26,VShortRepn,&nl);
    if (nl < 1) continue;   /* no voxels in this stratum found */

    /* compute cluster sizes */
    int max_k = 0;
    for (i=0; i<nvox; i++) clustersize[i] = 0;
    for (i=0; i<nvox; i++) {
      b = VPixel(map,0,0,i,VShort);
      r = VPixel(map,0,1,i,VShort);
      c = VPixel(map,0,2,i,VShort);
      k = VPixel(label_image,b,r,c,VShort);
      if (k > 0) clustersize[k]++;
      if (k > max_k) max_k = k;
    }

    /* compute contribution of the cluster to the voxels in it */
    for (i=0; i<=max_k; i++)
      cluster_cont[i] = sqrt(clustersize[i]) * strat_sq[j];
    
    /* add cluster contribution to each voxel */
    for (i=0; i<nvox; i++) {
      b = VPixel(map,0,0,i,VShort);
      r = VPixel(map,0,1,i,VShort);
      c = VPixel(map,0,2,i,VShort);
      k = VPixel(label_image,b,r,c,VShort);
      VPixel(tfce,b,r,c,VFloat) += cluster_cont[k];
    }
  }
  
  /* thickness correction */
  for (i=0; i<nvox; i++) {
    b = VPixel(map,0,0,i,VShort);
    r = VPixel(map,0,1,i,VShort);
    c = VPixel(map,0,2,i,VShort);
    k = VPixel(label_image,b,r,c,VShort);
    VPixel(tfce,b,r,c,VFloat) *= step;
  }
  
  VDestroyImage(bin_image);
  VDestroyImage(label_image);
  VFree(clustersize);
  VFree(cluster_cont);
  VFree(stratum);
  VFree(strat_sq);
    
  return tfce;
}
