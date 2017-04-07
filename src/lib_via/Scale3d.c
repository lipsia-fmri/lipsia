/*! \file
  3D scaling using trilinear interpolation.

Scale a 3D image using trilinear interpolation.


\par Author:
Gabriele Lohmann, MPI-CBS
*/


/* From the Vista library: */
#include <viaio/Vlib.h>
#include <viaio/mu.h>
#include <viaio/os.h>
#include <viaio/VImage.h>

#include <stdio.h>
#include <math.h>


/*!
\fn VImage VTriLinearScale3d (VImage src,VImage dest,int dst_nbands,int dst_nrows,int dst_ncols,
                              float shift[3],float scale[3])
\brief 3D scaling using trilinear interpolation, where Ax+b = y, and A is the scaling matrix.
\param src        input image (any repn)
\param dest       output image (any repn)
\param dst_nbands number of output slices
\param dst_nrows  number of output rows
\param dst_ncols  number of output columns
\param shift[3]   translation vector (band,row,column)
\param scale[3]   scaling vector (band,row,column)
*/
VImage 
VTriLinearScale3d (VImage src,VImage dest,int dst_nbands,int dst_nrows,int dst_ncols,
		   float shift[3],float scale[3])
{
  int   b,r,c;
  float bp,rp,cp;
  int   sx, sy, sz;   /* origin of subcube    */
  float px, py, pz;   /* fractions of subcube */
  float qx, qy, qz;   /* fractions of subcube */
  int   lx, ly, lz;   /* lengths */
  int   ox, oy, oz;   /* offsets */
  float val;
  int   src_nrows,src_ncols,src_nbands;
  float xscale,yscale,zscale;
  VRepnKind repn;


  /* Extract data from source image */
  src_nrows  = VImageNRows(src);
  src_ncols  = VImageNColumns(src);
  src_nbands = VImageNBands(src);
  repn   = VPixelRepn(src);


  dest = VSelectDestImage("VTriLinearScale3d",dest,dst_nbands,dst_nrows,dst_ncols,repn);
  if (! dest) return NULL;
  VFillImage(dest,VAllBands,0);




#define GetValues(type) \
{ \
  type *src_pp; \
  src_pp = (type *) VPixelPtr (src, sz, sy, sx); \
  val += (float) pz * py * px * *src_pp; src_pp += ox; \
  val += (float) pz * py * qx * *src_pp; src_pp += oy; \
  val += (float) pz * qy * px * *src_pp; src_pp += ox; \
  val += (float) pz * qy * qx * *src_pp; src_pp += oz; \
  val += (float) qz * py * px * *src_pp; src_pp += ox; \
  val += (float) qz * py * qx * *src_pp; src_pp += oy; \
  val += (float) qz * qy * px * *src_pp; src_pp += ox; \
  val += (float) qz * qy * qx * *src_pp; \
  VPixel(dest,b,r,c,type) = val; \
}
    

  zscale = 1.0 / scale[2];
  yscale = 1.0 / scale[1];
  xscale = 1.0 / scale[0];

  for (b = 0; b < dst_nbands; b++) {
    bp = zscale * ((float) b - shift[2]);
    if (bp < 0 || bp > src_nbands) continue;

    for (r = 0; r < dst_nrows; r++) {
      rp = yscale * ((float) r - shift[1]);
      if (rp < 0 || rp > src_nrows) continue;

      for (c = 0; c < dst_ncols; c++) {
	cp = xscale * ((float) c - shift[0]);
	if (cp < 0 || cp > src_ncols) continue;

	/* compute origin of subcube */
	sx = (int) (cp);
	sy = (int) (rp);
	sz = (int) (bp);

	/* check subcube */
	if ((sx < -1) || (sx >= src_ncols)) continue;
	if ((sy < -1) || (sy >= src_nrows)) continue;
	if ((sz < -1) || (sz >= src_nbands)) continue;

	/* compute fractions of subcube */
	qx = cp - sx; px = 1 - qx;
	qy = rp - sy; py = 1 - qy;
	qz = bp - sz; pz = 1 - qz;

	/* compute lengths and offsets */
	lx = 1;
	ly = src_ncols;
	lz = src_nrows * src_ncols;
	if (sx == -1) {sx = 0; lx = 0;};
	if (sy == -1) {sy = 0; ly = 0;};
	if (sz == -1) {sz = 0; lz = 0;};
	if (sx == src_ncols  - 1) lx = 0;
	if (sy == src_nrows  - 1) ly = 0;
	if (sz == src_nbands - 1) lz = 0;
	ox = lx;
	oy = ox + ly - 2 * lx;
	oz = oy + lz - 2 * ly;
	val = 0;

	switch(repn) {
	case VShortRepn:
	  GetValues(VShort);
	  break;
	case VUByteRepn:
	  GetValues(VUByte);
	  break;
	case VFloatRepn:
	  GetValues(VFloat);
	  break;
	case VSByteRepn:
	  GetValues(VSByte);
	  break;
	default:
	  VError(" illegal pixel repn");
	}
      }
    }
  }

  VCopyImageAttrs (src, dest);
  return dest;
}
