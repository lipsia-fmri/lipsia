
#include <viaio/VImage.h>
#include <viaio/Volumes.h>

#ifdef __cplusplus
extern "C" {
#endif



/*!
  \struct VoxelList
  \brief data struct for representing voxels.
  \param short <b>b</b> slice (band) index
  \param short <b>r</b> row index
  \param short <b>c</b> column index
  
 \par Author:
 Gabriele Lohmann, MPI-CBS
*/
typedef struct {
  short b;
  short r;
  short c;
} Voxel, *VoxelList;



/*!
  \struct PixelList
  \brief data struct for representing pixels.
  \param short <b>r</b> row index
  \param short <b>c</b> column index
  
 \par Author:
 Gabriele Lohmann, MPI-CBS
*/
typedef struct {
  short r;
  short c;
} Pixel, *PixelList;


/*!
  \struct VPoint
  \brief data struct for representing labelled voxels. The label is a float value. 
  \param short <b>b</b> slice (band) index
  \param short <b>r</b> row index
  \param short <b>c</b> column index
  \param float <b>val</b> label value
  
 \par Author:
 Gabriele Lohmann, MPI-CBS
*/
typedef struct {
  short b;
  short r;
  short c;
  float val;
} VPoint;


/*!
  \struct XPoint
  \brief data struct for representing 3D points.
  \param float <b>x</b> slice (band) index
  \param float <b>y</b> row index
  \param float <b>z</b> column index
  
 \par Author:
 Gabriele Lohmann, MPI-CBS
*/
typedef struct pointStruct{
  float x;
  float y;
  float z;
} XPoint;



/*
** access to a pixel
*/
#define VPixelIndex(image, i) \
   ((VPointer) ((char *) (VImageData(image) + (i) * VPixelSize (image))))

/*
** fast implementation of rounding
*/
#define VRintPos(x) ((int)((x) + 0.5))
#define VRintNeg(x) ((int)((x) - 0.5))
#define VRint(x)  ((x) >= 0 ? VRintPos(x) : VRintNeg(x))

#define VFloor(x) ((x) >= 0 ? (int)((x+0.000001)) : (int)((x) - 0.999999))

#define VCeilPos(x) (VRint((x)) == (x) ? (x) : (int)(x+0.99999))
#define VCeilNeg(x) (VRint((x)) == (x) ? (x) : (int)(x+0.00001))
#define VCeil(x) ((x) >= 0 ? VCeilPos(x) : VCeilNeg(x))


/*
** some other stuff
*/
#define VSqr(x) ((x) * (x))
#define VAbs(x) ((x) > 0 > ? (x) : -(x))


#ifdef __cplusplus
}                               /* extern "C" */
#endif                          /* __cplusplus */

