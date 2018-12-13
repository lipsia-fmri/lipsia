#include "viadata.h"

/* Utilities  */
extern VFloat VReadPixel(VImage,int,int,int); 
extern void   VWritePixel(VImage,int,int,int,VFloat); 
extern VFloat VGetPixelValue (VImage,int); 
extern void   VSetPixelValue (VImage,int,VFloat);
extern void   VQuaternionsRotation(float *,float **);
extern void   VRotationMatrix(double,double,double,double [3][3]);
extern void   VPoint_hpsort(unsigned long,VPoint[]);

extern VImage VContrast(VImage,VImage,VRepnKind,VFloat,VFloat);
extern VImage VContrastUByte(VImage,VImage,VFloat,VFloat);
extern VImage VContrastShort(VImage,VImage,VFloat,VFloat);
extern VImage VContrastShortUByte(VImage, int, int, float, float);
extern VImage VMapImageRange(VImage,VImage,VRepnKind);
extern VImage VMaskRange(VImage,VImage,VFloat,VFloat);
extern VImage VMask(VImage,VImage,VImage);

/* 3D binary morphology */
extern VImage VDilateImage3d (VImage, VImage,VoxelList,int);
extern VImage VErodeImage3d(VImage,VImage,VoxelList,int);
extern VImage VGenSphere3d(VShort);
extern VImage VGenSphere2d(VShort);
extern VoxelList VConvertSE3d(VImage,int *);

/* 3D greylevel morphology  */
extern VImage VGreyDilation3d(VImage,VImage,VImage);
extern VImage VGreyErosion3d(VImage,VImage,VImage);

/* morphological operators by thresholding distance transform */
extern VImage VDTClose(VImage,VImage,VDouble);
extern VImage VDTOpen(VImage,VImage,VDouble);
extern VImage VDTErode(VImage,VImage,VDouble);
extern VImage VDTDilate(VImage,VImage,VDouble);

/* distance transforms */
extern VImage VEuclideanDist3d(VImage,VImage,VRepnKind);
extern VImage VChamferDist3d(VImage,VImage,VRepnKind);
extern VImage VChamferDist2d(VImage,VImage,VBand);
extern VImage VCDT3d (VImage,VImage,VLong,VLong,VRepnKind);

/* connected components */
extern VImage VLabelImage2d(VImage,VImage,long,VRepnKind,long *);
extern VImage VLabelImage3d(VImage,VImage,long,VRepnKind,long *);
extern VImage VSelectBig (VImage,VImage);
extern VImage VDeleteSmall(VImage,VImage,long);


/* topological operators */
extern VImage VBorderImage3d (VImage,VImage);
extern int    VBorderPoint(VImage,int,int,int);
extern int    VSimplePoint(VImage,int,int,int,int);
extern VImage VTopoclass(VImage,VImage);
extern int    VGenusLee (VImage,VShort);
extern VImage VThin3d(VImage,VImage,int);
extern VImage VSkel3d(VImage,VImage);
extern VImage VSkel2d(VImage,VImage);


/* resampling, geometric transformations */
extern VImage VTriLinearSample3d(VImage,VImage,VImage,float,float,float,int,int,int);
extern VImage VNNSample3d(VImage,VImage,VImage,float,float,float,int,int,int);
extern VImage VCubicSplineSample3d(VImage,VImage,VImage,float,float,float,int,int,int);
extern VImage VTriLinearScale3d (VImage,VImage,int,int,int,float[3],float[3]);
extern VImage VNNScale3d (VImage,VImage,int,int,int,float[3],float[3]);
extern VImage VCubicSplineScale3d (VImage,VImage,int,int,int,float[3],float[3]);
extern VImage VBicubicScale2d(VImage,VImage,VFloat);
extern VImage VBiLinearScale2d(VImage,VImage,int,int,float[2],float[2]);
extern VImage VRotateImage2d(VImage,VImage,VBand,double);
extern VImage VRotateImage3d(VImage,VImage,VFloat,VShort);
extern VImage VShearImageX (VImage src, VImage dest, VBand band, double shear);
extern VImage VShearImageY (VImage src, VImage dest, VBand band, double shear);


/* filters */
extern VImage VAniso2d(VImage,VImage,VShort,VShort,VFloat,VFloat);
extern VImage VAniso3d(VImage,VImage,VShort,VShort,VFloat,VFloat);
extern VImage VConvolve3d(VImage,VImage,VImage);
extern VImage VConvolve2d(VImage,VImage,VImage);
extern VImage VConvolveCol(VImage,VImage,VImage);
extern VImage VConvolveRow(VImage,VImage,VImage);
extern VImage VConvolveBand(VImage,VImage,VImage);
extern VImage VLeeImage(VImage,VImage,VLong,VDouble,VLong,VLong);
extern VImage VSmoothImage3d(VImage,VImage,VLong,VLong);
extern VImage VFilterGauss2d(VImage,VImage,double);
extern VImage VFilterGauss3d(VImage,VImage,double *);
extern VImage VFilterBox3d(VImage,VImage,int);
extern VImage VMedianImage3d (VImage,VImage,int,VBoolean);
extern VImage VMedianImage2d (VImage,VImage,int,VBoolean);



/* segmentation, clustering, binarization */
extern VImage VBinarizeImage(VImage,VImage,VDouble,VDouble);
extern VImage VIsodataImage3d(VImage,VImage,VLong,VLong);

/* edge detection and curvature */
extern void   VCanny3d(VImage,int,VImage *,VImage *,VImage *);
extern void   VCanny2d(VImage,int,VImage *,VImage *);
extern void   VDeriche3d(VImage,VFloat,VImage *,VImage *,VImage *);
extern void   VDeriche2d(VImage,VFloat,VImage *,VImage *);
extern VImage VMagnitude3d(VImage,VImage,VImage,VImage);
extern VImage VMagnitude2d(VImage,VImage,VImage);
extern VImage VNonmaxSuppression(VImage,VImage,VImage,VImage);
extern VImage VCurvature(VImage,VImage,VFloat,VLong,VBoolean);

/* operations on volumes */
extern Volumes VImage2Volumes(VImage src);
extern VImage  Volumes2Image(Volumes,VRepnKind repn);
extern VImage  Volume2Bin(Volume);
extern Volume  VBin2Volume(VImage);
extern double  VolumeRadius(Volume,double *);
extern VBoolean VolumeInside(Volume,long,long,long);
extern VTrack  VolumeGetTrack(Volume,long,long,long);
extern Volume  VCopyVolume(Volume);


/* shape features: moments,center of gravity,size */
extern void   VolumeCentroid (Volume,double *);
extern double VolumeMoment (Volume,double *,long,long,long);
extern void   VBinCentroid (VImage, double *);
extern double VBinMoment (VImage,double *,long,long,long);
extern long   VolumeSize(Volume);
extern long   VBinSize(VImage);
