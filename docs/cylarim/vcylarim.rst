
vcylarim
======================================

**vcylarim** is the core tool of the Cylarim package.
It expects as input an activation map, a rim image and a metric image.
They can be generated using the cylarim-programs called **vrim** and **vmetric**.
Optionally, they can also be obtained using the software package LayNii.
Note that the rim image must use the following codes: 
1(CSF/GM boundary), 2:(WM/GM boundary), 3:(cortex interior).

The output of **vcylarim** is a file that contains four images. 
The first three images contain the GLM-coefficients estimated for the deep, middle and superficial layers.
The fourth image contains the intercept values. The GLM coefficients are estimated separately for each cylinder.
In most voxels, several cylinders overlap. In these voxels, the GLM coefficients 
estimated of the overlapping cylinders are averaged.




Example:
``````````

 :: 
 
   vcylarim -in zmap.v -out cylbeta.v -rim rim.v -metric metric.v -mask mask.v -radius 3
   


Here, the input is an activation map (zmap.v), a rim image (rim.v), a metric image (metric.v) and a mask image (mask.v).
The mask is optional. If present, it will restrict computations to some part of the input image. This may help
to reduce computation time. It should be a binary image.  The program **vcylarim_getmask** can help to
obtain useful masks.
The parameter "-radius" specifies the radii of the cylinders in mm.
The parameter "-equivol" specifies whether or not equivolume correction should be applied to each cylinder.
**vcylarim** implements parallel computations uses openmp. The number of processors to use can be specified with 
the parameter "-j". The output (cylbeta.v) contains four images containing the estimated GLM coefficients.
This file should be analyzed in subsequent processing steps using **cylarim_stats**.





Parameters of 'vcylarim'
````````````````````````````````

 -help     Prints usage information.
 -in       Input file.
 -out      Output file.
 -mask     Region of interest mask (Optional)
 -metric   Metric image (required)
 -rim      Rim image (required)
 -radius   Cylinder radius in mm. Default: 2
 -equivol  Whether to apply quivolume correction [ true | false ]. Default: false
 -j        Number of processors to use, '0' to use all. Default: 0



.. index:: cylarim

