
vcylarim
--------------------------------------

**vcylarim** is the core tool of the Cylarim package.
It expects as input an activation map, a rim image and a metric image.
They can be generated using the cylarim-programs called **vrim** and **vmetric**.
Optionally, they can also be obtained using the software package LayNii.
Note that the rim image must use the following codes: 
1:(CSF/GM boundary), 2:(WM/GM boundary), 3:(cortex interior).

The output of **vcylarim** is a file that contains seven images. 
The first three images contain the GLM-coefficients estimated for the deep, middle and superficial layers.
The fourth image contains the intercept values. The following three images contain
z-values of the laminar contrasts (deep-middle), (deep-superficial), (middle-superficial).
The output of **vcylarim**  is subsequently used by **vcylarim_stats** to compute various laminar statistics.

The GLM coefficients are estimated separately for each cylinder.
In most voxels, several cylinders overlap. In these voxels, the GLM coefficients and results of t-tests between
slabs are averaged.



Example:
``````````

 :: 
 
   vcylarim -in zmap.v -out cylstats.v -rim rim.v -metric metric.v -mask mask.v -type glm -nperm 1000 -radius 3
   

 
Here, the input is an activation map (zmap.v), a rim image (rim.v), a metric image (metric.v) and a mask image (mask.v).
The mask is optional. If present, it will restrict computations to some part of the input image. This may help
to reduce computation time. It should be a binary image.  The program **vcylarim_getmask** can help to
obtain useful masks.
The parameter "-radius" specifies the radii of the cylinders in mm.
The parameter "-equivol" specifies whether or not equivolume correction should be applied to each cylinder.

The type of model to be used can be specified by the parameter "-type". If set to "glm", a GLM is applied.
Alternatively, this parameter can be set to "-type slabs". In this case, no GLM is fitted. Instead
a simple averaging within three slabs (deep,middle,superficial) is done. This is easier and computationally
less involved.

The parameter "-nperm 1000" specifies that permutation tests with 1000 permutations are done. These tests
compare the three layers (deep,middle, superficial) against each other. During the permutations,
a null distribution is estimated by randomly permuting the cortical depth values.
If "-nperm" is set to zero, then no permutations are done, and a crude approximations are returned.

**vcylarim** implements parallel computations uses openmp. The number of processors to use can be specified with 
the parameter "-j". The output (cylbeta.v) contains several images as explained above.
This file can be analyzed in subsequent processing steps using **cylarim_stats**.




Parameters of 'vcylarim'
````````````````````````````````
::

  -help     Prints usage information.
  -in       Input file.
  -out      Output file.
  -mask     Region of interest mask (Optional)
  -metric   Metric image (required)
  -rim      Rim image (required)
  -radius   Cylinder radius in mm. Default: 2
  -equivol  Whether to apply equivolume correction [ true | false ]. Default: false
  -nperm    Number of permutations. Default: 1000
  -type     Type of model [ glm | slabs ]. Default: glm
  -j        Number of processors to use, '0' to use all. Default: 0



.. index:: cylarim
