vcylarim_getmask
======================================

The program **vcylarim_getmask** can be used to obtain an initial brain mask.
The idea is to exclude brain regions where values in the activation map 
are too small to be relevant for laminar-specific analysis.

Just like the main tool **vcylarim** the program begins by generating cylinders.
It then computes the largest and the smallest value of the activation map
within each cylinder.
In voxels where cylinders overlap, the obtained min/max z-values are averaged.

The program outputs a map of min/max z-values. This map can then be thresholded
to obtain a mask to be used in subsequent calls to **vcylarim**.



Example:
``````````

 :: 
 
   vcylarim_getmask -in zmap.v -out cylmask.v -rim rim.v -metric metric.v -radius 3


 
Parameters of 'vcylarim_getmask'
````````````````````````````````

 -help     Prints usage information.
 -in       Input file.
 -out      Output file.
 -mask     Region of interest mask (Optional)
 -metric   Metric image (required)
 -rim      Rim image (required)
 -radius   Cylinder radius in mm. Default: 2
 -j        Number of processors to use, '0' to use all. Default: 0



.. index:: cylarim_getmask
