vcylarim_seed
---------------------

The program **vcylarim_seed** is designed to allow inspection of individual cylinders.
The cylinder is specified using the parameter **-seed**. The program locates the cylinder
that is closest to the given seed. The seed must be specified using voxel coordinates (x,y,z)
(not MNI coordinates). 
The output image shows the selected cylinder filled with voxel values of the input zmap.
All other voxels of the output image are set to zero.

In addition, txt files are produced depending on the parameter **-type**.

If "-type 3bins", then three output txt files will be produced called "deep.txt",
"middle.txt", "superficial.txt" containing one row per voxel depending on membership
in one of these three layers. Each row contains two values: the normalized cortical depth and the
corresponding z-value of the input zmap.

If "-type chebyshev", only one txt file will be produced. Its name is specified using
the parameter **-report**. It has one row for each voxel in the selected
cylinder. And each row contains three values:  1) the normalized cortical depth,
2) the corresponding z-value of the input zmap, 3) the value of the 
fitted chebyshev polynomial of order 4 at that point.


Input file formats can be either nifti (.nii) or vista (.v). The output image
is in vista-format, and can be converted to nifti using the program *vnifti**.


Example:
``````````

 :: 
 
   vcylarim_seed -in zmap.v -metric metric.v -rim rim.v -rad 2 -seed 20 30 43 -type 3bins -out z.v

 

Parameters of 'vcylarim_seed'
````````````````````````````````

-help     Prints usage information.
-in       Input zmap file.
-out      Output image file.
-report   Output txt file for type 'chebyshev'.
-metric   Image containing the normalized cortical depths.
-rim      Rim image. Required.
-radius   Cylinder radius in mm. Default: 2
-type     Type of output [ 3bins | chebyshev ]. Default: 3bins



.. index:: cylarim_seed
