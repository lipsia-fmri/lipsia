vcylarim_roi
---------------------

The **vcylarim_roi** program is designed to plot layer-specific activations within a defined Region of Interest (ROI),
providing an alternative analysis pipeline. Its output is a text file, with each row corresponding to a specific bin.
For each bin, the file details four key metrics: cortical depth, the mean activation,
the standard deviation of activation, and the total number of voxels included.
Input file formats can be either nifti (.nii) or vista (.v).


Example:
``````````

 :: 
 
   vcylarim_roi -in zmap.v -roi roi.nii -metric metric.v -nbins 20 -out list.txt

 

Parameters of 'vcylarim_roi'
````````````````````````````````

-help     Prints usage information.
-in       Input zmap file.
-out      Output txt file.
-metric   Image containing the normalized cortical depths.
-roi      Region of interest. Required.
-nbins    Number of bins. Default: 15


.. index:: cylarim_roi
