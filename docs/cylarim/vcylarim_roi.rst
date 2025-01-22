vcylarim_roi
=====================

The inputs into the program  **vcylarim_roi** are a file containing
the GLM-coefficients computed by **vcylarim**, and also a region of interest (ROI).
It produces a txt file as output that can be used for plotting laminar-specific results.

The output txt file contains three rows, one for each layer.

The first row begins with an index '0', and shows results for the deep layer.

The second row begins with an index '1', and shows results for the middle layer.

The second row begins with an index '2', and shows results for the superficial layer.

The second column in each row is the median value of all beta values within the given ROI
of that layer. The third column is the lower quartile, the fourth colum is the upper quartile.



Example:
``````````

 :: 
 
   vcylarim_roi -in cylbeta.v -roi roi.nii -out list.txt


 
 
Example output:
``````````
   0   -0.02029  -0.21332   0.13843
   1    0.13023  -0.05394   0.29205
   2    0.13844   0.01691   0.27366


 

Parameters of 'vcylarim_roi'
````````````````````````````````

-help     Prints usage information.
-in       Input file.
-out      Output txt file.
-roi      Region of interest. Required.



.. index:: cylarim_roi
