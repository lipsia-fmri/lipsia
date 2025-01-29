vcylarim_plot
=====================

The inputs into the program  **vcylarim_plot** are a file containing
the GLM-coefficients computed by **vcylarim**, and also a region of interest (ROI).
The original spatial resolution (prior to upsampling) should also be provided.
The program produces a txt file as output that can be used for plotting laminar-specific results.

The output txt file contains three rows, one for each layer.

The first row begins with an index '0', and shows results for the deep layer.

The second row begins with an index '1', and shows results for the middle layer.

The second row begins with an index '2', and shows results for the superficial layer.

The second column in each row is the mean value of all beta values within the given ROI
of that cortical layer. The third and fourth columns are the lower and upper bounds
of the 95 percent confidence interval. The degrees of freedom required for this purpose
are approximated by estimating the number of non-zero voxels in the ROI (prior to upsampling),
divided by 4, which is generally a very conservative estimate.



Example:
``````````

 :: 
 
   vcylarim_plot -in cylbeta.v -roi roi.nii -resolution 0.8 0.8 0.8 -out list.txt


 
 
Example output:
``````````


.. list-table:: 
   :widths: 20 20 20 20
      
   * - layer ID
     - median
     - lower quartile
     - upper quartile
   * - 0
     - 0.072
     - 0.052
     - 0.081
   * - 1	       
     - 0.104
     - 0.086
     - 0.125    
   * - 2
     - 0.048
     - 0.032
     - 0.059


 

Parameters of 'vcylarim_plot'
````````````````````````````````

-help     Prints usage information.
-in       Input file.
-out      Output txt file.
-roi      Region of interest. Required.



.. index:: cylarim_plot
