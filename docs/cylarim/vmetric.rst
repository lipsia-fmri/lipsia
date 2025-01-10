vmetric
=====================

The program  **vmetric** expects a high-resolution segmentation into grey/white matter and CSF as input.
The spatial resolution should be around 0.2mm (not necessarily isotropic).
Data with lower resolution should be upsampled accordingly.
The segmentation image should encode the tissue classes with some integer values.
The codes for the tissue classes are specified using the parameters "-gm", "-wm", "-csf".

The output is a "metric" image representing cortical depth, scaled to the range [0,1].




Example:
``````````

 :: 
 
   vmetric -in segmentation.v -out metric.v -wm 1 -gm 2 -csf 0



 
Parameters of 'vmetric'
````````````````````````````````

-help     Prints usage information.
-in       Input file.
-out      Output file.
-gm       Code for tissue class 'grey matter'
-wm       Code for tissue class 'white matter'
-csf      Code for tissue class 'CSF'


.. index:: metric
