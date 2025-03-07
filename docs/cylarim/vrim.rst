vrim
---------------------

The program  **vrim** expects a high-resolution segmentation into grey/white matter and CSF as input.
The spatial resolution should be around 0.2mm (not necessarily isotropic).
Data with lower resolution should be upsampled accordingly.
The segmentation image should encode the tissue classes with some integer values.
The codes for the tissue classes are specified using the parameters "-gm", "-wm", "-csf".

The output is a "rim" image that shows the borders between tissues classes.
Specifically, it encodes the CSF/GM boundary with the value '1',
the WM/GM boundary with '2', and the interior of GM with '3'.



Example:
``````````

 :: 
 
   vrim -in segmentation.v -out rim.v -wm 1 -gm 2 -csf 0



 
Parameters of 'vrim'
````````````````````````````````

-help     Prints usage information.
-in       Input file.
-out      Output file.
-gm       Code for tissue class 'grey matter'
-wm       Code for tissue class 'white matter'
-csf      Code for tissue class 'CSF'


.. index:: rim
