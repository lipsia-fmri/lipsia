vfdr
---------------------

The program  **vfdr** is used for controlling the 
false discovery rate (FDR) based on the Benjamini-Hochberg procedure.

The input image must be a map containing z-values that result from a statistical
analysis procedure such as a standard GLM contrast. In the context of cylarim,
the input into vfdr typically results from permutation testing of laminar profiles.
The output is a 3D image corrected for multiple comparisons using FDR.

The FDR significance level can be adjusted using the parameter **-alpha**. If **-alpha** is 
set to '1', no correction is performed.

If the parameter **-dependence** is set to "true", then the output is modified
to account for spatial dependencies in the data using
the Benjamini-Yekutieli Procedure. In this context of cylarim, this is usually not needed.



Example:
``````````

 :: 
 
   vcylarim -in zmap.v -metric metric.v -rim rim.v -radius 2 -nperm 1000 -3bins true -out outfilename
   
   vcylarim_probe -in outfilename_3bins_zvals.v -out middle_dominant.v -type top_m
 
   vfdr -in middle_dominant.v -out middle_dominant_corrected.v -alpha 0.01



 
Parameters of 'vfdr'
````````````````````````````````

-help     Prints usage information.
-in       Input file.
-out      Output file.
-alpha    FDR significance level. Default: 0.05
-dependence Dependence correction [true | false]. Default: false


.. index:: vfdr
