
Cylarim overview
======================================


Cylarim is an innovative software package for laminar fMRI analysis that
addresses a major limitation of many existing approaches: their reliance on pre-defined regions of interest (ROIs).
The definition of ROIs often involves arbitrary choices,
making them difficult to justify and prone to bias.
To overcome this problem, Cylarim utilizes small, overlapping cylinders that traverse the cortical ribbon,
enabling a comprehensive mapping of laminar-specific activity across extensive cortical regions.




Initial steps
```````````````````

Cylarim assumes the availability of high-precision segmentations of grey matter, white matter,
and cerebrospinal fluid (CSF) with a spatial resolution of approximately 0.2mm.
Data with lower resolution should be upsampled accordingly.

The initial processing steps are similar to the LayNii software package (Huber et al., 2021).
Specifically, we extract cortical rims delineating the boundaries between CSF and white matter.
Subsequently, we derive a normalized 'metric' image representing cortical depth, scaled to the range [0,1].
Cylrim offers two programs for this purpose: **vrim** and **vmetric**.
The outputs of LayNii may be used as well. Note that it is mandatory that the resulting "rim" image
uses the following codes: 1(CSF-GM boundary), 2(WM-GM boundary), 3(GM interior),
with GM(grey matter), WM(white matter), and CSF(cerebral spinal fluid).





The main steps
`````````````````````

**vcylarim** is the core tool in this package.
It first constructs cortical cylinders by connecting pairs of voxels:
one on the gray matter GM/WM rim and the other on the GM/CSF rim.
Voxel pairs are identified by finding the nearest neighbor on the opposite rim for each voxel.
After excluding redundant pairs, cylinders are defined by dilating the line connecting each pair of
voxels with a user-specified radius (typically 2mm).

After the cylinders have been constructed, vcylarim analysis activation profiles using
various methods. For example, vcylarinm can be instructed to compute three depth strata (deep, middle, superficial),
and perform pairwise t-tests to compare them. Alternatively, it can be instructed to identify
peaks in the cortical profiles, or identify regions where the profile is concave up or down.
The output of **vcylarim** is a set of files depending on the various type of analysis selected by the user.

The output files are 4D images that can be further processed using the program **vcylarim_probe**.
Alternatively, they can be visualized using out tool **vini*, or converted to Nifti-format and analysed using
other imaging software.


A typical workflow is shown below. The final step instructs **vcylarim_probe**
to identify predominant middle layer activations ('-type top_m') using a conjunction analysis.





Example workflow
```````````````````````

 ::

   vrim -in segmentation.v -out rim.v

   vmetric -in rim.v -out metric.v
 
   vcylarim -in zmap.v -rim rim.v -metric metric.v -radius 2 -3bins true -nperm 1000 -out cyl
   
   vcylarim_probe -in cyl_3bins_zvals.v -out result.v -type top_m
   
 


Contents
```````````````````````

.. toctree::
   :maxdepth: 1

   vcylarim
   vcylarim_probe
   vrim
   vmetric
   
