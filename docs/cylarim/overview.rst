
Overview
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
voxels with a user-specified radius (typically 3mm).

After the cylinders have been constructed, vcylarim applies a general linear model (GLM) with
three Gaussian basis functions representing deep, middle and superficial activations.
The GLM coefficients are estimated using a pseudoinverse, and can then be used to identify
layer specific activations within each cylinder.

The output of **vcylarim** is a file that contains the estimated GLM-coefficients.
The first three images contain the GLM-coefficients estimated for the deep, middle and superficial layers.
The fourth image contains the intercept values. The GLM coefficients are estimated separately for each cylinder.
In most voxels, several cylinders overlap. In these voxels, the GLM coefficients of the overlapping cylinders are averaged.

Subsequent analysis of the output is performed using **vcylarim_stats**.
This program generates maps that visualize various contrasts between the beta-images,
as defined by the user.

A typical workflow is shown below.






Example workflow
```````````````````````

 ::

   vrim -in segmentation.v -out rim.v

   vmetric -in rim.v -out metric.v
 
   vcylarim -in zmap.v -out cylbeta.v -rim rim.v -metric metric.v -radius 3
   
   vcylarim_stats -in cylbeta.v -out result1.v -type middle
   
   vcylarim_stats -in cylbeta.v -out result2.v -type notmiddle
 


