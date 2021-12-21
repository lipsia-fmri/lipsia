
LIPSIA 3.1.0 (May 13, 2019): fMRI analysis tools
======================================

Lipsia is a collection of tools for the analysis of fMRI data. Its main focus is on new algorithms
such as statistical inference (LISA), Eigenvector centrality mapping (ECM) and network detection in task-fMRI (TED).

The development of the Lipsia open access software package was financially supported by the Horizon 2020/CDS-QUAMRI/634541 project. This support is gratefully acknowledged.

Below, a brief description follows. For further details see `documentation`_.



Installation
```````````````````````
Lipsia currently supports Linux and all other operating systems via Docker. Follow the instructions here: `install`_.


Documentation
```````````````````````
Find the full lipsia documentation here: `documentation`_.


Statistical inference (LISA) in examples:
```````````````````````````````````````````````````

*Onesample test at the 2nd level* (`vlisa_onesample`_). 
Example: the input is a set of contrast maps called "data_*.nii.gz"::

  vlisa_onesample -in data_*.nii.gz -mask mask.nii -out result.v
  vnifti -in result.v -out result.nii


*Twosample test at the 2nd level* (`vlisa_twosample`_). 
Example: input are two sets of contrast maps called "data1_*.nii.gz" and "data2_*.nii.gz"::

  vlisa_twosample -in1 data1_*.nii.gz -in2 data2_*.nii.gz -mask mask.nii -out result.v
  vnifti -in result.v -out result.nii


*Single subject test (1st level)* (`vlisa_prewhitening`_). 
Example: input are two runs acquired in the same session called "run1.nii.gz" and "run2.nii.gz".
Preprocessing should include a correction for baseline drifts!::


  vlisa_prewhitening -in run1.nii.gz run2.nii.gz -design des1.txt des2.txt -mask mask.nii -out result.v 
  vnifti -in result.v -out result.nii



Eigenvector centrality mapping (ECM) in examples:
```````````````````````````````````````````````````

Example: input is an fMRI data set called "data.nii.gz" and a brain mask called "mask.nii.gz".::

  vecm -in data.nii.gz -mask mask.nii.gz -j 0 -out ecm.v
  vnifti -in ecm.v -out ecm.nii




Lipsia file format
```````````````````````````````````````
Lipsia uses its own data format, which is called vista (extension *.v).
Many lipsia programs also accept gzipped files or nifti-files as input (*.v.gz or *.nii.gz).
The output is always in unzipped vista-format.
You can easily convert your nifti data from and to lipsia with the program *vnifti*::

  vnifti -in data.nii -out data.v
  vnifti -in data.nii.gz -out data.v
  vnifti -in result.v -out result.nii


Alternatively, you can import a folder with DICOM files into the vista format::

  vdicom -in dir_dicom



Preprocessing
```````````````````````
The current release contains only a rudimentary set of preprocessing tools.
Preprocessing should therefore be performed beforehand using other software packages.
Note that some lipsia algorithms require that the preprocessing pipeline
contains a removal of baseline drifts.
This step can be done using the lipsia program "vpreprocess" if it was omitted
in the initial preprocessing.



.. _install: INSTALL.rst
.. _documentation: docs/index_github.rst


.. _vlisa_onesample: docs/stats/vlisa_onesample.rst
.. _vlisa_twosample: docs/stats/vlisa_twosample.rst
.. _vlisa_prewhitening: docs/stats/vlisa_prewhitening.rst
