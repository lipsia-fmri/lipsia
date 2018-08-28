
LIPSIA 3.0: fMRI analysis tools
======================================

Lipsia is a collection of tools for the analysis of fMRI data. Its main focus is on new algorithms
such as statistical inference (LISA), Eigenvector centrality mapping (ECM) and network detection in task-fMRI (TED).
Below, a brief description follows. For further details see `documentation`_.



Installation
```````````````````````
Lipsia currently supports Linux and MacOS. Follow the instructions here: `install`_.


Documentation
```````````````````````
Find the full lipsia documentation here: `documentation`_.


Statistical inference (LISA) in examples:
```````````````````````````````````````````````````

*Onesample test at the 2nd level* `vlisa_onesample`_: vlisa_onesample
Example: the input is a set of contrast maps called "data_*.nii.gz"::

  vlisa_onesample -in data_*.nii.gz -out result.v
  vnifti -in result.v -out result.nii


*Twosample test at the 2nd level* `vlisa_twosample`_: vlisa_twosample
Example: input are two sets of contrast maps called "data1_*.nii.gz" and "data2_*.nii.gz"::

  vlisa_twosample -in1 data1_*.nii.gz -in2 data2_*.nii.gz -out result.v
  vnifti -in result.v -out result.nii


*Single subject test (1st level)* `vlisa_prewhitening`_: vlisa_prewhitening
Example: input are two runs acquired in the same session called "run1.nii.gz" and "run2.nii.gz".
Preprocessing should include a correction for baseline drifts!::


  vlisa_prewhitening -in run1.nii.gz run2.nii.gz -design des1.txt des2.txt -out result.v
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


.. _vlisa_onesample: stats/vlisa_onesample.rst
.. _vlisa_twosample: stats/vlisa_twosample.rst
.. _vlisa_prewhitening: stats/vlisa_prewhitening.rst
