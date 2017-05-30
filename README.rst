
LIPSIA 3.0: fMRI analysis tools
======================================

Lipsia is a collection of tools for the analysis of fMRI data. Its main focus is on new algorithms
such as TED (Task-related edge density) and ECM (Eigenvector centrality mapping).

The current release contains only a rudimentary set of preprocessing tools.
Preprocessing should therefore be performed beforehand using other software packages.
Note that lipsia connectivity algorithms require that the preprocessing pipeline
contains a removal of baseline drifts. 
This step can be done using the lipsia program "vpreprocess" if it was omitted
in the initial preprocessing.


Installation
```````````````````````
Follow the instructions here: `install`_.

Documentation
```````````````````````
Find the lipsia documentation here: `documentation`_.


Lipsia file format
```````````````````````````````````````
Lipsia uses its own data format, which is called vista (extension *.v*).
You can easily convert your nifti data from and to lipsia with the program *vnifti*::

  vnifti -in data.nii -out data.v
  vnifti -in result.v -out result.nii

Alternatively, you can import a folder with DICOM files into the vista format::

  vdicom -in dir_dicom


Getting started: 30 seconds to Lipsia
```````````````````````````````````````
Convert fmri data file into lipsia format::

  vnifti -in fmri.nii -out fmri.v


Convert rfegion of interest mask file into lipsia format::

  vnifti -in mask.nii -out mask.v


Compute Eigenvector Centrality Map (ECM)::

  vecm -in fmri.v -mask mask.v -j 0 -out ecm.v


	
 
.. _install: INSTALL.rst
.. _documentation: docs/index_github.rst

