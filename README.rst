
LIPSIA 3.0: fMRI analysis tools
======================================

* vista file format



Installation
```````````````````````
Follow the instructions here: `install`_.

Documentation
```````````````````````
Find the lipsia documentation here: `documentation`_.


Getting started: 30 seconds to Lipsia
```````````````````````````````````````
Convert fmri data file into lipsia format::

  vnifti -in fmri.nii -out fmri.v


Convert mask data file into lipsia format::

  vnifti -in mask.nii -out mask.v


Compute Eigenvector Centrality Map (ECM)::

  vecm -in fmri.v -mask mask.v -j 0 -out ecm.v



.. _install: INSTALL.rst
.. _documentation: docs/index.rst
