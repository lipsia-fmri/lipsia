
LIPSIA 3.1.0: fMRI analysis tools
======================================

Lipsia is a collection of tools for the analysis of fMRI data. Its main focus is on new algorithms
such as statistical inference (LISA), semi-blind machine learning (SML) and Eigenvector centrality mapping (ECM).

The development of the Lipsia open access software package was financially supported by the Horizon 2020/CDS-QUAMRI/634541 project. 
This support is gratefully acknowledged.

Below, a brief description follows. For further details see `documentation`_.



Installation
```````````````````````
Lipsia currently supports Linux and all other operating systems via Docker. Follow the instructions here: `install`_.


Documentation
```````````````````````
Find the full lipsia documentation here: `documentation`_.



Semi-blind machine learning (SML) in an example:
```````````````````````````````````````````````````

  SML is implemented in the program *vsml*. It expects connectome data for all subjects of the training and the test set.
  It is assumed that the connectomes have been precomputed and exist in csv-format. The first step is to
  convert the connectomes into the lipsia-format. This is done using the lipsia program *vreadconnectome*
  as shown in the example below. Here the training set consists of 400 subjects, the test set has 100
  subjects. 
  
  The information about the target variable of interest (e.g.IQ) must be supplied as a text-file ("IQ_train.txt")
  Likewise, information about supplementary non-imaging information 
  (e.g. educational levels) must be supplied as text-files ("Edu_train.txt", "Edu_test"). Each row of those text-file contains
  the value (e.g. IQ or educational level) as a number.
  The ordering of the rows of the text-files must correspond to the ordering in which
  the connectomes are input into *vsml*.  
  
  The output of *vsml* is a text file showing the predictions of the target variable
  for the subjects of the test set.
  
 
  for i in {1...400}; do
    vreadconnectome -in traindata_${i}.csv -out traindata_${i}.v -ncomponents 100; done
  
  for i in {1...100}; do
    vreadconnectome -in testdata_${i}.csv -out testdata_${i}.v -ncomponents 100; done
  


  vsml -train train_*.v -test test_*.v -ytrain IQ_train.txt -ytest IQ_test.txt -xtrain Edu_train.txt -xtest Edu_test.txt \
      -dimX 800 -npls 10 -nensembles 1000 -seed 12345 -out results.txt




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
