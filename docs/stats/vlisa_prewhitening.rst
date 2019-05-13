LISA single subject analysis using prewhitening (vlisa_prewhitening)
=====================================================================


The program **vlisa_prewhitening** performs a single-subject analysis
including a correction for multiple comparisons using the LISA algorithm [2018_Lohmann]_.
The program expects one or more preprocessed 4D fMRI files as input.
The preprocessing pipeline must include a correction for baseline drifts.
The files typically represent runs of an experiment acquired within the same session.

The output is a map thresholded such that FDR < alpha for every voxel.
The default is alpha=0.05.
The resulting image shows (1-FDR) so that larger values indicate higher significance.


Temporal autocorrelations are accounted for using
the prewhitening approach [2002_Worsley]_.

Design files in txt-format are required as input. The number of design files must match the number of runs.
The design files should not include any hemodynamic modelling.
The design file is a txt-file containing one row per trial. Each row must have four values.
The format of each row is: label, onset in seconds, duration in seconds, amplitude.
The labels are integers starting from 1, and represent the label (event type) of that trial.
For example, if there are three different event types, then their labels should be one of 1,2,3.
For a detailed description see  `designformat_first`_.

The user must supply a contrast vector whose dimension matches that of the design files.

Hemodynamic modelling is specified using the parameter '-hemo'.
Possible values are "gamma_0" (canonical function), "gamma_1" (canonical function plus its first derivative),
"gamma_2" (canonical function plus its first and second derivatives), or "gauss" (Gaussian function).

The option '-col1' determines whether the hemodynamic modelling procedure should add a 
a first (leftmost) column of the design file containing the constant value '1'.
This column represents a global offset. If set to 'true' (default)
then the contrast vector must contain an additional first entry which is usually set to '0'.

A txt-file containing additional nuisance covariates (e.g. motion parameters) may be added as a separate file.
It must have one line per time volume and one column per covariate.
The covariates are not subjected to any form of hemodynamic modelling,
and are exempt from random permutations.
These nuisance covariates need not be specified in the contrast vector.
Nuisance files from several runs must be concatenated into one single file. The number of lines in this file must match
the total number of volumes of all runs.
By default, the mean of each regressor is subtracted ('-demean true').

Optionally, the design file can be obtained as a txt-file using the option '-plot'.
The txt-file that will be output using this option contains one column for each regressor.
The task regressors are convolved with a hemodynamic model, the nuisance covariates are normalized.
The file can be visualized using various plotting tools (e.g. vdesignplot or gnuplot).

Note that a region-of-interest mask is required. The mask should exclude non-brain voxels,
and it may cover the entire brain. In the example below, the mask is in the file "braimmask.nii".



**Caution:**
The preprocessing pipeline must include a correction for baseline drifts.



Examples:

::

   vlisa_prewhitening -in run_*.v -design des_*.txt -contrast 0 1 -1  -mask brainmask.nii -out result.v -plotdesign X.txt
   vdesignplot -in X.txt


Sometimes, it may be useful to obtain an output image that is not FDR-thresholded at a pre-defined alpha level.
This can be achieved by setting the alpha-parameter to 1, i.e. '-alpha 1.0'.
The resulting image can later be thresholded using the program "vlisa_applythreshold". Example:

::

   vlisa_prewhitening -in run_*.v -design des_*.txt -contrast 0 1 -1  -mask brainmask.nii -out result.v -alpha 1.0
   vlisa_applythreshold -in result.v -out thresholded.v -threshold 0.05

::


Note that 'vlisa_prewhitening' also accepts input images in Nifti format ("run_*.nii" or "run_*.nii.gz"),
but the output is always in vista format.
To convert the output to the Nifti format, use the following command:


::


  vnifti -in result.v -out result.nii





Parameters of 'vlisa_prewhitening':
```````````````````````````````````

    -help    Prints usage information.
    -in      Input file(s).
    -out     Output file.
    -design   Design files.
    -contrast   Contrast vector.
    -nuisance   Nuisance regressors (optional).
    -demean  Whether to subtract mean in nuisance regressors. Default: true
    -plotdesign    Filename for plotting design matrix X (optional).
    -order   Order of AR model. Default: 1
    -hemo    Hemodynamic model [ gamma_0 | gamma_1 | gamma_2 | gauss ]. Default: gamma_0
    -col1    Whether to add a constant first column. Default: true
    -alpha   FDR significance level. Default: 0.05
    -perm    Number of permutations. Default: 5000
    -mask    Region of interest mask.
    -minval  Signal threshold.
    -seed    Seed for random number generation. Default: 99402622
    -radius  Bilateral parameter (radius in voxels). Default: 2
    -rvar    Bilateral parameter (radiometric). Default: 2.0
    -svar    Bilateral parameter (spatial). Default: 2
    -filteriterations   Bilateral parameter (number of iterations). Default: 2
    -cleanup  Whether to delete isolated voxels. Default: true
    -j        Number of processors to use, '0' to use all. Default: 0




References
^^^^^^^^^^^^^^^^^^^^^^^

.. [2002_Worsley] Worsley et al. (2002). A general statistical analysis for fMRI data. NeuroImage 15(1):1-15. `(link) <http://www.ncbi.nlm.nih.gov/pubmed/11771969>`_
.. [2018_Lohmann] Lohmann G., Stelzer J., Lacosse E., Kumar V.J., Mueller K., Kuehn E., Grodd W., Scheffler K. (2018). LISA improves statistical analysis for fMRI. Nature Communications 9:4014. `(link) <https://www.nature.com/articles/s41467-018-06304-z>`_



.. _designformat_first: designformat_first.rst



.. index:: lisa_prewhitening
