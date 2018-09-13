LISA single subject analysis using precoloring (vlisa_precoloring)
====================================================================


The program **vlisa_precoloring** performs a single-subject analysis
including a correction for multiple comparisons using the LISA algorithm [2018_Lohmann]_.
The program expects one or more preprocessed 4D fMRI files as input.
The preprocessing pipeline must include a correction for baseline drifts.
The files typically represent runs of an experiment acquired within the same session.

The output is a map thresholded such that FDR < alpha for every voxel.
The default is alpha=0.05.
The resulting image shows (1-FDR) so that larger values indicate higher significance.


Temporal autocorrelations are accounted for using the precoloring approach [1995_Worsley]_.

Design files in txt-format are also required as input. The number of design files must match the number of runs.
The design files should not include any hemodynamic modelling.
The design file is a txt-file containing one row per trial. Each row must have four values.
The format of each row is: label, onset in seconds, duration in seconds, amplitude.
The labels are integers starting from 1, and represent the label (event type) of that trial.
For example, if there are three different event types, then their labels should be one of 1,2,3.
For a detailed description see  `designformat`_.

The program 'vlisa_precoloring' automatically adds a column of constant values '1' to account for offsets.
It is hidden from view and need not be specified in the design file.

The user must supply a contrast vector whose dimension matches that of the design files
(excluding the hidden constant column).

Hemodynamic modelling is specified using the parameter '-hemo'.
Possible values are "gamma_0" (canonical function), "gamma_1" (canonical function plus its first derivative),
"gamma_2" (canonical function plus its first and second derivatives), or "gauss" (Gaussian function).

A txt-file containing additional covariates (e.g. motion parameters) may be added as a separate file.
It must have one line per time volume and one column per covariate.
The covariates are not subjected to any form of hemodynamic modelling,
and are exempt from random permutations. These covariates are normalized (subtraction of mean).
These nuisance covariates need not be specified in the contrast vector.


The optional parameter '-minval' can be used to specify a threshold that separates brain from non-brain
voxels (brain mask). Voxels whose values in the first volume are below this threshold are discarded from
further analysis. If no value is given, then a heuristic is used to set this threshold automatically.


**Caution:**
The preprocessing pipeline must include a correction for baseline drifts.



Example:

::

   vlisa_precoloring -in run_*.v -design des_*.txt -contrast 1 0 -1 0
	-hemo gamma_1  -out result.v 



Note that this program also accepts input images in Nifti format ("run_*.nii" or "run_*.nii.gz"), 
but the output is always in vista format.
To convert the output to the Nifti format, use the following command:


::


  vnifti -in result.v -out result.nii





Parameters of 'vlisa_precoloring':
```````````````````````````````````

    -help    Prints usage information.
    -in      Input file(s).
    -out     Output file.
    -design   Design files.
    -nuisance   Nuisance covariates (optional).
    -contrast   Contrast vector.
    -hemo    Hemodynamic model [ gamma_0 | gamma_1 | gamma_2 | gauss ]. Default: gamma_0
    -alpha   FDR significance level. Default: 0.05
    -perm    Number of permutations. Default: 5000
    -minval  Signal threshold.
    -seed    Seed for random number generation. Default: 99402622
    -radius  Bilateral parameter (radius in voxels). Default: 2
    -rvar    Bilateral parameter (radiometric). Default: 2.0
    -svar    Bilateral parameter (spatial). Default: 2
    -filteriterations   Bilateral parameter (number of iterations). Default: 2
    -cleanup  Whether to delete isolated voxels. Default: true
    -j        Number of processors to use, '0' to use all. Default: 0





.. index:: lisa_precoloring



References
^^^^^^^^^^^^^^^^^^^^^^^

.. [1995_Worsley] Worsley, K.J. et al (1995). Analysis of fMRI time-series revisited - again. `(link) <https://www.ncbi.nlm.nih.gov/pubmed/9343600>`_
.. [2018_Lohmann] Lohmann G., Stelzer J., Lacosse E., Kumar V.J., Mueller K., Kuehn E., Grodd W., Scheffler K. (in press). LISA improves statistical analysis for fMRI. Nature Communications. `(link) <http://www.kyb.tuebingen.mpg.de/nc/de/employee/details/lohmann.html>`_





.. _designformat: designformat.rst
