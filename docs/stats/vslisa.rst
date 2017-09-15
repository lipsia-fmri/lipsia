LISA single subject analysis
=======================================


The program **vslisa** performs a single-subject analysis
including a correction for multiple comparisons using the LISA algorithm.
The program expects one or more preprocessed 4D fMRI files as input.
The files typically represent runs of an experiment acquired within the same session.
Design files in txt-format are also required as input. The number of design files must match the number of runs.
They should not include any hemodynamic modelling.
The user must supply a contrast vector whose dimension matches that of the design files.

A txt-file containing additional covariates (e.g. motion parameters) may be added as a separate file.
It must have one line per time volume and one column per covariate. 
The covariates are added to the end of the design matrix
"as is". They are not subjected to any form of hemodynamic modelling,
and are exempt from random permutations.
If such covariates are added, the contrast vector must be extended accordingly.

By default, the global mean is also regressed out. The last value in the contrast vector corresponds
to the global mean, and should generally be set to zero.

The output is a map that gives the voxel-wise statistical significance of 
an effect (FDR corrected).
Examples:

::

   vslisa -in run_*.v -design des_*.txt -contrast 0 1 -1 -out zmap.v -alpha 0.05 -globalmean false

::

In the following example, the global mean is regressed out. Therefore, the contrast vector 
has one additional value.

::

   vslisa -in run_*.v -design des_*.txt -contrast 0 1 -1 0 -out zmap.v -alpha 0.05 -globalmean true

::


**Reference:**
*Lohmann et al (2017),
"Inflated False Negative Rates Undermine Reproducibility In Task-Based fMRI",
bioRxiv, https://doi.org/10.1101/122788*


Parameters of 'vslisa':
```````````````````````````````````

    -help    Prints usage information.
    -in      Input file(s).
    -out     Output file.
    -design  Design files.
    -covariates  Additional covariates (optional).
    -contrast Contrast vector.
    -hemo [ gamma_0 | gamma_1 | gamma_2 | gauss ]. Hemodynamic model. Default: gamma_0
    -alpha   FDR significance level. Default: 1
    -perm    Number of permutations. Default: 0
    -minval  Signal threshold. Default: -1e+12
    -seed    Seed for random number generation. Default: 99402622
    -radius  Neighbourhood radius in voxels. Default: 2
    -rvar    Bilateral parameter (radiometric). Default: 2.0
    -svar    Bilateral parameter (spatial). Default: 2.0
    -numiter Number of iterations in bilateral filter. Default: 2
    -cleanup  Whether to delete isolated voxels. Default: false
    -globalmean  Whether to regress out global mean. Default: true
    -fdrfile  Name of output fdr txt-file. Default: 
    -j        Number of processors to use, '0' to use all. Default: 10


.. index:: slisa
