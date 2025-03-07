LISA second-level GLM (vlisa_2ndlevel)
-----------------------------------------------


The program **vlisa_2ndlevel** performs a 2nd-level (group level) inference using the general linear model
including a correction for multiple comparisons using the LISA algorithm [2018_Lohmann]_.
The input is a list of 3D images and a 2nd-level design matrix.
The output is a map thresholded such that FDR < alpha for every voxel. The default is alpha=0.05.
The resulting image shows (1-FDR) so that larger values indicate higher significance.

The design file is a txt-file that follows standard conventions for 2nd-level analyses.
It may for instance be generated using the tool `GLM <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/GLM>`_ of the FSL software package,
alternatively the design textfile can also be written manually.
The number of rows in the design file must match the number of input images.
For a detailed description see `designformat_higher`_.

It is important to make sure that the order in which the input images are
specified on the command line matches the order of the rows in the design file.

In some cases - for example repeated measures designs - permutations must be constrained to ensure exchangeability.
These constraints can be specified using a file containing group IDs (option '-grp').
Each input image is represented by one row of this file and contains
an integer ID (>=1) that represents the group to which the image belongs.
Images within the same group are assumed to be exchangeable.
This file corresponds to the file ".grp" of FSL-Glm.
If no group file is specified then it is assumed that no exchangeability constraints exist.

Optionally, nuisance regressors can be supplied in a separate txt-file using the option '-nuisance'.
Each column in this file corresponds to one regressor. The number of rows must match the number of input images (scans).
The nuisance regressors are exempt from random permutations, and they are not included in the contrast vector.
By default, the mean of each regressor is subtracted ('-demean true').

The parameter '-contrast' is used to specify the contrast vector.
The length of the contrast vector must match the number of columns of the design file.
Each element of the contrast vector corresponds to a column in the design file.

Note that a region-of-interest mask is required. The mask should exclude non-brain voxels,
and it may cover the entire brain. In the example below, the mask is in the file "braimmask.nii".



Example:

::

   vlisa_2ndlevel -in images_*.nii.gz -out result.v -mask brainmask.nii -design design.mat
        -grp design.grp -contrast 1 0 0 0 -nuisance regressors.txt





Parameters of 'vlisa_2ndlevel':
```````````````````````````````````

    -help    Prints usage information.
    -in      Input files.
    -out     Output file.
    -design  Design file.
    -grp     Group Ids needed for exchangeability (optional).
    -contrast  Contrast vector.
    -nuisance  Nuisance regressors (optional).
    -demean  Whether to subtract mean in nuisance regressors. Default: true
    -alpha   FDR significance level. Default: 0.05
    -perm    Number of permutations. Default: 5000
    -mask    Region of interest mask.
    -seed    Seed for random number generation. Default: 99402622
    -radius  Bilateral parameter (radius in voxels). Default: 2
    -rvar    Bilateral parameter (radiometric). Default: 2.0
    -svar    Bilateral parameter (spatial). Default: 2
    -filteriterations   Bilateral parameter (number of iterations). Default: 2
    -cleanup  Whether to delete isolated voxels. Default: true
    -j        Number of processors to use, '0' to use all. Default: 0



.. _designformat_higher: designformat_higher.rst



.. index:: lisa_2ndlevel



References
^^^^^^^^^^^^^^^^^^^^^^^

.. [2018_Lohmann] Lohmann G., Stelzer J., Lacosse E., Kumar V.J., Mueller K., Kuehn E., Grodd W., Scheffler K. (2018). LISA improves statistical analysis for fMRI. Nature Communications 9:4014. `(link) <https://www.nature.com/articles/s41467-018-06304-z>`_
