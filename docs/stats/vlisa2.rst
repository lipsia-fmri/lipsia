LISA twosample test (vlisa2)
=============================


The program **vlisa2** performs a two-sample t-test on two groups of images
including a correction for multiple comparisons using the LISA algorithm.
For ease of use, the input images may be specified using wildcards as shown below.
The output is a map thresholded such that FDR < alpha for every voxel.
Voxels that survive the threshold have values of (1-FDR) so that larger values indicate
higher significance.


The type of test can be specified using the option '-test'.
The default ('-test ttest') is a twosample t-test based on pooled variance.
The welch test ('-test welch') is a t-test applicable for unequal variances or unequal sample sizes.
The paired test ('-test paired') can be used if the samples are paired.
If wildcards are used to specify the input files, it is important to make sure that
the pairs match.

::

   vlisa2 -in1 images1_*.v -in2 images2_*.v -out result.v -test ttest



**Reference:**
*Lohmann et al (2017),
"Inflated False Negative Rates Undermine Reproducibility In Task-Based fMRI",
bioRxiv, https://doi.org/10.1101/122788*


Parameters of 'vlisa2':
```````````````````````````````````

    -help    Prints usage information.
    -in1     Input files 1.
    -in2     Input files 2.
    -out     Output file.
    -alpha   FDR significance level. Default: 0.05
    -perm    Number of permutations. Default: 2000
    -test    Type of test to be performed [ ttest | paired | welch ].
    -seed    Seed for random number generation. Default: 99402622
    -radius  Bilateral parameter (radius in voxels). Default: 2
    -rvar    Bilateral parameter (radiometric). Default: 2.0
    -svar    Bilateral parameter (spatial). Default: 2
    -filteriterations   Bilateral parameter (number of iterations). Default: 2
    -cleanup  Whether to remove isolated voxels. Default: true
    -j        Number of processors to use, '0' to use all. Default: 0


.. index:: lisa2
