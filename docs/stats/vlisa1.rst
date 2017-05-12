LISA onesample test
=======================


The program **vlisa1** performs a one-sample t-test on a group of images
including a correction for multiple comparisons using the LISA algorithm.
For ease of use, the input images may be specified using wildcards as shown below.
The output is a map that gives the voxel-wise statistical significance of 
an effect (FDR corrected).


::

   vlisa1 -in images_*.v -out zmap.v -alpha 0.05



**Reference:**
*Lohmann et al (2017),
"Inflated False Negative Rates Undermine Reproducibility In Task-Based fMRI",
bioRxiv, https://doi.org/10.1101/122788*


Parameters of 'vlisa1':
```````````````````````````````````

    -help    Prints usage information.
    -in      Input files.
    -out     Output file.
    -alpha   FDR significance level. Default: 1
    -perm    Number of permutations. Default: 0
    -seed    Seed for random number generation. Default: 99402622
    -radius  Neighbourhood radius in voxels. Default: 2
    -rvar    Bilateral parameter (radiometric). Default: 2.0
    -svar    Bilateral parameter (spatial). Default: 2
    -numiter Number of iterations in bilateral filter. Default: 2
    -cleanup  Whether to delete isolated voxels. Default: false
    -filename Name of output fdr txt-file. Default: 
    -j        Number of processors to use, '0' to use all. Default: 10


.. index:: lisa1
