LISA generic framework
=======================


The program **vlisa0** implements a generic framework for statistical inference of fMRI data.
The user must supply two files as input. The first file is an uncorrected map
in which each voxel contains some test statistic uncorrected for multiple comparisons.
The second file contains a long list of permuted maps that resulting from random permutations.
These permuted maps are subsequently used to derive a null distribution.
The output is a map that gives the voxel-wise statistical significance of 
an effect (FDR corrected).


::

   vlisa0 -in zmap.v -permutations permfile.v -out corrected.v -alpha 0.05



**Reference:**
*Lohmann et al (2017),
"Inflated False Negative Rates Undermine Reproducibility In Task-Based fMRI",
bioRxiv, https://doi.org/10.1101/122788*


Parameters of 'vlisa0':
```````````````````````````````````

    -help    Prints usage information.
    -in      Input map.    
    -permutations  Input file containing a list of permuted maps.
    -out     Output file.
    -alpha   FDR significance level. Default: 1
    -seed    Seed for random number generation. Default: 99402622
    -radius  Neighbourhood radius in voxels. Default: 2
    -rvar    Bilateral parameter (radiometric). Default: 2.0
    -svar    Bilateral parameter (spatial). Default: 2
    -numiter Number of iterations in bilateral filter. Default: 2
    -centering  Whether to do mode centering. Default: false
    -cleanup  Whether to delete isolated voxels. Default: true
    -fdrfile Name of output fdr txt-file. Default: 
    -j        Number of processors to use, '0' to use all. Default: 0


.. index:: lisa0
