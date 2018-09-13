LISA twosample test (vlisa_twosample)
=========================================


The program **vlisa_twosample** performs a two-sample t-test on two groups of images
including a correction for multiple comparisons using the LISA algorithm [2018_Lohmann]_.
For ease of use, the input images may be specified using wildcards as shown below.
The output is a map thresholded such that FDR < alpha for every voxel. The default is alpha=0.05.
The resulting image shows (1-FDR) so that larger values indicate higher significance.


The type of test can be specified using the option '-test'.
The default ('-test pooled') is a twosample t-test based on pooled variance.
The welch test ('-test welch') is a t-test applicable for unequal variances or unequal sample sizes.
The paired test ('-test paired') can be used if the samples are paired.
If wildcards are used to specify the input files, it is important to make sure that
the pairs match.

::

   vlisa_twosample -in1 images1_*.v -in2 images2_*.v -out result.v -test welch



Note that this program also accepts input images in Nifti format ("*.nii" or "*.nii.gz"), 
but the output is always in vista format.
To convert the output to the Nifti format, use the following command:


::


  vnifti -in result.v -out result.nii




Parameters of 'vlisa_twosample':
```````````````````````````````````

    -help    Prints usage information.
    -in1     Input files 1.
    -in2     Input files 2.
    -out     Output file.
    -alpha   FDR significance level. Default: 0.05
    -perm    Number of permutations. Default: 5000
    -test    Type of test to be performed [ pooled | paired | welch ]. Default: pooled
    -seed    Seed for random number generation. Default: 99402622
    -radius  Bilateral parameter (radius in voxels). Default: 2
    -rvar    Bilateral parameter (radiometric). Default: 2.0
    -svar    Bilateral parameter (spatial). Default: 2
    -filteriterations   Bilateral parameter (number of iterations). Default: 2
    -cleanup  Whether to remove isolated voxels. Default: true
    -j        Number of processors to use, '0' to use all. Default: 0


.. index:: lisa_twosample


References
^^^^^^^^^^^^^^^^^^^^^^^

.. [2018_Lohmann] Lohmann G., Stelzer J., Lacosse E., Kumar V.J., Mueller K., Kuehn E., Grodd W., Scheffler K. (in press). LISA improves statistical analysis for fMRI. Nature Communications. `(link) <http://www.kyb.tuebingen.mpg.de/nc/de/employee/details/lohmann.html>`_
