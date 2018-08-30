LISA onesample test (vlisa_onesample)
=======================================


The program **vlisa_onesample** performs a one-sample t-test on a group of images
including a correction for multiple comparisons using the LISA algorithm [2018_Lohmann]_.
The input is a list of 3D images resulting from a 1st level computation.
For ease of use, the input images may be specified using wildcards as shown below.
The output is a map thresholded such that FDR < alpha for every voxel. The default is alpha=0.05.
The resulting image shows (1-FDR) so that larger values indicate higher significance.



::

   vlisa_onesample -in images_*.v -out result.v





Parameters of 'vlisa_onesample':
```````````````````````````````````

    -help    Prints usage information.
    -in      Input files.
    -out     Output file.
    -alpha   FDR significance level. Default: 0.05
    -perm    Number of permutations. Default: 5000
    -seed    Seed for random number generation. Default: 99402622
    -radius  Bilateral parameter (radius in voxels). Default: 2
    -rvar    Bilateral parameter (radiometric). Default: 2.0
    -svar    Bilateral parameter (spatial). Default: 2
    -filteriterations   Bilateral parameter (number of iterations). Default: 2
    -cleanup  Whether to remove isolated voxels. Default: true
    -j        Number of processors to use, '0' to use all. Default: 0


.. index:: lisa_onesample

References
^^^^^^^^^^^^^^^^^^^^^^^

.. [2018_Lohmann] Lohmann G., Stelzer J., Lacosse E., Kumar V.J., Mueller K., Kuehn E., Grodd W., Scheffler K. (in press). LISA improves statistical analysis for fMRI. Nature Communications. `(link) <http://www.kyb.tuebingen.mpg.de/nc/de/employee/details/lohmann.html>`_
