LISA onesample test (vlisa_onesample)
=======================================


The program **vlisa_onesample** performs a one-sample t-test on a group of images
including a correction for multiple comparisons using the LISA algorithm [2018_Lohmann]_.
The input is a list of 3D images resulting from a 1st level computation.
For ease of use, the input images may be specified using wildcards as shown below.
The output is a map thresholded such that FDR < alpha for every voxel. The default is alpha=0.05.
The resulting image shows (1-FDR) so that larger values indicate higher significance.

Note that a region-of-interest mask is required. The mask should exclude non-brain voxels,
and it may cover the entire brain. In the example below, the mask is in the file "braimmask.nii".


::

   vlisa_onesample -in images_*.v -mask brainmask.nii -out result.v



Note that this program also accepts input images in Nifti format ("images_*.nii" or "images_*.nii.gz"), 
but the output is always in vista format.
To convert the output to the Nifti format, use the following command:


::


  vnifti -in result.v -out result.nii





Parameters of 'vlisa_onesample':
```````````````````````````````````

    -help    Prints usage information.
    -in      Input files.
    -out     Output file.
    -alpha   FDR significance level. Default: 0.05
    -perm    Number of permutations. Default: 5000
    -mask    Region of interest mask.
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

.. [2018_Lohmann] Lohmann G., Stelzer J., Lacosse E., Kumar V.J., Mueller K., Kuehn E., Grodd W., Scheffler K. (2018). LISA improves statistical analysis for fMRI. Nature Communications 9:4014. `(link) <https://www.nature.com/articles/s41467-018-06304-z>`_
