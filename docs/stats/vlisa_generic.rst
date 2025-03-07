LISA generic framework (vlisa_generic)
----------------------------------------


The program **vlisa_generic** implements a generic framework for statistical inference of fMRI data [2018_Lohmann]_.
This program allows the user to supply his/her own permutation images. This can be used to do statistical inference
in arbitrary scenarios for which no dedicated Lisa-Program exists.

The user must supply two files as input. The first file is an uncorrected map
in which each voxel contains some test statistic uncorrected for multiple comparisons.
The second file is a 4D file containing permuted maps supplied by the user. Each "time point"
corresponds to one permutation.

These permuted maps are subsequently used by Lisa to  derive statistical significance including multiple comparison correction.
The output is a map thresholded such that FDR < alpha for every voxel. The default is alpha=0.05.
The resulting image shows (1-FDR) so that larger values indicate higher significance.


::

   vlisa_generic -in zmap.v -permutations permfile.v -out result.v



Note that this program also accepts input images in Nifti format ("*.nii" or "*.nii.gz").
In this case, the permutations file is a 4D image in which each volume represents a permutation.
The output is in vista format.
To convert the output to the Nifti format, use the following command:


::


  vnifti -in result.v -out result.nii




Parameters of 'vlisa_generic':
```````````````````````````````````

    -help    Prints usage information.
    -in      Input map.
    -permutations  Input file containing permutations.
    -out     Output file.
    -alpha   FDR significance level. Default: 0.05
    -seed    Seed for random number generation. Default: 99402622
    -radius  Bilateral parameter (radius in voxels). Default: 2
    -rvar    Bilateral parameter (radiometric). Default: 2.0
    -svar    Bilateral parameter (spatial). Default: 2
    -filteriterations   Bilateral parameter (number of iterations). Default: 2
    -cleanup  Whether to remove isolated voxels. Default: true
    -j        Number of processors to use, '0' to use all. Default: 0


.. index:: lisa_generic

References
^^^^^^^^^^^^^^^^^^^^^^^

.. [2018_Lohmann] Lohmann G., Stelzer J., Lacosse E., Kumar V.J., Mueller K., Kuehn E., Grodd W., Scheffler K. (2018). LISA improves statistical analysis for fMRI. Nature Communications 9:4014 `(link) <https://www.nature.com/articles/s41467-018-06304-z> `_
