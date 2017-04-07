
Connectivity concordance mapping (CCM)
=======================================


The program **vccm** computes connectivity concordance maps (CCM).
The main idea is to assign a label to each voxel based on the intersubject reproducibility of
its whole-brain pattern of connectivity. Specifically, we compute the correlations of time-courses
of each voxel with every other voxel for each subject. Voxels whose correlation pattern is
consistent across data sets receive high values. 
Concordance across data sets is measured using either Kendall's W or the overall concordance measure
(OCCC).


Preprocessed functional data sets serve as input to 'vccm'. It is important that the preprocessing
pipeline includes a baseline correction (detrending).
The region of interest mask must be geometrically compatible with the
functional data and cover the desired portions of the brain (or the entire brain). 
In particular, the mask must have the same spatial resolution, the same image matrix size and
the same orientation as the functional data.

Note that the number of voxels in the mask determines the computational burden.
The exact number of voxels in the mask can be obtained
by the program :doc:`volumeinfo <../via/volumeinfo>`.


Example:
``````````

 :: 
 
   vccm -in func1.v func2.v -out result.v -mask mask.v -first 20 -len 80 -type kendall


This call uses computes the concordance between the connectivity profiles of the two
input files. For example, the two input files might represent measurements acquired of
the a subject several weeks apart to monitor recovery after stroke.
Here, connectivity is computed using subsequences consisting of 80 time points starting 
at time point 20 for each voxel covered by the region of interest mask. 




**Reference:**
*Lohmann et al (2012),
"Connectivity Concordance Mapping: A New Tool for Model-Free Analysis of fMRI Data of the Human Brain."
Frontiers in Systems Neuroscience. 2012;6:13. doi:10.3389/fnsys.2012.00013.*



Parameters of 'vccm'
````````````````````````````````

 -help    Prints usage information.
 -in      Input file.
 -out     Output file.
 -mask    Region of interest mask.
 -minval  Signal threshold. Default: 0
 -first   First timestep to use. Default: 0
 -length  Length of time series to use, '0' to use full length. Default: 0
 -type    Concordance metric [ kendall | occc ]. Default: kendall


.. index:: ccm

