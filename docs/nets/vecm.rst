
Eigenvector centrality mapping (ECM)
======================================

The program **vecm** computes eigenvector centrality maps (ECM) [2010_Lohmann]_.
ECM attributes a value to each voxel in the brain 
such that a voxel receives a large value if it is strongly correlated with many 
other nodes that are themselves central within the network. 
Google's PageRank algorithm is a variant of eigenvector centrality.
ECM produces a map in which each voxel has a value that indicates
its centrality. These maps can be used for subsequent statistical tests much like
contrast images in standard GLM analyses.

The program expects a preprocessed functional data set as input.
It is very important that the preprocessing
pipeline includes a correction for baseline drifts (detrending).
The region of interest mask must be geometrically compatible with the
functional data and cover the desired portions of the brain (or the entire brain). 
In particular, the mask must have the same spatial resolution, the same image matrix size and
the same orientation as the functional data.

ECM is only guaranteed to produce unique results if the correlation values are non-negative.
Below, several correlation metrics that produce non-negative values are listed. Here, 
*corr* denotes Pearson's linear correlation coefficient. 

0. rlc : Rectified linear unit correlation (RLC) as described in [2018b_Lohmann]_
1. add : add +1 to *corr*.
2. pos : set negative values of *corr* to zero.
3. abs : take absolute value of *corr*.
4. neg : set positive values of *corr* to zero, and invert all negative correlations to make them positive.
5. gauss: standard Gaussian of the Euclidean distance between two time courses


The computational burden depends on the type of metric used. For two metrics (rlc, add),
very fast and memory-efficient implementations are provided. The other metrics require that
the correlation matrix is stored in main memory. The amount of memory needed depends quadratically
on the number of voxels in the ROI mask. For standard PCs, ROI masks containing up to 80,000 voxels
are manageable (around 12 GByte of main memory).

If the number of voxels in the ROI mask is so large that the correlation matrix does not fit into main memory,
there are two options: either use one of the two memory-efficient metrics or 
use matrix projections. Matrix projections for ECM are described in [2018b_Lohmann]_  and can 
be used with *vecm* by setting the option '-project true'.



Example:
``````````

 :: 
 
   vecm -in func.v -out result.v -mask mask.v -first 20 -len 80 -metric rlc


This call uses a subsequence consisting of 80 time points starting at time point 20 for each voxel
covered by the region of interest mask. It computes a connectivity matrix that
relates every voxel to every other voxel in the mask, and then determines hubs, i.e.
voxels that are well connected with many other voxels that are themselves well connected.



References
^^^^^^^^^^^^^^

.. [2010_Lohmann] Lohmann, G. et al (2010), Eigenvector centrality mapping for analyzing connectivity patterns in fMRI data of the human brain. PLoS ONE 5(4): e10232. doi:10.1371/journal.pone.0010232. `(link) <https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0010232>`_

.. [2018b_Lohmann] Lohmann, G. et al (2018), Eigenvector centrality mapping for ultrahigh resolution fMRI data of the human brain.
		   bioRxiv  `(link) <https://www.biorxiv.org/content/early/2018/12/13/494732>`_




Parameters of 'vecm'
````````````````````````````````

 -help    Prints usage information.
 -in      Input file.
 -out     Output file.
 -mask    Region of interest mask.
 -first   First timestep to use. Default: 0
 -length  Length of time series to use, '0' to use full length. Default: 0
 -metric  Type of correlation metric [ rlc | add | pos | abs | neg | corr | gauss]. Default: rlc
 -iterations   Maximum number of power iterations. Default: 20
 -project  Whether to do matrix projection. Default: false
 -seed    Seed for random number generation (only for ECM-project). Default: 99402622
 -j       Number of processors to use, '0' to use all. Default: 0


.. index:: ecm

