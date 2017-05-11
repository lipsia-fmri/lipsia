
Eigenvector centrality mapping (ECM)
======================================

The program **vecm** computes eigenvector centrality maps (ECM).
ECM attributes a value to each voxel in the brain 
such that a voxel receives a large value if it is strongly correlated with many 
other nodes that are themselves central within the network. 
Google's PageRank algorithm is a variant of eigenvector centrality.
ECM produces a map in which each voxel has a value that indicates
its centrality. These maps can be used for subsequent statistical tests much like
contrast images in standard GLM analyses.

The program expects a preprocessed functional data set as input.
It is important that the preprocessing
pipeline includes a baseline correction (detrending).
The region of interest mask must be geometrically compatible with the
functional data and cover the desired portions of the brain (or the entire brain). 
In particular, the mask must have the same spatial resolution, the same image matrix size and
the same orientation as the functional data.

ECM requires that the connectivity metric yields only non-negative values.
Since linear correlation values range between -1 and +1, we first need to
rescale these values. There are three possible ways to do this:

1. set all negative correlations to zero.
2. add +1 to all  correlations.
3. take absolute values.
4. evaluate negative correlations only, i.e. ignore all positive correlations, and invert all negative correlations to make them positive.


ECM uses the second option (add +1) as the default. But the other three options
are just as valid.

Note that the number of voxels in the mask determines the computational burden.
On typical computer hardware (8 GByte of main memory), the mask may contain
up to 55,000 voxels. Masks containing more voxels are possible provided
sufficient computer memory exist. Note that the amount of memory needed scales up
quadratically with the number of voxels.
The exact number of voxels in the mask can be obtained
by the program :doc:`volumeinfo <../via/volumeinfo>`.


Example:
``````````

 :: 
 
   vecm -in func.v -out result.v -mask mask.v -first 20 -len 80 -type abs


This call uses a subsequence consisting of 80 time points starting at time point 20 for each voxel
covered by the region of interest mask. It first computes a connectivity matrix that
relates every voxel to every other voxel in the mask, and then determines hubs, i.e.
voxels that are well connected with many other voxels that are themselves well connected.

|

**Reference:**
*Lohmann et al (2010),
"Eigenvector centrality mapping for analyzing connectivity patterns in fMRI data of the human brain."
PLoS ONE 5(4): e10232. doi:10.1371/journal.pone.0010232*



Parameters of 'vecm'
````````````````````````````````

 -help    Prints usage information.
 -in      Input file.
 -out     Output file.
 -mask    Region of interest mask.
 -minval  Signal threshold. Default: 0
 -first   First timestep to use. Default: 0
 -length  Length of time series to use, '0' to use full length. Default: 0
 -type    Type of scaling [ pos | add | abs | neg ]. Default: add
 -j       Number of processors to use, '0' to use all. Default: 10


.. index:: ecm

