
Hubness maps resulting from TED
-----------------------------------

The program **vhubness** produces a hubness map resulting from the output of :doc:`vted <vted>`.
A voxel in the hubness map records the number of edges for which this voxel serves as an endpoint.
Voxels in which many edges accumulate may be viewed as hubs in a task-specific network,
and the number of edges meeting in a voxel is a measure of the voxel's hubness.

The details of the algorithm are described in Lohmann et al (2016) (see reference below).

An example calling sequence is shown below:


Example:
``````````

 ::

   vted -in1 A*.v -in2 B*.v -mask mask.v -perm   0 -q 0.99 -hist realhist.txt -out edgelist.v 
   vted -in1 A*.v -in2 B*.v -mask mask.v -perm 100 -q 0.99 -hist nullhist.txt
   vtedfdr -real realhist.txt -null nullhist.txt -out fdr.txt -alpha 0.05
   vhubness -in edgelist.v -out image.v -min 0.3
	 

The first call produces a list of candidate voxel pairs and stores them in the file "edgelist.v".
It also produces a histogram file "realhist.txt" which is later used for statistical inference.
The second call produces a null distribution which is obtained using 200 random permutations.
It yields a histogram file "nullhist.txt".
The third call to 'vtedfdr' uses the two histogram files as input and produces the txt-file "fdr.txt"
as output. The file "fdr.txt" can be used to determine a cutoff so that edges with edge densities
that exceed this cutoff have a sufficiently low false discovery rate.

In our example, edge densities larger than 0.3 have a false discovery rate of less than 0.05.
This cutoff is now used to produce a hubness map using the program 'vhubness'.
This voxel map highlights voxels that serve as an endpoint in at least one of the significant edges.


**Reference:**
Lohmann G, Stelzer J, Zuber V, Buschmann T, Margulies D, et al. (2016): 
Task-Related Edge Density (TED) - A New Method for Revealing Dynamic Network Formation in fMRI Data of the Human Brain. PLoS ONE 11(6): e0158185. doi: 10.1371/journal.pone.0158185 



Parameters of 'vhubness'
````````````````````````````````

 -help    Prints usage information.
 -in1     Input files 1.
 -in2     Input files 2.
 -out     Output file.
 -mask    Region of interest mask.
 -perm    Number of permutations. Default: 0
 -qthreshold  Initial quantile threshold. Default: 0.99
 -histogram    Output histogram filename.
 -first   First timestep to use. Default: 0
 -length  Length of time series to use, '0' to use full length. Default: 0
 -adj     Definition of adjacency [6 | 18 | 26]. Default: 26
 -edgelength   Minimum edge length in voxels. Default: 5
 -type    Type of metric [ SNR | median ]. Default: SNR
 -metric  Correlation metric [ pearson | spearman ]. Default: pearson
 -j       Number of processors to use, '0' to use all. Default: 10


.. index:: hubness

