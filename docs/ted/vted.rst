
Task-based edge density (TED)
-----------------------------------

The program **vted** performs a task-based connectivity analysis at the voxel level.
It compares two groups of data sets corresponding to two experimental conditions
A and B. As output it yields a list of voxel pairs that show a stronger correlation in A than in B.
The details of the algorithm are described in Lohmann et al (2016) (see reference below).

The program expects two sets of preprocessed functional data sets as input.
It is important that the preprocessing pipeline includes a baseline correction (detrending).
The input files may correspond to sets of trials. In this case,
the program :doc:`vcuttrials <vcuttrials>` must be called beforehand. The program
'vcuttrials' cuts a functional data set into several pieces where each piece represents one trial.

The program 'vted' also needs a region-of-interest mask as input.
This mask must be geometrically compatible with the
functional data and cover the desired portions of the brain (or the entire brain).
In particular, the mask must have the same spatial resolution, the same image matrix size and
the same orientation as the functional data.
Note that the number of voxels in the mask determines the computational burden.

The program 'vted' should be called in two stages.
The first stage yields a list of voxel pairs that show a stronger correlation in A than in B.
The second stage performs a permutation test to produce a null distribution. This null distribution is needed in
order to assess which of the voxel pairs produced in the first stage show a statistically significant
effect. The first stage does not do any permutations (-perm 0), the number of permutations to be used in
the second stage can be specified using the parameter -perm (e.g. -perm 200).

An example calling sequence is shown below:



Example:
````````````

 ::

   vted -in1 A*.v -in2 B*.v -mask mask.v -perm   0 -q 0.99 -hist realhist.txt -out edgelist.v
   vted -in1 A*.v -in2 B*.v -mask mask.v -perm 100 -q 0.99 -hist nullhist.txt
   vtedfdr -real realhist.txt -null nullhist.txt -out fdr.txt -alpha 0.05
   vhubness -in edgelist.v -out image.v -min 0.3

The first stage has two outputs: a file "edgelist.v" containing a list of candidate edges (voxels pairs) and a histogram file "realhist.txt", which is later used for statistical inference. It also produces a histogram file "realhist.txt"
The second stage produces only one output named "nullhist.txt".
Finally, to asses statistical significance, the program 'vtedfdr' must be called.

The initial threshold is set to "-q 0.99" which means that only the top one percent of all
voxel pairs are considered for subsequent processing.
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



Parameters of 'vted'
`````````````````````````

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


.. index:: ted
