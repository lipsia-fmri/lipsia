
False discovery rates in TED
===================================

The program **vtedfdr** performs a statistical inference based on the output produced by the
program :doc:`vted <vted>`.
The details of the algorithm are described in Lohmann et al (2016) (see reference below).

The program expects two histogram files as input representing the permuted and non-permuted
outputs of "vted". It produces as output a txt-file containing the false discovery rates (Fdr)
for very edge density value. The cutoff above which Fdr falls below the predetermined significance
level alpha is also reported.

An example calling sequence is shown below:


Example:
``````````

 ::

   vted -in1 A*.v -in2 B*.v -mask mask.v -perm   0 -q 0.99 -hist realhist.txt -out edgelist.v
	 vted -in1 A*.v -in2 B*.v -mask mask.v -perm 100 -q 0.99 -hist nullhist.txt
	 vtedfdr -real realhist.txt -null nullhist.txt -out fdr.txt -alpha 0.05


The first two calls produce the two histograms needed as input into **vtedfdr**.
The third call to 'vtedfdr' uses the two histogram files as input and produces the txt-file "fdr.txt"
as output. The file "fdr.txt" can be used to determine a cutoff so that edges with edge densities
that exceed this cutoff have a sufficiently low false discovery rate. This cutoff is both reported as terminal output and additionally in the output file fdr.txt.


**Reference:**
Lohmann G, Stelzer J, Zuber V, Buschmann T, Margulies D, et al. (2016):
Task-Related Edge Density (TED) - A New Method for Revealing Dynamic Network Formation in fMRI Data of the Human Brain. PLoS ONE 11(6): e0158185. doi: 10.1371/journal.pone.0158185



Parameters of 'vtedfdr'
````````````````````````````````

 -help    Prints usage information.
 -null    Null histogram.
 -real    Non-permutated real histogram.
 -out     Output txt file containing false discovery rates.
 -alpha   Significance level. Default: 0.05

.. index:: tedfdr
