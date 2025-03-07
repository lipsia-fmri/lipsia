
Bipartite connectivity mapping (BCM)
--------------------------------------

The program **vbcm** implements a new approach for 
analyzing connectivity in rs-fMRI between two brain regions.
The two regions can be of any size, but should be non-overlapping. 
The main idea is to represent connectivity between the regions as a bipartite graph and analyse
it using bipartite network projections. The advantage of this approach is that it allows to fully
preserve spatial precision in both ROIs. This sets it apart from traditional seed-based connectivity mapping
where the seed ROI must be averaged so that its spatial information is lost.



Example:
^^^^^^^^^^

 :: 
 
   vbcm -in func.v -roi1 roi1.v -roi2 roi2.v -out result.v







Parameters of 'vbcm'
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

 -help    Prints usage information.
 -in      Input file.
 -out     Output file.
 -roi1    ROI 1
 -roi2    ROI 2


.. index:: bcm

References
^^^^^^^^^^^^^^

.. [2018a_Lohmann] Lohmann, G. et al (2018), 6th Biennial Conf on Resting State and Brain Connectivity, Montreal, Canada, 


