
vcylarim
--------------------------------------

The program **vcylarim** is the core tool of the Cylarim package.
It requires three primary input images: an **activation map**, a **rim image**
that delineates gray matter (GM) and its boundaries with white matter (WM) and cerebrospinal fluid (CSF),
and a **cortical depth map**. The depth map normalizes cortical depth values to a range of :math:`[0,1]`,
where '0' represents the GM/WM border and '1' signifies the GM/CSF border.
Additionally, a user-defined mask can be provided to restrict the analysis to specific brain regions.
Note that the rim image must use the following codes: 
1:(CSF/GM boundary), 2:(WM/GM boundary), 3:(GM).


An example of a typical call to **vcylarim** is shown below.
In this example, ``zmap.v`` is an fMRI activation map, ``rim.v`` is the image encoding the cortical rim,
and ``metric.v`` is the file that encodes cortical depths. Cylinders of a 2mm radius will be used.

.. code-block::

   vcylarim -in zmap.v -metric metric.v -rim rim.v -radius 2 -nperm 1000 -3bins true -peak true -out outfilename

In this command, Cylarim is instructed to perform two separate analyses: **'3bins'** where vcylarim will
use three depth strata ('deep','middle','superficial'), and conduct pairwise t-tests in which these strata will be
compared against each other. The other analysis  **'peak'** yields a 4D image with two volumes where the first volume
will contain valleys of the cortical depth profiles, and the second volumes contains profile peaks.

Specifically, the above call generates the following output files:

* ``outfilename_3bins_coeff.v``: A 4D image where the three volumes represent the median values of the three depth strata,
  spatially resolved within each image.
* ``outfilename_3bins_zvals.v``: A 4D file containing the statistical maps derived from the 1000 random permutations,
  with three volumes ("deep-middle", "deep-superficial", "middle-superficial").
* ``outfilename_peak_coeff.v``: A 4D file with two volumes, one showing depth profile valleys and the other showing depth profile peaks.

These output files can be visualized directly using Lipsia's own viewer **vini**.
Alternatively, they can be converted to NIfTI format using the program **vnifti** for use with other standard imaging software.

If the option '-3bins true' is set to 'true', then '-nperm' specifies the the number of random permutations will be performed
to obtain z-values for the tests to compare the three layers (deep, middle, superficial) against each other.
During the permutations, a null distribution is estimated by randomly permuting the cortical depth values.
If "-nperm" is set to zero, then no permutations are done, and crude approximations to these z-values are returned.

Other types of analyses can be computed as well. For example, by setting  **-concave true**, a 3D image is computed
which shows the fraction of points per cylinder where the depth profile is "concave up".

A complete list of analysis options is given below.


Further Analyses
----------------

The output files from **vcylarim** can be further processed to perform additional analyses.
For example, to perform a conjunction analysis, the following command using the program **vcylarim_probe** can be used:

.. code-block::

   vcylarim_probe -in outfilename_3bins_zvals.v -out middle_dominant.v -type top_m

Here, the parameter ``-type top_m`` instructs the software to identify predominant middle layer activations using a conjunction analysis.
For further information, see the documentation of **vcylarim_probe**.



Parameters of 'vcylarim'
````````````````````````````````
::

  -help     Prints usage information.
  -in       Input file.
  -out      Output file.
  -mask     Region of interest mask (Optional)
  -metric   Metric image (required)
  -rim      Rim image (required)
  -radius   Cylinder radius in mm. Default: 2
  -equivol  Whether to apply equivolume correction [ true | false ]. Default: false
  -nperm    Number of permutations. Default: 1000
  -3bins    Compute three depth strata (bins) [ true | false ]. Default: false
  -nbins    Compute N depth strata (bins). Default: 0
  -peak     Compute peaks and valleys [ true | false ]. Default: false
  -concave  Compute concave up. Default: false
  -linear   Compute linear slopes [ true | false ]. Default: false
  -R2       Compute model fits [ true | false ]. Default: false
  -maxabs   Compute max absolute z-values [ true | false ]. Default: false
  -seed     Seed for random number generator. Default: 5555
  -j        Number of processors to use, '0' to use all. Default: 0



.. index:: cylarim
