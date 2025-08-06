
vcylarim_probe
--------------------------------------

The program **vcylarim_probe** is used to further analyse the output produced by **vcylarim**.
The coefficients corresponding to the three depth strata are called d (deep), m (middle), s (superficial), respectively.
Examples of a typical call to **vcylarim_probe** are shown below.


.. code-block::

   vcylarim_probe -in outfilename_3bins_zvals.v -out middle_dominant.v -type top_m
   vcylarim_probe -in outfilename_nbins_coeff.v -out xmax.v -type max_id
   vcylarim_probe -in outfilename_peak_coeff.v -out peaks.v -volume 1


   
In the first example, the parameter ``-type top_m`` instructs the software to identify predominant middle layer activations
using a conjunction analysis. The second example outputs the code for the largest bin in each voxel. The smallest code
has the value '1' and corresponds to the depth stratum closest to the GM/WM border. The largest code value is 'N'
corresponding to the number of bins chosen by the user. It represents the bin closest to the GM/CSF border.
The third example extracts the volume containing profile peaks. Note that that volume numbers begin at '0', so that
the second volume in the file is addressed as '-volume 1'. The first volume contains valleys (not peaks).
 
Other types of analyses can also be computed as detailed in the list of parameters below.



Parameters of 'vcylarim_probe'
````````````````````````````````
::

  -help     Prints usage information.
  -in       Input file.
  -out      Output file.
  -slice    Slice to be selected. Default: 0
  -xinvert  Whether to invert the selected slice[ true | false ]. Default: false
  -type     Type of output [ slice | top_d | top_m | top_s | bottom_d | bottom_m | bottom_s | max_id | min_id | maxabs ]. Default: slice

.. index:: cylarim_probe
