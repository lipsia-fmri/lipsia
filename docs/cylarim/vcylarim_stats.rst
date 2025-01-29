vcylarim_stats
=====================

The program  **vcylarim_stats** expects the GLM-coefficients computed by **vcylarim** as its input,
and produces a laminar-specific map. The GLM coefficients corresponding to the three layers
are called d (deep), m (middle), s (superficial), respectively.
Below is a list of output options.


 - deep ( d)
 - middle ( m )
 - superficial ( s )
 - xdeep  ( 2d - m - s ) 
 - xmiddle  ( 2m - d - s ) 
 - xsuperficial   ( 2s - m - s ) 
 - zdeep  ( conjunction((d-m),(d-s)) ) 
 - zmiddle  ( conjunction((m-d),(m-s)) ) 
 - zsuperficial  ( conjunction((s-m),(s-d)) ) 



Examples:
``````````

 :: 
 
   vcylarim_stats -in cylbeta.v -out result1.v -type deep

   vcylarim_stats -in cylbeta.v -out result2.v -type xmiddle


 

Parameters of 'vcylarim_stats'
````````````````````````````````

-help     Prints usage information.
-in       Input file.
-out      Output file.
-type     Type of output (see list above)



.. index:: cylarim_stats
