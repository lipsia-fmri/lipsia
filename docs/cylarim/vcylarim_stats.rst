vcylarim_stats
=====================

The program  **vcylarim_stats** expects the GLM-coefficients computed by **vcylarim** as its input,
and produces a laminar-specific map. The GLM coefficients corresponding to the three layers
are called dx, mx, sx, respectively.
Below is a list of output options.


 - deep ( dx)
 - middle ( mx )
 - superficial ( sx )
 - xdeep  ( 2dx - mx - sx ) 
 - xmiddle  ( 2mx - dx - sx ) 
 - xsuperficial   ( 2sx - mx - sx ) 



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
