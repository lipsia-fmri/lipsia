vcylarim_stats
---------------------

The program  **vcylarim_stats** expects the GLM-coefficients computed by **vcylarim** as its input,
and produces a laminar-specific map. The GLM coefficients corresponding to the three layers
are called d (deep), m (middle), s (superficial), respectively.
Below is a list of output options.


 - d   ( deep )
 - m   ( middle )
 - s   ( superficial  )
 - d-m    ( z-stats of twosample test 'deep-middle' )
 - d-s    ( z-stats of twosample test 'deep-superficial' )
 - m-s    ( z-stats of twosample test 'middle-superficial' )
 - m-d    ( z-stats of twosample test 'middle-deep' )
 - s-d    ( z-stats of twosample test 'superficial-deep' )
 - s-m    ( z-stats of twosample test 'superficial-middle' )
 - top_d      ( conjunction, d>m & d>s)
 - top_m      ( conjunction, m>d & m>s )
 - top_s      ( conjunction, s>d & s>m)
 - max_id            ( id of largest beta-value, 1(d), 2(m), 3(s))
 - min_id            ( id of smallest beta-value, 1(d), 2(m), 3(s))
 - max     ( maximum of d,m,s )
 - min     ( minimum of d,m,s )
 - maxabs    ( maximum of absolute values of d,m,s )
   


Examples:
``````````

 :: 
 
   vcylarim_stats -in cylbeta.v -out result1.v -type d

   vcylarim_stats -in cylbeta.v -out result2.v -type top_m
   


 

Parameters of 'vcylarim_stats'
````````````````````````````````

-help     Prints usage information.
-in       Input file.
-out      Output file.
-type     Type of output (see list above)



.. index:: cylarim_stats
