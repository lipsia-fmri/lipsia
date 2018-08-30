
Design file format
=====================


Statistical evaluations require a specification of the experimental design which must
be supplied by the user as a text file. There are two fundamentally different types of
designs: first level designs and higher-level designs.



First level design files
``````````````````````````
First level design files describe the experimental setup within one session in a simple txt format.
They specify onset times and duration of the various experimental events.
They can be used as input into single-subject LISA (**vlisa_prewhitening**).

The design specification of first level designs must be given as a text file in the following format.
Each line of the file corresponds to one trial. Each line must have 4 entries:
1. event type specified by an integer, 2. onset time in seconds, 3. trial duration in seconds, 4. amplitude.
Comment lines begin with one of the following characters '%', '#', '/'. Here are some examples:

 ::

   % Example 1
   % event     onset      duration   amplitude
       1       18.00        3.00        1.00
       1      150.00        3.00        1.00
       1      168.00        2.00        1.00
       2       24.00        1.00        1.00
       2       90.00        1.00        1.00


This design file describes two event types where event 1 occurs at 18,150 and 168 sec., and event 2 occurs
at 24 and 90 sec. The trial lengths of event '1' are 2 or 3 seconds. Trials of event '2' have a duration of
only one second. Both events have a standard amplitude of 1. 


 ::
 
    % Example 2
    % event     onset      duration   amplitude
        1       18.00        3.00        1.00
        1      150.00        3.00        1.00
        1      168.00        2.00        1.00
        2       24.00        1.00        1.00
        2       90.00        1.00        1.00
        3       24.00        1.00        5.80
        3       90.00        1.00        7.26


This is essentially the same design as in example 1. However, event 2
has an additional *parametric* component which is specified as event '3'.
Note that the event therefore occurs twice, once without and once with
the additional parameter.


Higher level design files
``````````````````````````

Higher-level designs are used to compute higher-level statistics such as two-sample or paired t-tests
so that groups of data sets can be statistically assessed. 


For a higher level analysis,
a design file has to be created that can be used in 'vlisa_2ndlevel'.
For all higher level design files, the *number of rows* must coincide
with the number of input files to 'vlisa_2ndlevel'.

Note that convenience programs exist that can be used instead of the more complicated way of
specifying high-level statistics via 'vlisa_2ndlevel'. These programs are:

 - vlisa_onesample
 - vlisa_twosample


However, higher-level design files are needed in more complex situations. 
Below, some simple design files are listed for tutorial purposes.

One-sample test:

 ::

        % 
        % second level design file, one-sample test
        %
        1
        1
        1
        1
        1
        1
        1
        1

This example design can be used to perform a one-sample t-test for 8 contrast images.
Each row in this file corresponds to one contrast map. 
The contrast vector in this case is '1'.
An easier way to do a one-sample test is to use  **vlisa_onesample**.



Two-sample test:

 ::

        % 
        % second level design file, two-sample test
        %
        1 0 
        1 0  
        1 0 
        1 0 
        1 0 
        1 0 
        1 0 
        0 1 
        0 1 
        0 1 
        0 1 
        0 1 
        0 1 
        0 1 
        
				
				
This example design can be used to perform a two-sample test
for comparing two groups of 7 subjects each (e.g. left handed
vs. right handed). The contrast vector comparing the two groups should be '1 -1'. 
An easier way to do a two-sample test is to use  **vlisa_twosample**.


Paired t-test:

 ::

       % 
       % second level design file, paired test
       %
       1 0 1 0 0 0 0 0 0 0 1
       0 1 1 0 0 0 0 0 0 0 1
       1 0 0 1 0 0 0 0 0 0 1
       0 1 0 1 0 0 0 0 0 0 1
       1 0 0 0 1 0 0 0 0 0 1
       0 1 0 0 1 0 0 0 0 0 1
       1 0 0 0 0 1 0 0 0 0 1
       0 1 0 0 0 1 0 0 0 0 1
       1 0 0 0 0 0 1 0 0 0 1
       0 1 0 0 0 0 1 0 0 0 1
       1 0 0 0 0 0 0 1 0 0 1
       0 1 0 0 0 0 0 1 0 0 1
       1 0 0 0 0 0 0 0 1 0 1
       0 1 0 0 0 0 0 0 1 0 1
       1 0 0 0 0 0 0 0 0 1 1
       0 1 0 0 0 0 0 0 0 1 1




This example design can be used to perform a paired test. It is used
for 8 subjects which have been scanned twice (e.g. in the morning and in the afternoon).
If you are interested in differences within the pairs (e.g. morning vs. afternoon), 
the contrast vector is '1 -1 0 0 0 0 0 0 0 0 0'.  
An easier way to do a paired test is to use  **vlisa_twosample -test paired**.

.. index:: design files
