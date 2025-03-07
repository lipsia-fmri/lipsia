
Cut functional data into trials
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The program **vcuttrials** takes a preprocessed functional data set as input and cuts it up
into small pieces along trial borders. This step is needed as a preparation for a subsequent
call to **vted**. The programs reads the information about trial onsets and durations
from a design file. The format of the design file is described elsewhere.
The start and end of the pieces relative to the trial onsets are specified using the parameters
'-start' and 'length'. The output temporal resolution is set in seconds using the parameter '-resolution'.

An example calling sequence is shown below:




Example:
``````````

 :: 

   vcuttrials -in func.v -des design.txt -start 0 -length 25 -reso 1.5 -cond 1 -trial 0 -out A_0.v 
   vcuttrials -in func.v -des design.txt -start 0 -length 25 -reso 1.5 -cond 1 -trial 1 -out A_1.v
   vcuttrials -in func.v -des design.txt -start 0 -length 25 -reso 1.5 -cond 1 -trial 2 -out A_2.v
   
   vcuttrials -in func.v -des design.txt -start 0 -length 25 -reso 1.5 -cond 2 -trial 0 -out B_0.v
   vcuttrials -in func.v -des design.txt -start 0 -length 25 -reso 1.5 -cond 2 -trial 1 -out B_1.v
   vcuttrials -in func.v -des design.txt -start 0 -length 25 -reso 1.5 -cond 2 -trial 2 -out B_2.v
 

The above calls produce a total of six output files. The first three represent the first three trials
of condition 1. The last three represent the first three trials of condition 2.
The output files have names "A_*.v" and "B_*.v" corresponding to the two experimental conditions. 


Parameters of 'vcuttrials'
````````````````````````````````

 -help    Prints usage information.
 -in      Input file.
 -out     Output file.
 -design  Design file.
 -cond    Id of experimental condition. Default: 1
 -trial   Id of trial (starts at 0). Default: 0
 -reso    Output temporal resolution in secs. Default: 1
 -start   Start relative to beginning of trial in secs. Default: 0
 -length  Trial length in seconds. Default: 20

.. index:: cuttrials

