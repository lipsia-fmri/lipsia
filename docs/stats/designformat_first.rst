First level design files
============================
First level design files describe the experimental setup within one session in a simple txt format.
They specify onset times and duration of the various experimental events.
They can be used as input into single-subject LISA (`vlisa_precoloring`_ , `vlisa_prewhitening`_)

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


.. _vlisa_precoloring: vlisa_precoloring.rst

.. _vlisa_prewhitening: vlisa_prewhitening.rst



.. index:: design files
