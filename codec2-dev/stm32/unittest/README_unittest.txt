
This is the unittest system for the stm32 implementation of codec2/FreeDV

Objectives:

   It is important to have a robust set of tests for the functionality
   of any software project, especially when there are multiple developers.

   It is easier to test small units of code first as stand alone functions.

   The more automated the test system is the easier it is to run and thus
   the more likely people are to run it.


Directory Structure:

   scripts     Where scripts for this unittest system are found

   lib         Where includable files for this unittest system are found
      /python     python library files
      /octave     octave library files

   test_run    Each test is run in a subdirectory here.


Test Run Scripts:

   The basics of running each test are included as comments in the test's
   source code.  However some test code is used be multiple tests and 
   some tests have several actions to setup their data or check their
   output.  So there are scripts to run these tests.

   For example "run_tst_ofdm_demod_ideal" (in the scripts directory).

   These scripts use some common options:

       --clean     Clean out the test run directory (delete it and re-create)
       --setup     Create input files.
       --run       Run the test
       --check     Create reference data and check test output.

   If none of these options are given, then all of the steps will be 
   run, in order.


# vi:set ts=3 et sts=3:
