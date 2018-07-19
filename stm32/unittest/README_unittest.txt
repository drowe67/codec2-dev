!! This unittest structure and tools is an experiment and a proposal at this time !!


This is the unittest system for the stm32 implementation of codec2/FreeDV

Objectives:

   It is important to have a robust set of tests for the functionality
   of any software project, especially when there are multiple developers.

   It is easier to test small units of code first as stand alone functions.

   The more automated the test system is the easier it is to run and thus
   the more likely people are to run it.


Overview:

   This system uses a machine readable description of each test with a 
   common set of tools to simplify creating tests and make running them
   as standardized as possible.

   Hopefully there is enough latitude to support tests which have unique 
   requirements too.  The test description provides for skipping steps or
   substituting custom commands instead of the common parts.

   The basic command is

      run_unittest <name_of_test>

   If the name_of_test is "all", then all defined tests are run.

   Options can select the test steps to run or override default directories.


   NOTE: The process of running st-util and gdb is not yet working reliably,
         so that must be run manually at this time.


Test Description:

   This is a YAML formated file whose top level structure is a mapping
   (also called a dictionary, associateve array, etc.) with the name of
   the test as the key.

   The entry for each step can be of several forms:

     - Missing, this step will be skipped
     - A single string will be used as a shell command (after macro substitution).
     - A list will be used a sequence of shell comamnds.
     - A map will be expected to have at least a "type" entry whose value must
       be one of a recognized set of choices.

   Macros:

       Some macros are defined for file locations

       CODEC2_BIN    This is the location where compiled codec2 utility
                     programs are found.  (codec2-dev/build_linux/src)

       RUN_DIR       Each test will run in a subdir here
       CODEC2_UTST   Where codec2 unitest programs are found
       UTST_BIN      Where scripts and programs for this unittest system are found


   Example Test Description File:

      ---
      tst_ofdm_mod :
         binary    : tst_ofdm_mod.elf
         input     : ${CODEC2_BIN}/ofdm_get_test_bits stm_in.raw -f 1
         reference : ${codec2_bin}/ofdm_mod stm_in.raw ref_mod_out.raw
         test      : { type : manual }
         compare   :
      ...


Directory Structure:

   bin         Where scripts and programs for this unittest system are found

   lib         Where includable files for this unittest system are found
      /python     python library files
      /octave     octave library files

   test_run




# vi:set ts=3 et sts=3:
