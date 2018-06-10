Quickstart 1
------------

After installing the packages requires for your system in README.txt

Out-of-source builds are recommeneded (and often enforced) which prevents
accidental pollution of the source repository.

Make a build directory somewhere, ~/build/stm32 will be used in this example:

  $ mkdir -p ~/build/stm32
  $ cd ~/build/stm32


The STM32 Standard Preipheral Library is required and requires registration to
download. Save the zip file somewhere safe and then extract into the build
directory. Currently CMake looks for version 1.8.0:

  (still in build directory)
  $ unzip /path/to/en.stm32f4_dsp_stdperiph_lib.zip

  $ ls
  STM32F4xx_DSP_StdPeriph_Lib_V1.8.0


Configure the build system by running cmake and pointing to the svn checkout
of the codec2-dev/stm32 directory:

  $ cmake /path/to/codec2-dev/stm32
  $ make

To see all the details during compilation:

  $ make VERBOSE=1


