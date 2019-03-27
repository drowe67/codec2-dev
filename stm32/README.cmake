Quickstart
------------

Your system must be able to compile C code natively and have the prerequisites for compiling codec2-dev installed.
See the main INSTALL and README for more information.

In addition, you must have a working gcc-arm-none-eabi toolchain installed (see Readme.txt)

Out-of-source builds is enforced which prevents
accidental pollution of the source repository.

Make a build directory somewhere, ~/build/stm32 will be used in this example:

  $ mkdir -p ~/build/stm32
  $ cd ~/build/stm32


The STM32 Standard Peripheral Library is required. The download`requires a registration on the STM website.
Save the zip file somewhere safe and then extract it anywhere you like. You will have to tell cmake where the
unzipped library is by giving the variable PERIPHLIBDIR the location of top level directory, e.g. for version 1.8.0
this is STM32F4xx_DSP_StdPeriph_Lib_V1.8.0


Configure the build system by running cmake and pointing to the svn checkout
of the codec2-dev/stm32 directory and the StdPeriph lib:

  $ cmake /path/to/codec2-dev/stm32 -DCMAKE_TOOLCHAIN_FILE=/path/to/codec2-dev/stm32/cmake/STM32_Toolchain.cmake -DPERIPHLIBDIR=/path/to/unzipped/STM32F4xx_DSP_StdPeriph_Lib_Vx.x.x
  $ make

If your arm gcc is not in your path, you can specify the location of your compiler 
tool chain like this (Mind the / at the end of the bin path):

  $ cmake -DARM_GCC_BIN:STRING=/path/to/gcc/bin/  /path/to/codec2-dev/stm32 -DCMAKE_TOOLCHAIN_FILE=/path/to/codec2-dev/stm32/cmake/STM32_Toolchain.cmake -DPERIPHLIBDIR=/path/to/unzipped/STM32F4xx_DSP_StdPeriph_Lib_Vx.x.x

To see all the details during compilation:

  $ make VERBOSE=1
