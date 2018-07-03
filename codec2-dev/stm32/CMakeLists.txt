#
# stm32f4 Codec2 test programs
#
# CMake configuration contributed by Richard Shaw (KF5OIM)
# Please report questions, comments, problems, or patches to the freetel
# mailing list: https://lists.sourceforge.net/lists/listinfo/freetel-codec2
#
cmake_minimum_required(VERSION 2.8)

###################################################

set(FLOAT_TYPE "hard" CACHE STRING "")

###################################################

# Set default C flags.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -std=gnu11 -mlittle-endian -mthumb -mthumb-interwork -nostartfiles -mcpu=cortex-m4")

add_definitions(-DSTM32F40_41xxx -DCORTEX_M4 -D__EMBEDDED__)

if(FLOAT_TYPE STREQUAL "hard")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsingle-precision-constant -Wdouble-promotion -mfpu=fpv4-sp-d16 -mfloat-abi=hard -D__FPU_PRESENT=1 -D__FPU_USED=1")
    #CFLAGS += -fsingle-precision-constant
else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msoft-float")
endif()

option(VERIFY_OPT "Enable this for dump files to help verify optimization" OFF)
if(VERIFY_OPT)
    add_definitions(-DDUMP)
endif()

option(SEMIHOSTING "Enable for standard arm semihosting." OFF)
if(SEMIHOSTING)
    link_libraries(rdimon)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --specs=rdimon.specs")
else()
    link_libraries(nosys)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --specs=nosys.specs")
endif()

# Set default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
endif()
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    link_libraries(g m)
else()
    link_libraries(c m)
endif()

include(cmake/STM32_Toolchain.cmake)
set(CMAKE_REQUIRED_DEFINITIONS "--specs=nosys.specs")
project(stm32f4 C ASM)

include_directories(../src ../unittest inc)


###################################################

# STM32F4 Standard Peripheral Library

include(cmake/STM32_Lib.cmake)


###################################################

# Codec 2

# Output path is such that #include <codec2/version.h> in codec2.h works
set(CODEC2_VERSION_PATH "${PROJECT_BINARY_DIR}/codec2")
configure_file ("${PROJECT_SOURCE_DIR}/../cmake/version.h.in"
                "${CODEC2_VERSION_PATH}/version.h" )
include_directories(${PROJECT_BINARY_DIR})

set(CODEC2_SRC ../src)
set(CODEC2_SRCS
${CODEC2_SRC}/lpc.c
${CODEC2_SRC}/nlp.c
${CODEC2_SRC}/postfilter.c
${CODEC2_SRC}/sine.c
${CODEC2_SRC}/codec2.c
${CODEC2_SRC}/codec2_fft.c
${CODEC2_SRC}/gp_interleaver.c
${CODEC2_SRC}/interldpc.c
${CODEC2_SRC}/kiss_fft.c
${CODEC2_SRC}/kiss_fftr.c
${CODEC2_SRC}/interp.c
${CODEC2_SRC}/lsp.c
${CODEC2_SRC}/mpdecode_core.c
${CODEC2_SRC}/ofdm.c
${CODEC2_SRC}/phase.c
${CODEC2_SRC}/quantise.c
${CODEC2_SRC}/pack.c
${CODEC2_SRC}/codebook.c
${CODEC2_SRC}/codebookd.c
${CODEC2_SRC}/codebookjvm.c
${CODEC2_SRC}/codebookge.c
${CODEC2_SRC}/dump.c
${CODEC2_SRC}/fdmdv.c
${CODEC2_SRC}/freedv_api.c
${CODEC2_SRC}/varicode.c
${CODEC2_SRC}/golay23.c
${CODEC2_SRC}/fsk.c
${CODEC2_SRC}/fmfsk.c
${CODEC2_SRC}/freedv_vhf_framing.c
${CODEC2_SRC}/freedv_data_channel.c
)


set(FFT_TEST_SRCS 
${DSPLIB}/Examples/arm_fft_bin_example/GCC/arm_fft_bin_data.c
src/fft_test.c
src/startup_stm32f4xx.s
src/system_stm32f4xx.c
src/stm32f4_machdep.c
src/gdb_stdio.c
../src/kiss_fft.c
)

add_executable(fft_test ${FFT_TEST_SRCS})
target_link_libraries(fft_test stm32f4 CMSIS)
elf2bin(fft_test)


# Rule for programming the SM1000
#%.pgm: %.bin
#	$(SUDO) dfu-util -d 0483:df11 -c 1 -i 0 -a 0 -s 0x08000000 -D $<

# Rule for programming the SM1000
#%.pgm: %.bin
#	$(SUDO) dfu-util -d 0483:df11 -c 1 -i 0 -a 0 -s 0x08000000 -D $<

####################################################

set(CODEC2_PROFILE_SRCS
src/codec2_profile.c
src/gdb_stdio.c
src/stm32f4_machdep.c
src/startup_stm32f4xx.s
src/init.c
src/system_stm32f4xx.c
)

list(APPEND CODEC2_PROFILE_SRCS ${CODEC2_SRCS})

add_executable(codec2_profile ${CODEC2_PROFILE_SRCS})
target_link_libraries(codec2_profile stm32f4 CMSIS)
target_compile_definitions(codec2_profile PRIVATE PROFILE)
elf2bin(codec2_profile)

set(DAC_UT_SRCS
src/dac_ut.c
../src/fifo.c
src/stm32f4_dac.c
src/debugblinky.c
src/system_stm32f4xx.c
src/startup_stm32f4xx.s
src/init.c
)

add_executable(dac_ut ${DAC_UT_SRCS})
target_link_libraries(dac_ut stm32f4)
elf2bin(dac_ut)

set(FAST_DAC_UT_SRCS
src/fast_dac_ut.c
../src/fifo.c
src/iir_duc.c
src/gdb_stdio.c
src/stm32f4_dacduc.c
src/debugblinky.c
src/system_stm32f4xx.c
src/startup_stm32f4xx.s
src/init.c)

add_executable(fast_dac_ut ${FAST_DAC_UT_SRCS})
target_link_libraries(fast_dac_ut stm32f4)
target_compile_options(fast_dac_ut PRIVATE "-O3")
elf2bin(fast_dac_ut)

set(ADCDAC_UT_SRCS
src/adcdac_ut.c
../src/fifo.c
src/stm32f4_dac.c
src/stm32f4_adc.c
src/sm1000_leds_switches.c
src/debugblinky.c
src/system_stm32f4xx.c
src/startup_stm32f4xx.s
src/init.c)

add_executable(adcdac_ut ${ADCDAC_UT_SRCS})
target_link_libraries(adcdac_ut stm32f4)
target_compile_options(adcdac_ut PRIVATE "-O0")
elf2bin(adcdac_ut)


set(DAC_PLAY_SRCS
src/dac_play.c
../src/fifo.c
src/gdb_stdio.c
src/stm32f4_dac.c
src/debugblinky.c
src/system_stm32f4xx.c
src/startup_stm32f4xx.s
src/init.c)

add_executable(dac_play ${DAC_PLAY_SRCS})
target_link_libraries(dac_play stm32f4)
target_compile_options(dac_play PRIVATE "-O0")
elf2bin(dac_play)

set(ADC_REC_SRCS
src/adc_rec.c
../src/fifo.c
src/gdb_stdio.c
src/stm32f4_adc.c
src/debugblinky.c
src/system_stm32f4xx.c
src/startup_stm32f4xx.s
src/init.c)

add_executable(adc_rec ${ADC_REC_SRCS})
target_link_libraries(adc_rec stm32f4)
elf2bin(adc_rec)

set(ADC_SD_SRCS
src/adc_sd.c
../src/fifo.c
src/gdb_stdio.c
src/stm32f4_adc.c
src/stm32f4_dac.c
src/debugblinky.c
src/system_stm32f4xx.c
src/startup_stm32f4xx.s
src/init.c)

add_executable(adc_sd ${ADC_SD_SRCS})
target_link_libraries(adc_sd stm32f4)
elf2bin(adc_sd)

set(PWM_UT_SRCS
src/gdb_stdio.c
src/stm32f4_pwm.c
src/system_stm32f4xx.c
src/startup_stm32f4xx.s
src/init.c)

add_executable(pwm_ut ${PWM_UT_SRCS})
target_link_libraries(pwm_ut stm32f4)
elf2bin(pwm_ut)

set(POWER_UT_SRCS
src/power_ut.c
src/gdb_stdio.c
../src/fifo.c
src/stm32f4_adc.c
src/stm32f4_dac.c
src/debugblinky.c
src/system_stm32f4xx.c
src/startup_stm32f4xx.s
src/init.c
src/stm32f4_machdep.c)

list(APPEND POWER_UT_SRCS ${CODEC2_SRCS})

# stm32/src/power_ut.c:124:5: error: too few arguments to function 'adc_open'
# stm32/src/power_ut.c:125:5: error: too few arguments to function 'dac_open'
#add_executable(power_ut ${POWER_UT_SRCS})
#target_link_libraries(power_ut stm32f4)

set(USB_VCP
usb_conf/usb_bsp.c
usb_conf/usbd_desc.c
usb_conf/usbd_usr.c
usb_lib/cdc/usbd_cdc_core.c
usb_lib/cdc/usbd_cdc_vcp.c
usb_lib/core/usbd_core.c
usb_lib/core/usbd_ioreq.c
usb_lib/core/usbd_req.c
usb_lib/otg/usb_core.c
usb_lib/otg/usb_dcd.c
usb_lib/otg/usb_dcd_int.c)

set(USB_VCP_UT
src/usb_vcp_ut.c
src/stm32f4_usb_vcp.c
src/sm1000_leds_switches.c
src/system_stm32f4xx.c
src/startup_stm32f4xx.s
src/init.c)

list(APPEND USB_VCP_UT ${USB_VCP})

add_definitions(-DUSE_USB_OTG_FS -DUSE_ULPI_PHY)
include_directories(usb_conf usb_lib/cdc usb_lib/core usb_lib/otg)

add_executable(usb_vcp_ut ${USB_VCP_UT})
target_link_libraries(usb_vcp_ut stm32f4)
elf2bin(usb_vcp_ut)

set(ADC_REC_USB_SRCS
src/adc_rec_usb.c
../src/fifo.c
src/stm32f4_adc.c
src/stm32f4_usb_vcp.c
src/sm1000_leds_switches.c
src/debugblinky.c
src/system_stm32f4xx.c
src/startup_stm32f4xx.s
src/init.c)

list(APPEND ADC_REC_USB_SRCS ${USB_VCP})

add_executable(adc_rec_usb ${ADC_REC_USB_SRCS})
target_link_libraries(adc_rec_usb stm32f4)
elf2bin(adc_rec_usb)


set(FDMDV_PROFILE_SRCS
src/fdmdv_profile.c
src/gdb_stdio.c
src/system_stm32f4xx.c
src/startup_stm32f4xx.s
src/init.c
src/stm32f4_machdep.c)

list(APPEND FDMDV_PROFILE_SRCS ${CODEC2_SRCS})

# Doesn't compile because M_PITCH is not defined but M_PITCH_S is.
#add_executable(fdmdv_profile ${FDMDV_PROFILE_SRCS})
#target_link_libraries(fdmdv_profile stm32f4)
#target_compile_definitions(fdmdv_profile PRIVATE PROFILE)

set(SM1000_LEDS_SWITCHES_UT_SRCS
src/sm1000_leds_switches_ut.c
src/sm1000_leds_switches.c
src/system_stm32f4xx.c
src/startup_stm32f4xx.s
src/init.c)

add_executable(sm1000_leds_switches_ut ${SM1000_LEDS_SWITCHES_UT_SRCS})
target_link_libraries(sm1000_leds_switches_ut stm32f4)
elf2bin(sm1000_leds_switches_ut)


set(SM1000_SRCS
src/sm1000_main.c
src/tone.c
src/sfx.c
src/sounds.c
src/morse.c
src/menu.c
src/tot.c
src/sm1000_leds_switches.c
../src/fifo.c
src/debugblinky.c
src/system_stm32f4xx.c
src/startup_stm32f4xx.s
src/stm32f4_vrom.c
src/init.c)

list(APPEND SM1000_SRCS ${CODEC2_SRCS})

add_library(stm32f4_dac STATIC src/stm32f4_dac.c)

add_library(stm32f4_adc STATIC src/stm32f4_adc.c)

add_executable(sm1000 ${SM1000_SRCS})
target_link_libraries(sm1000 stm32f4_dac stm32f4_adc stm32f4 CMSIS)
target_compile_options(sm1000 PRIVATE "-O3")
elf2bin(sm1000)

set(FREEDV_TX_PROFILE_SRCS
src/freedv_tx_profile.c
src/stm32f4_machdep.c
src/gdb_stdio.c
src/system_stm32f4xx.c
src/startup_stm32f4xx.s
src/init.c)

list(APPEND FREEDV_TX_PROFILE_SRCS ${CODEC2_SRCS})

add_executable(freedv_tx_profile ${FREEDV_TX_PROFILE_SRCS})
target_link_libraries(freedv_tx_profile stm32f4 CMSIS)
target_compile_definitions(freedv_tx_profile PRIVATE PROFILE)
elf2bin(freedv_tx_profile)


set(FREEDV_RX_PROFILE_SRCS
src/freedv_rx_profile.c
src/stm32f4_machdep.c
src/gdb_stdio.c
src/system_stm32f4xx.c
src/startup_stm32f4xx.s
src/init.c)

list(APPEND FREEDV_RX_PROFILE_SRCS ${CODEC2_SRCS})

add_executable(freedv_rx_profile ${FREEDV_RX_PROFILE_SRCS})
target_link_libraries(freedv_rx_profile stm32f4 CMSIS)
elf2bin(freedv_rx_profile)

set(FDMDV_DUMP_RT_SRCS
src/fdmdv_dump_rt.c
src/sm1000_leds_switches.c
../src/fifo.c
src/debugblinky.c
#gdb_stdio.c
src/system_stm32f4xx.c
src/startup_stm32f4xx.s
src/init.c)

list(APPEND FDMDV_DUMP_RT_SRCS ${CODEC2_SRCS})

# Doesn't build due to FREEDV_NSAMPLES not being defined.
#add_executable(fdmdv_dump_rt ${FDMDV_DUMP_RT_SRCS})
#target_link_libraries(fdmdv_dump_rt stm32f4_dac src/stm32f4_adc stm32f4)
#target_compile_options(fdmdv_dump_rt PRIVATE "-O3")

# ---------------------------------------------------------------------------

set(TUNER_UT_SRCS
src/tuner_ut.c
src/gdb_stdio.c
../src/fifo.c
src/stm32f4_dac.c
src/iir_tuner.c
src/sm1000_leds_switches.c
src/debugblinky.c
src/system_stm32f4xx.c
src/startup_stm32f4xx.s
src/init.c
../src/fm.c)

# this needs to be compiled without the optimiser or ugly things happen
# would be nice to work out why as ISRs need to run fast

add_library(stm32f4_adc_tuner STATIC src/stm32f4_adc_tuner.c)
add_executable(tuner_ut ${TUNER_UT_SRCS})
target_compile_options(tuner_ut PRIVATE "-O3")
target_link_libraries(tuner_ut stm32f4_adc_tuner stm32f4)
elf2bin(tuner_ut)

# -----------------------------------------------------------------------------

set(ADC_SFDR_UT_SRCS
src/adc_sfdr_ut.c
src/gdb_stdio.c
../src/fifo.c
src/iir_tuner.c
src/sm1000_leds_switches.c
src/debugblinky.c
src/system_stm32f4xx.c
src/startup_stm32f4xx.s
src/init.c)

add_executable(adc_sfdr_ut ${ADC_SFDR_UT_SRCS})
target_link_libraries(adc_sfdr_ut stm32f4_adc_tuner stm32f4 CMSIS)
elf2bin(adc_sfdr_ut)

# -----------------------------------------------------------------------------

#set(FM_LODUC_PLAY_SRCS
#src/fm_loduc_play.c
#gdb_stdio.c
#../src/fifo.c
#../src/fm.c
#src/debugblinky.c
#src/system_stm32f4xx.c
#src/startup_stm32f4xx.s
#src/init.c)

#add_library(stm32f4_dacloduc STATIC src/stm32f4_dacloduc.c)

#add_executable(fm_loduc_play ${FM_LODUC_PLAY_SRCS})
#target_link_libraries(fm_loduc_play stm32f4_dacloduc)

# -----------------------------------------------------------------------------

set(SI5351_UT_SRCS
src/si5351_ut.c
src/new_i2c.c
src/si53xx.c
src/system_stm32f4xx.c
src/startup_stm32f4xx.s
src/init.c)

add_executable(si5351_ut ${SI5351_UT_SRCS})
target_link_libraries(si5351_ut stm32f4)
elf2bin(si5351_ut)

# -----------------------------------------------------------------------------

set(MCO_UT_SRCS
src/mco_ut.c
src/tm_stm32f4_mco_output.c
src/tm_stm32f4_gpio.c
src/system_stm32f4xx.c
src/startup_stm32f4xx.s
src/init.c)

add_executable(mco_ut ${MCO_UT_SRCS})
target_link_libraries(mco_ut stm32f4)
elf2bin(mco_ut)

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------

set(SM2000_RXDEMO_SRCS
src/sm2000_rxdemo.c
src/sm1000_leds_switches.c
src/debugblinky.c
src/new_i2c.c
src/si53xx.c
src/stm32f4_dac.c
src/stm32f4_adc.c
../src/fifo.c
src/system_stm32f4xx.c
src/startup_stm32f4xx.s
src/init.c)

list(APPEND SM2000_RXDEMO_SRCS ${CODEC2_SRCS})
add_executable(sm2000_rxdemo ${SM2000_RXDEMO_SRCS})
target_link_libraries(sm2000_rxdemo stm32f4 CMSIS)
elf2bin(sm2000_rxdemo)

# ---------------------------------------------------------------------------------

set(SM2000_STW_SRCS
src/sm2000_stw.c
src/sm1000_leds_switches.c
src/debugblinky.c
src/new_i2c.c
src/si53xx.c
src/stm32f4_dac.c
src/stm32f4_adc.c
../src/fifo.c
src/system_stm32f4xx.c
src/startup_stm32f4xx.s
src/init.c)

list(APPEND SM2000_STW_SRCS ${CODEC2_SRCS} ${USB_VCP})
add_executable(sm2000_stw ${SM2000_STW_SRCS})
target_link_libraries(sm2000_stw stm32f4 CMSIS)
elf2bin(sm2000_stw)

# -----------------------------------------------------------------------------

set(SM2000_ADCDUMP_SRCS
src/sm2000_adc_dump.c
src/sm1000_leds_switches.c
src/debugblinky.c
src/new_i2c.c
src/si53xx.c
src/stm32f4_dac.c
src/stm32f4_usb_vcp.c
src/stm32f4_adc.c
../src/fifo.c
src/system_stm32f4xx.c
src/startup_stm32f4xx.s
src/init.c)

list(APPEND SM2000_ADCDUMP_SRCS ${CODEC2_SRCS} ${USB_VCP})
add_executable(sm2000_adcdump ${SM2000_ADCDUMP_SRCS})
target_link_libraries(sm2000_adcdump stm32f4 CMSIS)
elf2bin(sm2000_adcdump)

# -----------------------------------------------------------------------------

# Objects that require the peripheral library
#src/sm1000_main.o: $(PERIPHLIBDIR)/.unpack
#src/codec2_profile.o: $(PERIPHLIBDIR)/.unpack

# -----------------------------------------------------------------------------
