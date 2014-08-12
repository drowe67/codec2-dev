# Makefile for stm32f4 Codec 2 test programs

###################################################

FLOAT_TYPE=hard

###################################################

BINPATH=~/gcc-arm-none-eabi-4_7-2013q1/bin
CC=$(BINPATH)/arm-none-eabi-gcc
OBJCOPY=$(BINPATH)/arm-none-eabi-objcopy
SIZE=$(BINPATH)/arm-none-eabi-size

###################################################

CFLAGS  = -std=gnu99 -g -Wall -Tstm32_flash.ld -DSTM32F4XX -DCORTEX_M4
CFLAGS += -mlittle-endian -mthumb -mthumb-interwork -nostartfiles -mcpu=cortex-m4

ifeq ($(FLOAT_TYPE), hard)
CFLAGS += -fsingle-precision-constant -Wdouble-promotion
#CFLAGS += -fsingle-precision-constant
CFLAGS += -mfpu=fpv4-sp-d16 -mfloat-abi=hard -D__FPU_PRESENT=1 -D__FPU_USED=1
else
CFLAGS += -msoft-float
endif

###################################################

# Definitions for the STM32F4 Standard Peripheral Library

PERIPHLIBURL    = http://www.st.com/st-web-ui/static/active/en/st_prod_software_internet/resource/technical/software/firmware/
PERIPHLIBZIP    = stm32f4_dsp_stdperiph_lib.zip
PERIPHLIBVER	= V1.1.0
PERIPHLIBNAME	= STM32F4xx_DSP_StdPeriph_Lib
PERIPHLIBDIR	= $(PERIPHLIBNAME)_$(PERIPHLIBVER)
CMSIS		= $(PERIPHLIBDIR)/Libraries/CMSIS
STM32F4LIB	= $(PERIPHLIBDIR)/Libraries/STM32F4xx_StdPeriph_Driver
STM32F4TEMPLATE	= $(PERIPHLIBDIR)/Project/STM32F4xx_StdPeriph_Templates
DSPLIB          = $(PERIPHLIBDIR)/Libraries/CMSIS/DSP_Lib

CFLAGS		+= -DUSE_STDPERIPH_DRIVER -I$(STM32F4LIB)/inc -I$(STM32F4TEMPLATE)
CFLAGS		+= -I$(CMSIS)/Include -I$(CMSIS)/Device/ST/STM32F4xx/Include
CFLAGS		+= -DARM_MATH_CM4

###################################################

# Codec 2

CODEC2_SRC=../src
CODEC2_SRCS=\
$(CODEC2_SRC)/lpc.c \
$(CODEC2_SRC)/nlp.c \
$(CODEC2_SRC)/postfilter.c \
$(CODEC2_SRC)/sine.c \
$(CODEC2_SRC)/codec2.c \
$(CODEC2_SRC)/kiss_fft.c \
$(CODEC2_SRC)/interp.c \
$(CODEC2_SRC)/lsp.c \
$(CODEC2_SRC)/phase.c \
$(CODEC2_SRC)/quantise.c \
$(CODEC2_SRC)/pack.c \
$(CODEC2_SRC)/codebook.c \
$(CODEC2_SRC)/codebookd.c \
$(CODEC2_SRC)/codebookjvm.c \
$(CODEC2_SRC)/codebookge.c \
$(CODEC2_SRC)/dump.c \
$(CODEC2_SRC)/fdmdv.c \
$(CODEC2_SRC)/freedv_api.c \
$(CODEC2_SRC)/varicode.c \
$(CODEC2_SRC)/golay23.c

CFLAGS += -D__EMBEDDED__ 

#enable this for dump files to help verify optimisation
#CFLAGS += -DDUMP

CFLAGS += -I../src
CFLAGS += -I../unittest
CFLAGS += -Iinc

FFT_TEST_SRCS = \
$(DSPLIB)/Examples/arm_fft_bin_example/arm_fft_bin_data.c \
fft_test.c \
src/startup_stm32f4xx.s \
stm32f4_timer.c \
gdb_stdio.c \
../src/kiss_fft.c

###################################################

vpath %.c src
vpath %.a lib

ROOT=$(shell pwd)

# Library paths

LIBPATHS = 

# Libraries to link

LIBS = libstm32f4.a -lg -lnosys -lm

# startup file

SRCS += src/startup_stm32f4xx.s src/init.c

OBJS = $(SRCS:.c=.o)

###################################################

all: libstm32f4.a codec2_profile.elf fft_test.elf dac_ut.elf dac_play.elf adc_rec.elf pwm_ut.elf fdmdv_profile.elf sm1000_leds_switches_ut.elf sm1000.elf adcdac_ut.elf

dl/$(PERIPHLIBZIP):
	mkdir -p dl
	cd dl; wget $(PERIPHLIBURL)/$(PERIPHLIBZIP)

$(PERIPHLIBDIR): dl/$(PERIPHLIBZIP)
	cd dl; unzip $(PERIPHLIBZIP)
	mv dl/$(PERIPHLIBDIR) $(PERIPHLIBDIR)

libstm32f4.a: $(PERIPHLIBDIR)
	$(MAKE) $(STM32F4TEMPLATE)/system_stm32f4xx.o
	for F in $(STM32F4LIB)/src/*.c ; do $(MAKE) $${F%.c}.o ; done
	for F in $(CMSIS)/DSP_Lib/Source/*/*.c ; do $(MAKE) $${F%.c}.o ; done
	find $(PERIPHLIBDIR) -type f -name '*.o' -exec $(AR) crs libstm32f4.a {} ";"	

####################################################

CODEC2_PROFILE_SRCS=\
src/codec2_profile.c \
src/gdb_stdio.c \
src/stm32f4_timer.c \
src/startup_stm32f4xx.s \
src/init.c \
src/system_stm32f4xx.c
CODEC2_PROFILE_SRCS += $(CODEC2_SRCS)

codec2_profile.elf: $(CODEC2_PROFILE_SRCS) 
	$(CC) $(CFLAGS) -DTIMER $^ -o $@ $(LIBPATHS) $(LIBS)

fft_test.elf: $(FFT_TEST_SRCS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBPATHS) $(LIBS)

DAC_UT_SRCS=\
src/dac_ut.c \
../src/fifo.c \
src/stm32f4_dac.c \
src/debugblinky.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c

dac_ut.elf: $(DAC_UT_SRCS)
	$(CC) $(CFLAGS) -O0 $^ -o $@ $(LIBPATHS) $(LIBS)

ADCDAC_UT_SRCS=\
src/adcdac_ut.c \
../src/fifo.c \
src/stm32f4_dac.c \
src/stm32f4_adc.c \
src/debugblinky.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c

adcdac_ut.elf: $(ADCDAC_UT_SRCS)
	$(CC) $(CFLAGS) -O0 $^ -o $@ $(LIBPATHS) $(LIBS)

DAC_PLAY_SRCS=\
src/dac_play.c \
../src/fifo.c \
gdb_stdio.c \
src/stm32f4_dac.c \
src/debugblinky.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c

dac_play.elf: $(DAC_PLAY_SRCS)
	$(CC) $(CFLAGS) -O0 $^ -o $@ $(LIBPATHS) $(LIBS)

ADC_REC_SRCS=\
src/adc_rec.c \
../src/fifo.c \
gdb_stdio.c \
src/stm32f4_adc.c \
src/debugblinky.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c

adc_rec.elf: $(ADC_REC_SRCS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBPATHS) $(LIBS)

PWM_UT_SRCS=\
gdb_stdio.c \
src/stm32f4_pwm.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c

pwm_ut.elf: $(PWM_UT_SRCS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBPATHS) $(LIBS)

POWER_UT_SRCS=\
src/power_ut.c \
gdb_stdio.c \
../src/fifo.c \
src/stm32f4_adc.c \
src/stm32f4_dac.c \
src/debugblinky.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c \
src/stm32f4_timer.c \

POWER_UT_SRCS += $(CODEC2_SRCS)

power_ut.elf: $(POWER_UT_SRCS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBPATHS) $(LIBS)

FDMDV_PROFILE_SRCS=\
src/fdmdv_profile.c \
gdb_stdio.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c \
src/stm32f4_timer.c

FDMDV_PROFILE_SRCS += $(CODEC2_SRCS)

fdmdv_profile.elf: $(FDMDV_PROFILE_SRCS)
	$(CC) $(CFLAGS) -DTIMER $^ -o $@ $(LIBPATHS) $(LIBS)

SM1000_LEDS_SWITCHES_UT_SRCS=\
src/sm1000_leds_switches_ut.c \
src/sm1000_leds_switches.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c

sm1000_leds_switches_ut.elf: $(SM1000_LEDS_SWITCHES_UT_SRCS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBPATHS) $(LIBS)

SM1000_SRCS=\
src/sm1000_main.c \
src/sm1000_leds_switches.c \
../src/fifo.c \
src/debugblinky.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c 

SM1000_SRCS += $(CODEC2_SRCS)

src/stm32f4_dac.o: src/stm32f4_dac.c
	$(CC) $(CFLAGS)  $^ -c -o $@ 

src/stm32f4_adc.o: src/stm32f4_adc.c
	$(CC) $(CFLAGS)  $^ -c -o $@ 

sm1000.elf: $(SM1000_SRCS) src/stm32f4_dac.o src/stm32f4_adc.o
	$(CC) $(CFLAGS) -O3 $^ -o $@ $(LIBPATHS) $(LIBS)

clean:
	rm -f *.o
	rm -f *.elf
	rm -f libstm32f4.a	
	find $(PERIPHLIBDIR) -type f -name '*.o' -exec rm {} \;
