
# Include local definitions if they exist.
-include local.mak

all: tst_codec2_fft_init.elf arm_fft_bin_example_f32.elf fft_bin_example_512.elf

###################################################

CROSS_COMPILE ?= arm-none-eabi-
CC=$(BINPATH)$(CROSS_COMPILE)gcc
AS=$(BINPATH)$(CROSS_COMPILE)as
LD=$(BINPATH)$(CROSS_COMPILE)ld
OBJCOPY=$(BINPATH)$(CROSS_COMPILE)objcopy
SIZE=$(BINPATH)$(CROSS_COMPILE)size
SUDO ?= sudo

###################################################

CFLAGS  = -std=gnu11 -O2 -g -Wall -DSTM32F40_41xxx -DCORTEX_M4
CFLAGS += -mlittle-endian -mthumb -mthumb-interwork -nostartfiles -mcpu=cortex-m4 -Wno-unused-function

CFLAGS += -fsingle-precision-constant -Wdouble-promotion
CFLAGS += -fdata-sections -ffunction-sections -Xlinker --gc-sections
CFLAGS += -mfpu=fpv4-sp-d16 -mfloat-abi=hard -D__FPU_PRESENT=1 -D__FPU_USED=1


PERIPHLIBDIR	?= STM32F4xx_DSP_StdPeriph_Lib
CMSIS		= $(PERIPHLIBDIR)/Libraries/CMSIS
STM32F4LIB	= $(PERIPHLIBDIR)/Libraries/STM32F4xx_StdPeriph_Driver
STM32F4TEMPLATE	= $(PERIPHLIBDIR)/Project/STM32F4xx_StdPeriph_Templates
DSPLIB          = $(PERIPHLIBDIR)/Libraries/CMSIS/DSP_Lib

CFLAGS		+= -DUSE_STDPERIPH_DRIVER -I$(STM32F4LIB)/inc -I$(STM32F4TEMPLATE)
CFLAGS		+= -I$(CMSIS)/Include -I$(CMSIS)/Device/ST/STM32F4xx/Include
CFLAGS		+= -DARM_MATH_CM4
CFLAGS		+= -DFDV_ARM_MATH
CFLAGS		+= -DSEMIHOST_USE_STDIO


# Precious files that should be preserved at all cost!
.PRECIOUS: dl/$(PERIPHLIBZIP)

STM32F4LIB_SRCS=\
$(STM32F4LIB)/src/misc.c\
$(STM32F4LIB)/src/stm32f4xx_adc.c\
$(STM32F4LIB)/src/stm32f4xx_can.c\
$(STM32F4LIB)/src/stm32f4xx_cec.c\
$(STM32F4LIB)/src/stm32f4xx_crc.c\
$(STM32F4LIB)/src/stm32f4xx_cryp_aes.c\
$(STM32F4LIB)/src/stm32f4xx_cryp.c\
$(STM32F4LIB)/src/stm32f4xx_cryp_des.c\
$(STM32F4LIB)/src/stm32f4xx_cryp_tdes.c\
$(STM32F4LIB)/src/stm32f4xx_dac.c\
$(STM32F4LIB)/src/stm32f4xx_dbgmcu.c\
$(STM32F4LIB)/src/stm32f4xx_dcmi.c\
$(STM32F4LIB)/src/stm32f4xx_dma2d.c\
$(STM32F4LIB)/src/stm32f4xx_dma.c\
$(STM32F4LIB)/src/stm32f4xx_exti.c\
$(STM32F4LIB)/src/stm32f4xx_flash.c\
$(STM32F4LIB)/src/stm32f4xx_flash_ramfunc.c\
$(STM32F4LIB)/src/stm32f4xx_fmpi2c.c\
$(STM32F4LIB)/src/stm32f4xx_fsmc.c\
$(STM32F4LIB)/src/stm32f4xx_gpio.c\
$(STM32F4LIB)/src/stm32f4xx_hash.c\
$(STM32F4LIB)/src/stm32f4xx_hash_md5.c\
$(STM32F4LIB)/src/stm32f4xx_hash_sha1.c\
$(STM32F4LIB)/src/stm32f4xx_i2c.c\
$(STM32F4LIB)/src/stm32f4xx_iwdg.c\
$(STM32F4LIB)/src/stm32f4xx_ltdc.c\
$(STM32F4LIB)/src/stm32f4xx_pwr.c\
$(STM32F4LIB)/src/stm32f4xx_qspi.c\
$(STM32F4LIB)/src/stm32f4xx_rcc.c\
$(STM32F4LIB)/src/stm32f4xx_rng.c\
$(STM32F4LIB)/src/stm32f4xx_rtc.c\
$(STM32F4LIB)/src/stm32f4xx_sai.c\
$(STM32F4LIB)/src/stm32f4xx_sdio.c\
$(STM32F4LIB)/src/stm32f4xx_spdifrx.c\
$(STM32F4LIB)/src/stm32f4xx_spi.c\
$(STM32F4LIB)/src/stm32f4xx_syscfg.c\
$(STM32F4LIB)/src/stm32f4xx_tim.c\
$(STM32F4LIB)/src/stm32f4xx_usart.c\
$(STM32F4LIB)/src/stm32f4xx_wwdg.c

STM32F4LIB_OBJS = $(STM32F4LIB_SRCS:.c=.o)


CMSIS_SRCS=\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_abs_f32.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_add_f32.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_dot_prod_f32.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_mult_f32.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_negate_f32.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_offset_f32.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_scale_f32.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_sub_f32.c\
$(CMSIS)/DSP_Lib/Source/CommonTables/arm_common_tables.c\
$(CMSIS)/DSP_Lib/Source/CommonTables/arm_const_structs.c\
$(CMSIS)/DSP_Lib/Source/CommonTables/arm_const_structs.c\
$(CMSIS)/DSP_Lib/Source/ComplexMathFunctions/arm_cmplx_mag_f32.c\
$(CMSIS)/DSP_Lib/Source/SupportFunctions/arm_copy_f32.c\
$(CMSIS)/DSP_Lib/Source/SupportFunctions/arm_fill_f32.c\
$(CMSIS)/DSP_Lib/Source/StatisticsFunctions/arm_max_f32.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_bitreversal.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_cfft_f32.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_cfft_radix2_f32.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_cfft_radix2_init_f32.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_cfft_radix4_f32.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_cfft_radix4_init_f32.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_cfft_radix8_f32.c\

CMSIS_OBJS = $(CMSIS_SRCS:.c=.o) $(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_bitreversal2.o

libstm32f4.a: $(CMSIS_OBJS) $(STM32F4LIB_OBJS)
	find -L $(PERIPHLIBDIR) -type f -name '*.o' -exec $(AR) crs libstm32f4.a {} ";"

LDFLAGS += -T stm32_flash.ld 

# Standard ARM semihosting
LIBS = -lg -lrdimon -lm --specs=rdimon.specs


OBJS =  semihosting.o \
	startup_stm32f4xx.o init.o stm32f4_machdep.o system_stm32f4xx.o

#################
tst_codec2_fft_init.elf: tst_codec2_fft_init.o  \
	$(OBJS) libstm32f4.a
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LIBPATHS) $(LIBS)
#
tst_codec2_fft_init.o: tst_codec2_fft_init.c

#################
arm_fft_bin_example_f32.elf: arm_fft_bin_example_f32.o arm_fft_bin_data.o \
	$(OBJS) libstm32f4.a
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LIBPATHS) $(LIBS)
#
arm_fft_bin_example_f32.o: arm_fft_bin_example_f32.c 
#
arm_fft_bin_data.o: arm_fft_bin_data.c

#################
fft_bin_example_512.elf: fft_bin_example_512.o arm_fft_bin_data.o \
	$(OBJS) libstm32f4.a
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LIBPATHS) $(LIBS)
#
fft_bin_example_512.o: fft_bin_example_512.c 
