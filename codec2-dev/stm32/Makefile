# Makefile for stm32f4 Codec 2 test programs

# Include local definitions if they exist.
-include local.mak

###################################################

FLOAT_TYPE=hard

###################################################

CROSS_COMPILE ?= arm-none-eabi-
CC=$(BINPATH)$(CROSS_COMPILE)gcc
AS=$(BINPATH)$(CROSS_COMPILE)as
OBJCOPY=$(BINPATH)$(CROSS_COMPILE)objcopy
SIZE=$(BINPATH)$(CROSS_COMPILE)size
SUDO ?= sudo

###################################################

CFLAGS  = -std=gnu11 -O0 -g -Wall -Tstm32_flash.ld -DSTM32F40_41xxx -DCORTEX_M4
CFLAGS += -mlittle-endian -mthumb -mthumb-interwork -nostartfiles -mcpu=cortex-m4 -Wno-unused-function

ifeq ($(FLOAT_TYPE), hard)
CFLAGS += -fsingle-precision-constant -Wdouble-promotion
CFLAGS += -fdata-sections -ffunction-sections -Xlinker --gc-sections
CFLAGS += -mfpu=fpv4-sp-d16 -mfloat-abi=hard -D__FPU_PRESENT=1 -D__FPU_USED=1
else
CFLAGS += -msoft-float
endif

###################################################

# Definitions for the STM32F4 Standard Peripheral Library

PERIPHLIBURL    = http://www.st.com/resource/en/firmware
PERIPHLIBZIP    = stm32f4_dsp_stdperiph_lib.zip
PERIPHLIBVER	= V1.8.0
PERIPHLIBNAME	= STM32F4xx_DSP_StdPeriph_Lib
PERIPHLIBDIR	= $(PERIPHLIBNAME)
CMSIS		= $(PERIPHLIBDIR)/Libraries/CMSIS
STM32F4LIB	= $(PERIPHLIBDIR)/Libraries/STM32F4xx_StdPeriph_Driver
STM32F4TEMPLATE	= $(PERIPHLIBDIR)/Project/STM32F4xx_StdPeriph_Templates
DSPLIB          = $(PERIPHLIBDIR)/Libraries/CMSIS/DSP_Lib

CFLAGS		+= -DUSE_STDPERIPH_DRIVER -I$(STM32F4LIB)/inc -I$(STM32F4TEMPLATE)
CFLAGS		+= -I$(CMSIS)/Include -I$(CMSIS)/Device/ST/STM32F4xx/Include
CFLAGS		+= -DARM_MATH_CM4

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

# Not compiling for now
# $(STM32F4LIB)/src/stm32f4xx_fmc.c

STM32F4LIB_OBJS = $(STM32F4LIB_SRCS:.c=.o)

CMSIS_SRCS=\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_abs_f32.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_abs_q15.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_abs_q31.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_abs_q7.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_add_f32.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_add_q15.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_add_q31.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_add_q7.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_dot_prod_f32.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_dot_prod_q15.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_dot_prod_q31.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_dot_prod_q7.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_mult_f32.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_mult_q15.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_mult_q31.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_mult_q7.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_negate_f32.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_negate_q15.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_negate_q31.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_negate_q7.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_offset_f32.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_offset_q15.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_offset_q31.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_offset_q7.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_scale_f32.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_scale_q15.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_scale_q31.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_scale_q7.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_shift_q15.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_shift_q31.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_shift_q7.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_sub_f32.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_sub_q15.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_sub_q31.c\
$(CMSIS)/DSP_Lib/Source/BasicMathFunctions/arm_sub_q7.c\
$(CMSIS)/DSP_Lib/Source/CommonTables/arm_common_tables.c\
$(CMSIS)/DSP_Lib/Source/CommonTables/arm_const_structs.c\
$(CMSIS)/DSP_Lib/Source/ComplexMathFunctions/arm_cmplx_conj_f32.c\
$(CMSIS)/DSP_Lib/Source/ComplexMathFunctions/arm_cmplx_conj_q15.c\
$(CMSIS)/DSP_Lib/Source/ComplexMathFunctions/arm_cmplx_conj_q31.c\
$(CMSIS)/DSP_Lib/Source/ComplexMathFunctions/arm_cmplx_dot_prod_f32.c\
$(CMSIS)/DSP_Lib/Source/ComplexMathFunctions/arm_cmplx_dot_prod_q15.c\
$(CMSIS)/DSP_Lib/Source/ComplexMathFunctions/arm_cmplx_dot_prod_q31.c\
$(CMSIS)/DSP_Lib/Source/ComplexMathFunctions/arm_cmplx_mag_f32.c\
$(CMSIS)/DSP_Lib/Source/ComplexMathFunctions/arm_cmplx_mag_q15.c\
$(CMSIS)/DSP_Lib/Source/ComplexMathFunctions/arm_cmplx_mag_q31.c\
$(CMSIS)/DSP_Lib/Source/ComplexMathFunctions/arm_cmplx_mag_squared_f32.c\
$(CMSIS)/DSP_Lib/Source/ComplexMathFunctions/arm_cmplx_mag_squared_q15.c\
$(CMSIS)/DSP_Lib/Source/ComplexMathFunctions/arm_cmplx_mag_squared_q31.c\
$(CMSIS)/DSP_Lib/Source/ComplexMathFunctions/arm_cmplx_mult_cmplx_f32.c\
$(CMSIS)/DSP_Lib/Source/ComplexMathFunctions/arm_cmplx_mult_cmplx_q15.c\
$(CMSIS)/DSP_Lib/Source/ComplexMathFunctions/arm_cmplx_mult_cmplx_q31.c\
$(CMSIS)/DSP_Lib/Source/ComplexMathFunctions/arm_cmplx_mult_real_f32.c\
$(CMSIS)/DSP_Lib/Source/ComplexMathFunctions/arm_cmplx_mult_real_q15.c\
$(CMSIS)/DSP_Lib/Source/ComplexMathFunctions/arm_cmplx_mult_real_q31.c\
$(CMSIS)/DSP_Lib/Source/ControllerFunctions/arm_pid_init_f32.c\
$(CMSIS)/DSP_Lib/Source/ControllerFunctions/arm_pid_init_q15.c\
$(CMSIS)/DSP_Lib/Source/ControllerFunctions/arm_pid_init_q31.c\
$(CMSIS)/DSP_Lib/Source/ControllerFunctions/arm_pid_reset_f32.c\
$(CMSIS)/DSP_Lib/Source/ControllerFunctions/arm_pid_reset_q15.c\
$(CMSIS)/DSP_Lib/Source/ControllerFunctions/arm_pid_reset_q31.c\
$(CMSIS)/DSP_Lib/Source/ControllerFunctions/arm_sin_cos_f32.c\
$(CMSIS)/DSP_Lib/Source/ControllerFunctions/arm_sin_cos_q31.c\
$(CMSIS)/DSP_Lib/Source/FastMathFunctions/arm_cos_f32.c\
$(CMSIS)/DSP_Lib/Source/FastMathFunctions/arm_cos_q15.c\
$(CMSIS)/DSP_Lib/Source/FastMathFunctions/arm_cos_q31.c\
$(CMSIS)/DSP_Lib/Source/FastMathFunctions/arm_sin_f32.c\
$(CMSIS)/DSP_Lib/Source/FastMathFunctions/arm_sin_q15.c\
$(CMSIS)/DSP_Lib/Source/FastMathFunctions/arm_sin_q31.c\
$(CMSIS)/DSP_Lib/Source/FastMathFunctions/arm_sqrt_q15.c\
$(CMSIS)/DSP_Lib/Source/FastMathFunctions/arm_sqrt_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_biquad_cascade_df1_32x64_init_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_biquad_cascade_df1_32x64_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_biquad_cascade_df1_f32.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_biquad_cascade_df1_fast_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_biquad_cascade_df1_fast_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_biquad_cascade_df1_init_f32.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_biquad_cascade_df1_init_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_biquad_cascade_df1_init_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_biquad_cascade_df1_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_biquad_cascade_df1_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_biquad_cascade_df2T_f32.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_biquad_cascade_df2T_f64.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_biquad_cascade_df2T_init_f32.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_biquad_cascade_df2T_init_f64.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_biquad_cascade_stereo_df2T_f32.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_biquad_cascade_stereo_df2T_init_f32.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_conv_f32.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_conv_fast_opt_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_conv_fast_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_conv_fast_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_conv_opt_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_conv_opt_q7.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_conv_partial_f32.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_conv_partial_fast_opt_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_conv_partial_fast_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_conv_partial_fast_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_conv_partial_opt_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_conv_partial_opt_q7.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_conv_partial_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_conv_partial_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_conv_partial_q7.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_conv_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_conv_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_conv_q7.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_correlate_f32.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_correlate_fast_opt_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_correlate_fast_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_correlate_fast_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_correlate_opt_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_correlate_opt_q7.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_correlate_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_correlate_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_correlate_q7.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_decimate_f32.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_decimate_fast_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_decimate_fast_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_decimate_init_f32.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_decimate_init_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_decimate_init_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_decimate_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_decimate_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_f32.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_fast_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_fast_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_init_f32.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_init_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_init_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_init_q7.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_interpolate_f32.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_interpolate_init_f32.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_interpolate_init_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_interpolate_init_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_interpolate_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_interpolate_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_lattice_f32.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_lattice_init_f32.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_lattice_init_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_lattice_init_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_lattice_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_lattice_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_q7.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_sparse_f32.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_sparse_init_f32.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_sparse_init_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_sparse_init_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_sparse_init_q7.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_sparse_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_sparse_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_fir_sparse_q7.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_iir_lattice_f32.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_iir_lattice_init_f32.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_iir_lattice_init_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_iir_lattice_init_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_iir_lattice_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_iir_lattice_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_lms_f32.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_lms_init_f32.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_lms_init_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_lms_init_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_lms_norm_f32.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_lms_norm_init_f32.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_lms_norm_init_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_lms_norm_init_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_lms_norm_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_lms_norm_q31.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_lms_q15.c\
$(CMSIS)/DSP_Lib/Source/FilteringFunctions/arm_lms_q31.c\
$(CMSIS)/DSP_Lib/Source/MatrixFunctions/arm_mat_add_f32.c\
$(CMSIS)/DSP_Lib/Source/MatrixFunctions/arm_mat_add_q15.c\
$(CMSIS)/DSP_Lib/Source/MatrixFunctions/arm_mat_add_q31.c\
$(CMSIS)/DSP_Lib/Source/MatrixFunctions/arm_mat_cmplx_mult_f32.c\
$(CMSIS)/DSP_Lib/Source/MatrixFunctions/arm_mat_cmplx_mult_q15.c\
$(CMSIS)/DSP_Lib/Source/MatrixFunctions/arm_mat_cmplx_mult_q31.c\
$(CMSIS)/DSP_Lib/Source/MatrixFunctions/arm_mat_init_f32.c\
$(CMSIS)/DSP_Lib/Source/MatrixFunctions/arm_mat_init_q15.c\
$(CMSIS)/DSP_Lib/Source/MatrixFunctions/arm_mat_init_q31.c\
$(CMSIS)/DSP_Lib/Source/MatrixFunctions/arm_mat_inverse_f32.c\
$(CMSIS)/DSP_Lib/Source/MatrixFunctions/arm_mat_inverse_f64.c\
$(CMSIS)/DSP_Lib/Source/MatrixFunctions/arm_mat_mult_f32.c\
$(CMSIS)/DSP_Lib/Source/MatrixFunctions/arm_mat_mult_fast_q15.c\
$(CMSIS)/DSP_Lib/Source/MatrixFunctions/arm_mat_mult_fast_q31.c\
$(CMSIS)/DSP_Lib/Source/MatrixFunctions/arm_mat_mult_q15.c\
$(CMSIS)/DSP_Lib/Source/MatrixFunctions/arm_mat_mult_q31.c\
$(CMSIS)/DSP_Lib/Source/MatrixFunctions/arm_mat_scale_f32.c\
$(CMSIS)/DSP_Lib/Source/MatrixFunctions/arm_mat_scale_q15.c\
$(CMSIS)/DSP_Lib/Source/MatrixFunctions/arm_mat_scale_q31.c\
$(CMSIS)/DSP_Lib/Source/MatrixFunctions/arm_mat_sub_f32.c\
$(CMSIS)/DSP_Lib/Source/MatrixFunctions/arm_mat_sub_q15.c\
$(CMSIS)/DSP_Lib/Source/MatrixFunctions/arm_mat_sub_q31.c\
$(CMSIS)/DSP_Lib/Source/MatrixFunctions/arm_mat_trans_f32.c\
$(CMSIS)/DSP_Lib/Source/MatrixFunctions/arm_mat_trans_q15.c\
$(CMSIS)/DSP_Lib/Source/MatrixFunctions/arm_mat_trans_q31.c\
$(CMSIS)/DSP_Lib/Source/StatisticsFunctions/arm_max_f32.c\
$(CMSIS)/DSP_Lib/Source/StatisticsFunctions/arm_max_q15.c\
$(CMSIS)/DSP_Lib/Source/StatisticsFunctions/arm_max_q31.c\
$(CMSIS)/DSP_Lib/Source/StatisticsFunctions/arm_max_q7.c\
$(CMSIS)/DSP_Lib/Source/StatisticsFunctions/arm_mean_f32.c\
$(CMSIS)/DSP_Lib/Source/StatisticsFunctions/arm_mean_q15.c\
$(CMSIS)/DSP_Lib/Source/StatisticsFunctions/arm_mean_q31.c\
$(CMSIS)/DSP_Lib/Source/StatisticsFunctions/arm_mean_q7.c\
$(CMSIS)/DSP_Lib/Source/StatisticsFunctions/arm_min_f32.c\
$(CMSIS)/DSP_Lib/Source/StatisticsFunctions/arm_min_q15.c\
$(CMSIS)/DSP_Lib/Source/StatisticsFunctions/arm_min_q31.c\
$(CMSIS)/DSP_Lib/Source/StatisticsFunctions/arm_min_q7.c\
$(CMSIS)/DSP_Lib/Source/StatisticsFunctions/arm_power_f32.c\
$(CMSIS)/DSP_Lib/Source/StatisticsFunctions/arm_power_q15.c\
$(CMSIS)/DSP_Lib/Source/StatisticsFunctions/arm_power_q31.c\
$(CMSIS)/DSP_Lib/Source/StatisticsFunctions/arm_power_q7.c\
$(CMSIS)/DSP_Lib/Source/StatisticsFunctions/arm_rms_f32.c\
$(CMSIS)/DSP_Lib/Source/StatisticsFunctions/arm_rms_q15.c\
$(CMSIS)/DSP_Lib/Source/StatisticsFunctions/arm_rms_q31.c\
$(CMSIS)/DSP_Lib/Source/StatisticsFunctions/arm_std_f32.c\
$(CMSIS)/DSP_Lib/Source/StatisticsFunctions/arm_std_q15.c\
$(CMSIS)/DSP_Lib/Source/StatisticsFunctions/arm_std_q31.c\
$(CMSIS)/DSP_Lib/Source/StatisticsFunctions/arm_var_f32.c\
$(CMSIS)/DSP_Lib/Source/StatisticsFunctions/arm_var_q15.c\
$(CMSIS)/DSP_Lib/Source/StatisticsFunctions/arm_var_q31.c\
$(CMSIS)/DSP_Lib/Source/SupportFunctions/arm_copy_f32.c\
$(CMSIS)/DSP_Lib/Source/SupportFunctions/arm_copy_q15.c\
$(CMSIS)/DSP_Lib/Source/SupportFunctions/arm_copy_q31.c\
$(CMSIS)/DSP_Lib/Source/SupportFunctions/arm_copy_q7.c\
$(CMSIS)/DSP_Lib/Source/SupportFunctions/arm_fill_f32.c\
$(CMSIS)/DSP_Lib/Source/SupportFunctions/arm_fill_q15.c\
$(CMSIS)/DSP_Lib/Source/SupportFunctions/arm_fill_q31.c\
$(CMSIS)/DSP_Lib/Source/SupportFunctions/arm_fill_q7.c\
$(CMSIS)/DSP_Lib/Source/SupportFunctions/arm_float_to_q15.c\
$(CMSIS)/DSP_Lib/Source/SupportFunctions/arm_float_to_q31.c\
$(CMSIS)/DSP_Lib/Source/SupportFunctions/arm_float_to_q7.c\
$(CMSIS)/DSP_Lib/Source/SupportFunctions/arm_q15_to_float.c\
$(CMSIS)/DSP_Lib/Source/SupportFunctions/arm_q15_to_q31.c\
$(CMSIS)/DSP_Lib/Source/SupportFunctions/arm_q15_to_q7.c\
$(CMSIS)/DSP_Lib/Source/SupportFunctions/arm_q31_to_float.c\
$(CMSIS)/DSP_Lib/Source/SupportFunctions/arm_q31_to_q15.c\
$(CMSIS)/DSP_Lib/Source/SupportFunctions/arm_q31_to_q7.c\
$(CMSIS)/DSP_Lib/Source/SupportFunctions/arm_q7_to_float.c\
$(CMSIS)/DSP_Lib/Source/SupportFunctions/arm_q7_to_q15.c\
$(CMSIS)/DSP_Lib/Source/SupportFunctions/arm_q7_to_q31.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_bitreversal.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_cfft_f32.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_cfft_q15.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_cfft_q31.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_cfft_radix2_f32.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_cfft_radix2_init_f32.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_cfft_radix2_init_q15.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_cfft_radix2_init_q31.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_cfft_radix2_q15.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_cfft_radix2_q31.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_cfft_radix4_f32.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_cfft_radix4_init_f32.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_cfft_radix4_init_q15.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_cfft_radix4_init_q31.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_cfft_radix4_q15.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_cfft_radix4_q31.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_cfft_radix8_f32.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_dct4_f32.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_dct4_init_f32.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_dct4_init_q15.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_dct4_init_q31.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_dct4_q15.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_dct4_q31.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_rfft_f32.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_rfft_fast_f32.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_rfft_fast_init_f32.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_rfft_init_f32.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_rfft_init_q15.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_rfft_init_q31.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_rfft_q15.c\
$(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_rfft_q31.c\


CMSIS_OBJS = $(CMSIS_SRCS:.c=.o) $(CMSIS)/DSP_Lib/Source/TransformFunctions/arm_bitreversal2.o

###################################################

# Codec 2

CODEC2_SRC=../src
CODEC2_SRCS=\
$(CODEC2_SRC)/lpc.c \
$(CODEC2_SRC)/nlp.c \
$(CODEC2_SRC)/postfilter.c \
$(CODEC2_SRC)/sine.c \
$(CODEC2_SRC)/codec2.c \
$(CODEC2_SRC)/codec2_fft.c \
$(CODEC2_SRC)/kiss_fft.c \
$(CODEC2_SRC)/kiss_fftr.c \
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
$(CODEC2_SRC)/golay23.c \
$(CODEC2_SRC)/fsk.c \
$(CODEC2_SRC)/fmfsk.c \
$(CODEC2_SRC)/freedv_vhf_framing.c \
$(CODEC2_SRC)/freedv_data_channel.c

CFLAGS += -D__EMBEDDED__

#enable this for dump files to help verify optimisation
#CFLAGS += -DDUMP

CFLAGS += -I../src
CFLAGS += -I../unittest
CFLAGS += -Iinc

FFT_TEST_SRCS = \
$(DSPLIB)/Examples/arm_fft_bin_example/GCC/arm_fft_bin_data.c \
fft_test.c \
src/startup_stm32f4xx.s \
src/system_stm32f4xx.c \
stm32f4_machdep.c \
gdb_stdio.c \
../src/kiss_fft.c

###################################################

vpath %.c src
vpath %.a lib

ROOT=$(shell pwd)

# Library paths

LIBPATHS =

# Libraries to link

LIBS = -lg -lnosys -lm 
#Uncomment for standard arm semihosting
#LIBS = -lg -lrdimon -lm --specs=rdimon.specs

# startup file

SRCS += src/startup_stm32f4xx.s src/init.c

OBJS = $(SRCS:.c=.o)

all: libstm32f4.a codec2_profile.bin fft_test.bin dac_ut.bin dac_play.bin adc_rec.bin pwm_ut.bin fdmdv_profile.bin sm1000_leds_switches_ut.bin sm1000.bin adcdac_ut.bin freedv_tx_profile.bin freedv_rx_profile.bin adc_sd.bin usb_vcp_ut.bin tuner_ut.bin fast_dac_ut.bin adc_sfdr_ut.bin adc_rec_usb.bin si5351_ut.bin mco_ut.bin sm2000_stw.bin sm2000_adcdump.bin sm2000_rxdemo.bin

# Rule for making directories automatically.
# Note we don't use -p as it's a GNU extension.
%/.md:
	parent=$(shell dirname $(@D) ); \
		[ -d $${parent} ] || $(MAKE) $${parent}/.md
	[ -d $(@D) ] || mkdir $(@D)
	touch $@

dl/$(PERIPHLIBZIP): dl/.md
	wget -O$@.part -c $(PERIPHLIBURL)/$(PERIPHLIBZIP)
	mv $@.part $@
	touch $@

$(PERIPHLIBDIR)/.unpack: dl/$(PERIPHLIBZIP)
	test ! -d $(PERIPHLIBDIR)_$(PERIPHLIBVER) || \
		rm -fr $(PERIPHLIBDIR)_$(PERIPHLIBVER)
	unzip dl/$(PERIPHLIBZIP)
	test ! -d $(PERIPHLIBDIR) || rm -fr $(PERIPHLIBDIR)
	mv $(PERIPHLIBDIR)_$(PERIPHLIBVER) $(PERIPHLIBDIR)
	touch $@

$(CMSIS_OBJS) $(STM32F4LIB_OBJS): $(PERIPHLIBDIR)/.unpack

libstm32f4.a: $(CMSIS_OBJS) $(STM32F4LIB_OBJS)
	find $(PERIPHLIBDIR) -type f -name '*.o' -exec $(AR) crs libstm32f4.a {} ";"

# Kludgy target to build a file with CFLAGS -O3
%.O3.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -O3 -c -o $@ $<

# Kludgy target to build a file with CFLAGS -DPROFILE
%.profile.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -DPROFILE -c -o $@ $<

# Rule for building .bin files from a .elf
%.bin: %.elf
	$(OBJCOPY) -O binary $< $@

# Rule for programming the SM1000
%.pgm: %.bin
	$(SUDO) dfu-util -d 0483:df11 -c 1 -i 0 -a 0 -s 0x08000000 -D $<

####################################################

CODEC2_PROFILE_SRCS=\
src/codec2_profile.c \
src/gdb_stdio.c \
src/stm32f4_machdep.c \
src/startup_stm32f4xx.s \
src/init.c \
src/system_stm32f4xx.c
CODEC2_PROFILE_SRCS += $(CODEC2_SRCS)

codec2_profile.elf: $(CODEC2_PROFILE_SRCS:.c=.profile.o) libstm32f4.a
	$(CC) $(CFLAGS) -DPROFILE $^ -o $@ $(LIBPATHS) $(LIBS)

fft_test.elf: $(FFT_TEST_SRCS:.c=.o) libstm32f4.a
	$(CC) $(CFLAGS) $^ -o $@ $(LIBPATHS) $(LIBS)

DAC_UT_SRCS=\
src/dac_ut.c \
../src/fifo.c \
src/stm32f4_dac.c \
src/debugblinky.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c

dac_ut.elf: $(DAC_UT_SRCS:.c=.o) libstm32f4.a
	$(CC) $(CFLAGS) -O0 $^ -o $@ $(LIBPATHS) $(LIBS)

FAST_DAC_UT_SRCS=\
src/fast_dac_ut.c \
../src/fifo.c \
src/iir_duc.c \
src/gdb_stdio.c \
src/stm32f4_dacduc.c \
src/debugblinky.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c

fast_dac_ut.elf: $(FAST_DAC_UT_SRCS:.c=.O3.o) libstm32f4.a
	$(CC) $(CFLAGS) -O3 $^ -o $@ $(LIBPATHS) $(LIBS)

ADCDAC_UT_SRCS=\
src/adcdac_ut.c \
../src/fifo.c \
src/stm32f4_dac.c \
src/stm32f4_adc.c \
src/sm1000_leds_switches.c \
src/debugblinky.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c

adcdac_ut.elf: $(ADCDAC_UT_SRCS:.c=.o) libstm32f4.a
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

dac_play.elf: $(DAC_PLAY_SRCS:.c=.o) libstm32f4.a
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

adc_rec.elf: $(ADC_REC_SRCS:.c=.o) libstm32f4.a
	$(CC) $(CFLAGS) $^ -o $@ $(LIBPATHS) $(LIBS)

ADC_SD_SRCS=\
src/adc_sd.c \
../src/fifo.c \
gdb_stdio.c \
src/stm32f4_adc.c \
src/stm32f4_dac.c \
src/debugblinky.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c

adc_sd.elf: $(ADC_SD_SRCS:.c=.o) libstm32f4.a
	$(CC) $(CFLAGS) $^ -o $@ $(LIBPATHS) $(LIBS)

PWM_UT_SRCS=\
gdb_stdio.c \
src/stm32f4_pwm.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c

pwm_ut.elf: $(PWM_UT_SRCS:.c=.o) libstm32f4.a
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
src/stm32f4_machdep.c \

POWER_UT_SRCS += $(CODEC2_SRCS)

power_ut.elf: $(POWER_UT_SRCS:.c=.o) libstm32f4.a
	$(CC) $(CFLAGS) $^ -o $@ $(LIBPATHS) $(LIBS)

USB_VCP=\
usb_conf/usb_bsp.c \
usb_conf/usbd_desc.c \
usb_conf/usbd_usr.c \
usb_lib/cdc/usbd_cdc_core.c \
usb_lib/cdc/usbd_cdc_vcp.c \
usb_lib/core/usbd_core.c \
usb_lib/core/usbd_ioreq.c \
usb_lib/core/usbd_req.c \
usb_lib/otg/usb_core.c \
usb_lib/otg/usb_dcd.c \
usb_lib/otg/usb_dcd_int.c

USB_VCP_UT=\
src/usb_vcp_ut.c \
src/stm32f4_usb_vcp.c \
src/sm1000_leds_switches.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c

USB_VCP_UT+=$(USB_VCP)

CFLAGS += -DUSE_USB_OTG_FS -DUSE_ULPI_PHY -Iusb_conf -Iusb_lib/cdc -Iusb_lib/core -Iusb_lib/otg

usb_vcp_ut.elf: $(USB_VCP_UT:.c=.o) libstm32f4.a
	$(CC) $(CFLAGS) $^ -o $@ $(LIBPATHS) $(LIBS)

ADC_REC_USB_SRCS=\
src/adc_rec_usb.c \
../src/fifo.c \
src/stm32f4_adc.c \
src/stm32f4_usb_vcp.c \
src/sm1000_leds_switches.c \
src/debugblinky.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c

ADC_REC_USB_SRCS+=$(USB_VCP)

adc_rec_usb.elf: $(ADC_REC_USB_SRCS:.c=.o) libstm32f4.a
	$(CC) $(CFLAGS) $^ -o $@ $(LIBPATHS) $(LIBS)


FDMDV_PROFILE_SRCS=\
src/fdmdv_profile.c \
gdb_stdio.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c \
src/stm32f4_machdep.c

FDMDV_PROFILE_SRCS += $(CODEC2_SRCS)

fdmdv_profile.elf: $(FDMDV_PROFILE_SRCS:.c=.profile.o) libstm32f4.a
	$(CC) $(CFLAGS) -DPROFILE $^ -o $@ $(LIBPATHS) $(LIBS)

SM1000_LEDS_SWITCHES_UT_SRCS=\
src/sm1000_leds_switches_ut.c \
src/sm1000_leds_switches.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c

sm1000_leds_switches_ut.elf: $(SM1000_LEDS_SWITCHES_UT_SRCS:.c=.o) \
		libstm32f4.a
	$(CC) $(CFLAGS) $^ -o $@ $(LIBPATHS) $(LIBS)

SM1000_SRCS=\
src/sm1000_main.c \
src/tone.c \
src/sfx.c \
src/sounds.c \
src/morse.c \
src/menu.c \
src/tot.c \
src/sm1000_leds_switches.c \
../src/fifo.c \
src/debugblinky.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/stm32f4_vrom.c \
src/init.c

SM1000_SRCS += $(CODEC2_SRCS)

src/stm32f4_dac.o: src/stm32f4_dac.c
	$(CC) $(CFLAGS)  $^ -c -o $@

src/stm32f4_adc.o: src/stm32f4_adc.c
	$(CC) $(CFLAGS)  $^ -c -o $@

sm1000.elf: $(SM1000_SRCS:.c=.O3.o) src/stm32f4_dac.O3.o \
		src/stm32f4_adc.O3.o libstm32f4.a
	$(CC) $(CFLAGS) -O3 $^ -o $@ $(LIBPATHS) $(LIBS)

FREEDV_TX_PROFILE_SRCS=\
src/freedv_tx_profile.c \
src/stm32f4_machdep.c \
gdb_stdio.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c

FREEDV_TX_PROFILE_SRCS += $(CODEC2_SRCS)

freedv_tx_profile.elf: $(FREEDV_TX_PROFILE_SRCS:.c=.profile.o) libstm32f4.a
	$(CC) $(CFLAGS) -DPROFILE $^ -o $@ $(LIBPATHS) $(LIBS)

FREEDV_RX_PROFILE_SRCS=\
src/freedv_rx_profile.c \
src/stm32f4_machdep.c \
gdb_stdio.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c

FREEDV_RX_PROFILE_SRCS += $(CODEC2_SRCS)

freedv_rx_profile.elf: $(FREEDV_RX_PROFILE_SRCS:.c=.o) libstm32f4.a
	$(CC) $(CFLAGS) $^ -o $@ $(LIBPATHS) $(LIBS)

FDMDV_DUMP_RT_SRCS=\
src/fdmdv_dump_rt.c \
src/sm1000_leds_switches.c \
../src/fifo.c \
src/debugblinky.c \
gdb_stdio.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c

FDMDV_DUMP_RT_SRCS += $(CODEC2_SRCS)

fdmdv_dump_rt.elf: $(FDMDV_DUMP_RT_SRCS:.c=.O3.o) \
		src/stm32f4_dac.O3.o src/stm32f4_adc.o libstm32f4.a
	$(CC) $(CFLAGS) -O3 $^ -o $@ $(LIBPATHS) $(LIBS)

# ---------------------------------------------------------------------------

TUNER_UT_SRCS=\
src/tuner_ut.c \
gdb_stdio.c \
../src/fifo.c \
src/stm32f4_dac.c \
src/iir_tuner.c \
src/sm1000_leds_switches.c \
src/debugblinky.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c \
../src/fm.c

# this needs to be compiled without the optimiser or ugly things happen
# would be nice to work out why as ISRs need to run fast

src/stm32f4_adc_tuner.o: src/stm32f4_adc_tuner.c
	$(CC) $(CFLAGS) $^ -c -o $@

tuner_ut.elf: $(TUNER_UT_SRCS:.c=.O3.o) \
		src/stm32f4_adc_tuner.o libstm32f4.a
	$(CC) $(CFLAGS) -O3 $^ -o $@ $(LIBPATHS) $(LIBS)

# ---------------------------------------------------------------------------

ADC_SFDR_UT_SRCS=\
src/adc_sfdr_ut.c \
gdb_stdio.c \
../src/fifo.c \
src/iir_tuner.c \
src/sm1000_leds_switches.c \
src/debugblinky.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c \

adc_sfdr_ut.elf: $(ADC_SFDR_UT_SRCS:.c=.O3.o) src/stm32f4_adc_tuner.o \
		libstm32f4.a
	$(CC) $(CFLAGS) -O3 $^ -o $@ $(LIBPATHS) $(LIBS)


FM_LODUC_PLAY_SRCS=\
src/fm_loduc_play.c \
gdb_stdio.c \
../src/fifo.c \
../src/fm.c \
src/debugblinky.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c

src/stm32f4_dacloduc.o: src/stm32f4_dacloduc.c
	$(CC) $(CFLAGS)  $^ -c -o $@ 

fm_loduc_play.elf: $(FM_LODUC_PLAY_SRCS) src/stm32f4_dacloduc.o
	$(CC) $(CFLAGS) -O3 $^ -o $@ $(LIBPATHS) $(LIBS)

# ---------------------------------------------------------------------------------

SI5351_UT_SRCS=\
src/si5351_ut.c \
src/new_i2c.c \
src/si53xx.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c \

si5351_ut.elf: $(SI5351_UT_SRCS:.c=.o) libstm32f4.a
	$(CC) $(CFLAGS) $^ -o $@ $(LIBPATHS) $(LIBS)

# ---------------------------------------------------------------------------------

MCO_UT_SRCS=\
src/mco_ut.c \
src/tm_stm32f4_mco_output.c \
src/tm_stm32f4_gpio.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c \

mco_ut.elf: $(MCO_UT_SRCS:.c=.o) libstm32f4.a
	$(CC) $(CFLAGS) $^ -o $@ $(LIBPATHS) $(LIBS)

# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------

SM2000_RXDEMO_SRCS=\
src/sm2000_rxdemo.c \
src/sm1000_leds_switches.c \
src/debugblinky.c \
src/new_i2c.c \
src/si53xx.c \
src/stm32f4_dac.c \
src/stm32f4_adc.c \
../src/fifo.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c 

SM2000_RXDEMO_SRCS+=$(CODEC2_SRCS)

sm2000_rxdemo.elf: $(SM2000_RXDEMO_SRCS:.c=.O3.o) libstm32f4.a
	$(CC) $(CFLAGS) $^ -o $@ $(LIBPATHS) $(LIBS)

# ---------------------------------------------------------------------------------

SM2000_STW_SRCS=\
src/sm2000_stw.c \
src/sm1000_leds_switches.c \
src/debugblinky.c \
src/new_i2c.c \
src/si53xx.c \
src/stm32f4_dac.c \
src/stm32f4_adc.c \
../src/fifo.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c 

SM2000_STW_SRCS+=$(CODEC2_SRCS)

#SM2000_STW_SRCS+=$(USB_VCP)

sm2000_stw.elf: $(SM2000_STW_SRCS:.c=.O3.o) libstm32f4.a
	$(CC) $(CFLAGS) $^ -o $@ $(LIBPATHS) $(LIBS)

# ---------------------------------------------------------------------------------

SM2000_ADCDUMP_SRCS=\
src/sm2000_adc_dump.c \
src/sm1000_leds_switches.c \
src/debugblinky.c \
src/new_i2c.c \
src/si53xx.c \
src/stm32f4_dac.c \
src/stm32f4_usb_vcp.c \
src/stm32f4_adc.c \
../src/fifo.c \
src/system_stm32f4xx.c \
src/startup_stm32f4xx.s \
src/init.c 

SM2000_ADCDUMP_SRCS_SRCS+=$(CODEC2_SRCS)

SM2000_ADCDUMP_SRCS+=$(USB_VCP)

sm2000_adcdump.elf: $(SM2000_ADCDUMP_SRCS:.c=.O3.o) libstm32f4.a
	$(CC) $(CFLAGS) $^ -o $@ $(LIBPATHS) $(LIBS)

# ---------------------------------------------------------------------------------

# Objects that require the peripheral library
src/sm1000_main.o: $(PERIPHLIBDIR)/.unpack
src/codec2_profile.o: $(PERIPHLIBDIR)/.unpack

# ---------------------------------------------------------------------------------

clean:
	rm -f *.elf *.bin
	rm -f libstm32f4.a
	find . ../src -type f -name '*.o' | xargs rm -f
