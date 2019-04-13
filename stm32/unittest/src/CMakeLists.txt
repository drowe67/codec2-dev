set(STM32F4_SYSTEM_SRCS
../../src/system_stm32f4xx.c
../../src/stm32f4_machdep.c
startup_stm32f4xx.s
)

list(APPEND SEMIHOSTING_SRCS semihosting.c ${STM32F4_SYSTEM_SRCS})
list(APPEND SEMIHOSTING_PROFILE_LIBS codec2_prof stm32f4 CMSIS rdimon)
list(APPEND SEMIHOSTING_LIBS codec2 stm32f4 CMSIS rdimon)

macro(profiledSemihostedBin target)
	add_mapped_executable(${target} ${target}.c ${SEMIHOSTING_SRCS})
	target_link_libraries(${target} ${SEMIHOSTING_PROFILE_LIBS})
	target_compile_definitions(${target} PRIVATE "-DPROFILE -DSEMIHOST_USE_STDIO --specs=rdimon.specs")
	elf2bin(${target})
endmacro()

macro(semihostedBin target)
	add_mapped_executable(${target} ${target}.c ${SEMIHOSTING_SRCS})
	target_link_libraries(${target} ${SEMIHOSTING_LIBS})
	target_compile_definitions(${target} PRIVATE "-DSEMIHOST_USE_STDIO --specs=rdimon.specs")
	elf2bin(${target})
endmacro()

semihostedBin(tst_api_tx)
semihostedBin(tst_codec2_enc)
semihostedBin(tst_codec2_dec)
semihostedBin(tst_ofdm_mod)
profiledSemihostedBin(tst_ofdm_demod) 
semihostedBin(tst_ldpc_enc)
semihostedBin(tst_ldpc_dec)
semihostedBin(tst_api_mod)
semihostedBin(tst_api_demod)
semihostedBin(tst_semihost) 
semihostedBin(tst_codec2_fft_init)
