set(STM32F4_SYSTEM_SRCS
../../src/system_stm32f4xx.c
../../src/init.c
startup_stm32f4xx.s
stm32f4_machdep.c
)

list(APPEND SEMIHOSTING_SRCS semihosting.c ${STM32F4_SYSTEM_SRCS})
list(APPEND SEMIHOSTING_PROFILE_LIBS codec2_prof stm32f4 CMSIS rdimon)

macro(profiledSemihostedBin target)
	add_executable(${target} ${target}.c ${SEMIHOSTING_SRCS})
	target_link_libraries(${target} ${SEMIHOSTING_PROFILE_LIBS})
	target_compile_definitions(${target} PRIVATE "-DPROFILE -DSEMIHOST_USE_STDIO --specs=rdimon.specs")
	elf2bin(${target})
endmacro()

profiledSemihostedBin(tst_api_tx)
profiledSemihostedBin(tst_codec2_enc)
profiledSemihostedBin(tst_codec2_dec)
profiledSemihostedBin(tst_ofdm_mod)
profiledSemihostedBin(tst_ofdm_demod) 
profiledSemihostedBin(tst_ldpc_enc)
profiledSemihostedBin(tst_ldpc_dec)
profiledSemihostedBin(tst_api_mod)
profiledSemihostedBin(tst_api_demod)
profiledSemihostedBin(tst_semihost) 
profiledSemihostedBin(tst_codec2_fft_init)
