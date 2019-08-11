add_definitions(-DFLOATING_POINT -DVAR_ARRAYS)
include_directories(../src)

add_executable(tfdmdv tfdmdv.c ../src/fdmdv.c ../src/kiss_fft.c ../src/octave.c)
target_link_libraries(tfdmdv codec2)

add_executable(tcohpsk tcohpsk.c ../src/cohpsk.c ../src/octave.c)
target_link_libraries(tcohpsk codec2)

add_executable(tfsk tfsk.c ../src/kiss_fft.c ../src/kiss_fftr.c ../src/octave.c ../src/modem_probe.c)
target_link_libraries(tfsk m)

add_executable(tfmfsk tfmfsk.c ../src/octave.c ../src/modem_probe.c)
target_link_libraries(tfmfsk m)

add_executable(tdeframer tdeframer.c)
target_link_libraries(tdeframer m codec2)

add_definitions(-DMODEMPROBE_ENABLE -DXXXXX)

add_executable(tofdm tofdm.c ../src/ofdm.c ../src/octave.c ../src/kiss_fft.c ../src/modem_probe.c ../src/mpdecode_core.c ../src/phi0.c ../src/filter.c ../src/HRA_112_112.c)
target_link_libraries(tofdm m)

add_executable(tfifo tfifo.c ../src/codec2_fifo.c)
target_link_libraries(tfifo codec2 ${CMAKE_THREAD_LIBS_INIT})

add_executable(fdmdv_mem fdmdv_mem.c)

add_executable(ofdm_mem ofdm_mem.c ../src/ofdm.c ../src/octave.c ../src/kiss_fft.c ../src/modem_probe.c ../src/mpdecode_core.c ../src/phi0.c ../src/filter.c)
target_link_libraries(ofdm_mem m)

add_library(function_trace STATIC ../unittest/function_trace.c)

add_executable(ofdm_stack ofdm_stack.c ../src/ofdm.c ../src/octave.c ../src/kiss_fft.c ../src/modem_probe.c ../src/mpdecode_core.c ../src/phi0.c ../src/filter.c)
if (CMAKE_C_COMPILER MATCHES "gcc$")
    target_link_libraries(ofdm_stack function_trace m -no-pie "-Wl,-Map=ofdm_stack.map")
    target_compile_options(ofdm_stack PUBLIC -fstack-usage -finstrument-functions -no-pie)
else()
    target_link_libraries(ofdm_stack function_trace m -no-pie)
    target_compile_options(ofdm_stack PUBLIC -finstrument-functions -no-pie)
endif()
add_definitions(-D__UNITTEST__)

add_executable(tnewamp1 tnewamp1.c ../src/quantise.c ../src/newamp1.c ../src/mbest.c ../src/kiss_fft.c ../src/sine.c ../src/nlp.c ../src/dump.c ../src/octave.c ${CODEBOOKS})
target_link_libraries(tnewamp1 codec2)

add_executable(compare_ints compare_ints.c)

add_executable(compare_floats compare_floats.c)

add_executable(test_phi0 test_phi0.c ../src/phi0.c)
target_link_libraries(test_phi0 m)

add_executable(tst_codec2_fft_init tst_codec2_fft_init.c)
target_link_libraries(tst_codec2_fft_init m codec2)
