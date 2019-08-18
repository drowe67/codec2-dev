#
# Generated sources
#

set(D ${CMAKE_CURRENT_SOURCE_DIR}/codebook)

# lsp quantisers

set(CODEBOOKS
    ${D}/lsp1.txt
    ${D}/lsp2.txt
    ${D}/lsp3.txt
    ${D}/lsp4.txt
    ${D}/lsp5.txt
    ${D}/lsp6.txt
    ${D}/lsp7.txt
    ${D}/lsp8.txt
    ${D}/lsp9.txt
    ${D}/lsp10.txt
)

# lspd quantisers

set(CODEBOOKSD
    ${D}/dlsp1.txt
    ${D}/dlsp2.txt
    ${D}/dlsp3.txt
    ${D}/dlsp4.txt
    ${D}/dlsp5.txt
    ${D}/dlsp6.txt
    ${D}/dlsp7.txt
    ${D}/dlsp8.txt
    ${D}/dlsp9.txt
    ${D}/dlsp10.txt
)

set(CODEBOOKSJVM
    ${D}/lspjvm1.txt
    ${D}/lspjvm2.txt
    ${D}/lspjvm3.txt
)

set(CODEBOOKSMEL
    ${D}/mel1.txt
    ${D}/mel2.txt
    ${D}/mel3.txt
    ${D}/mel4.txt
    ${D}/mel5.txt
    ${D}/mel6.txt
)

set(CODEBOOKSLSPMELVQ
    ${D}/lspmelvq1.txt
    ${D}/lspmelvq2.txt
    ${D}/lspmelvq3.txt
)

set(CODEBOOKSGE ${D}/gecb.txt)

set(CODEBOOKSNEWAMP1
    ${D}/train_120_1.txt
    ${D}/train_120_2.txt
)

set(CODEBOOKSNEWAMP1_ENERGY
    ${D}/newamp1_energy_q.txt
)

set(CODEBOOKSNEWAMP2
    ${D}/codes_450.txt
)

set(CODEBOOKSNEWAMP2_ENERGY
    ${D}/newamp2_energy_q.txt
)

# when crosscompiling we need a native executable
if(CMAKE_CROSSCOMPILING)
    include(ExternalProject)
    ExternalProject_Add(codec2_native
       SOURCE_DIR ${CMAKE_SOURCE_DIR}
       BUILD_COMMAND $(MAKE) generate_codebook
       INSTALL_COMMAND ${CMAKE_COMMAND} -E copy src/generate_codebook ${CMAKE_CURRENT_BINARY_DIR}
    )
    add_executable(generate_codebook IMPORTED)
    set_target_properties(generate_codebook PROPERTIES
        IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/generate_codebook)
    add_dependencies(generate_codebook codec2_native)
else(CMAKE_CROSSCOMPILING)
# Build code generator binaries. These do not get installed.
    # generate_codebook
    add_executable(generate_codebook generate_codebook.c)
    target_link_libraries(generate_codebook m ${CMAKE_REQUIRED_LIBRARIES})
    # Make native builds available for cross-compiling.
    export(TARGETS generate_codebook
        FILE ${CMAKE_BINARY_DIR}/ImportExecutables.cmake)
endif(CMAKE_CROSSCOMPILING)


# codebook.c
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/codebook.c
    COMMAND generate_codebook lsp_cb ${CODEBOOKS} > ${CMAKE_CURRENT_BINARY_DIR}/codebook.c
    DEPENDS generate_codebook ${CODEBOOKS}
)

# codebookd.c
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/codebookd.c
    COMMAND generate_codebook lsp_cbd ${CODEBOOKSD} > ${CMAKE_CURRENT_BINARY_DIR}/codebookd.c
    DEPENDS generate_codebook ${CODEBOOKSD}
)

# codebookjvm.c
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/codebookjvm.c
    COMMAND generate_codebook lsp_cbjvm ${CODEBOOKSJVM} > ${CMAKE_CURRENT_BINARY_DIR}/codebookjvm.c
    DEPENDS generate_codebook ${CODEBOOKSJVM}
)


# codebookmel.c
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/codebookmel.c
    COMMAND generate_codebook mel_cb ${CODEBOOKSMEL} > ${CMAKE_CURRENT_BINARY_DIR}/codebookmel.c
    DEPENDS generate_codebook ${CODEBOOKSMEL}
)

# codebooklspmelvq.c
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/codebooklspmelvq.c
    COMMAND generate_codebook lspmelvq_cb ${CODEBOOKSLSPMELVQ} > ${CMAKE_CURRENT_BINARY_DIR}/codebooklspmelvq.c
    DEPENDS generate_codebook ${CODEBOOKSLSPMELVQ}
)

# codebookge.c
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/codebookge.c
    COMMAND generate_codebook ge_cb ${CODEBOOKSGE} > ${CMAKE_CURRENT_BINARY_DIR}/codebookge.c
    DEPENDS generate_codebook ${CODEBOOKSGE}
)

# codebooknewamp1.c
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/codebooknewamp1.c
    COMMAND generate_codebook newamp1vq_cb ${CODEBOOKSNEWAMP1} > ${CMAKE_CURRENT_BINARY_DIR}/codebooknewamp1.c
    DEPENDS generate_codebook ${CODEBOOKSNEWAMP1}
)

# codebooknewamp1_energy.c
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/codebooknewamp1_energy.c
    COMMAND generate_codebook newamp1_energy_cb ${CODEBOOKSNEWAMP1_ENERGY} > ${CMAKE_CURRENT_BINARY_DIR}/codebooknewamp1_energy.c
    DEPENDS generate_codebook ${CODEBOOKSNEWAMP1_ENERGY}
)

# codebooknewamp2.c
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/codebooknewamp2.c
    COMMAND generate_codebook newamp2vq_cb ${CODEBOOKSNEWAMP2} > ${CMAKE_CURRENT_BINARY_DIR}/codebooknewamp2.c
    DEPENDS generate_codebook ${CODEBOOKSNEWAMP2}
)

# codebooknewamp2_energy.c
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/codebooknewamp2_energy.c
    COMMAND generate_codebook newamp2_energy_cb ${CODEBOOKSNEWAMP2_ENERGY} > ${CMAKE_CURRENT_BINARY_DIR}/codebooknewamp2_energy.c
    DEPENDS generate_codebook ${CODEBOOKSNEWAMP2_ENERGY}
)

#
# codec2 library sources
#

set(CODEC2_SRCS
    dump.c
    lpc.c
    nlp.c
    postfilter.c
    sine.c
    codec2.c
    codec2_fft.c
    cohpsk.c
    codec2_fifo.c
    fdmdv.c
    fm.c
    fsk.c
    fmfsk.c
    tdma.c
    kiss_fft.c
    kiss_fftr.c
    linreg.c
    interp.c
    lsp.c
    mbest.c
    newamp1.c
    newamp2.c
    ofdm.c
    phase.c
    quantise.c
    pack.c
    codebook.c
    codebookd.c
    codebookjvm.c
    codebookmel.c
    codebooklspmelvq.c
    codebookge.c
    codebooknewamp1.c
    codebooknewamp1_energy.c
    codebooknewamp2.c
    codebooknewamp2_energy.c
    golay23.c
    freedv_api.c
    freedv_vhf_framing.c
    freedv_data_channel.c
    varicode.c
    modem_stats.c
    ofdm.c
    mpdecode_core.c
    phi0.c
    gp_interleaver.c
    interldpc.c
    filter.c
    horus_api.c
    horus_l2.c
    HRA_112_112.c
    HRAb_396_504.c
)

set(CODEC2_PUBLIC_HEADERS
    golay23.h
    codec2.h
    codec2_fdmdv.h
    codec2_cohpsk.h
    codec2_fm.h
    codec2_ofdm.h
    fsk.h
    codec2_fifo.h
    comp.h
    comp_prim.h
    modem_stats.h
    kiss_fft.h
    freedv_api.h
    varicode.h
    freedv_api_internal.h
    filter.h
    horus_api.h
    ${CODEC2_VERSION_PATH}/version.h
)

#
# Debug options
#
#add_definitions(-DDEBUG_ALLOC)

#
# Setup the codec2 library
#
add_library(codec2 ${CODEC2_SRCS})
if(UNIX)
    target_link_libraries(codec2 m)
endif(UNIX)
if(LPCNET AND lpcnetfreedv_FOUND)
    target_link_libraries(codec2 lpcnetfreedv)
endif()
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(codec2 PROPERTIES
    SOVERSION ${SOVERSION}
    PUBLIC_HEADER "${CODEC2_PUBLIC_HEADERS}"
)

# Required include directories for export
target_include_directories(codec2 INTERFACE
    $<INSTALL_INTERFACE:include/codec2>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
)

# Export config for import from full install
install(EXPORT codec2-config
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/codec2
)

# Export config for import from build tree.
export(TARGETS codec2
    FILE ${CMAKE_BINARY_DIR}/codec2.cmake
)

add_executable(c2demo c2demo.c)
target_link_libraries(c2demo ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(c2enc c2enc.c)
target_link_libraries(c2enc ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(c2dec c2dec.c)
target_link_libraries(c2dec ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(c2sim c2sim.c)
target_link_libraries(c2sim ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(fdmdv_get_test_bits fdmdv_get_test_bits.c fdmdv.c kiss_fft.c codec2_fft.c kiss_fftr.c)
target_link_libraries(fdmdv_get_test_bits m ${CMAKE_REQUIRED_LIBRARIES})

add_executable(fdmdv_mod fdmdv_mod.c fdmdv.c kiss_fft.c codec2_fft.c kiss_fftr.c)
target_link_libraries(fdmdv_mod m ${CMAKE_REQUIRED_LIBRARIES})

add_executable(fdmdv_demod fdmdv_demod.c fdmdv.c kiss_fft.c octave.c modem_stats.c codec2_fft.c kiss_fftr.c)
target_link_libraries(fdmdv_demod m ${CMAKE_REQUIRED_LIBRARIES})

add_executable(fdmdv_put_test_bits fdmdv_put_test_bits.c fdmdv.c kiss_fft.c codec2_fft.c kiss_fftr.c)
target_link_libraries(fdmdv_put_test_bits m ${CMAKE_REQUIRED_LIBRARIES})

add_executable(fdmdv_channel fdmdv_channel.c)
target_link_libraries(fdmdv_channel ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(insert_errors insert_errors.c)
target_link_libraries(insert_errors ${CMAKE_REQUIRED_LIBRARIES})

add_executable(freedv_tx freedv_tx.c)
target_link_libraries(freedv_tx ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(freedv_rx freedv_rx.c)
target_link_libraries(freedv_rx ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(fsk_mod fsk_mod.c)
target_link_libraries(fsk_mod ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(fsk_mod_ext_vco fsk_mod_ext_vco.c)
target_link_libraries(fsk_mod_ext_vco m ${CMAKE_REQUIRED_LIBRARIES})

add_executable(fsk_demod fsk_demod.c modem_probe.c octave.c)
target_link_libraries(fsk_demod ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(fsk_get_test_bits fsk_get_test_bits.c)
target_link_libraries(fsk_get_test_bits)

add_executable(fsk_put_test_bits fsk_put_test_bits.c)
target_link_libraries(fsk_put_test_bits ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(fm_demod fm_demod.c fm.c)
target_link_libraries(fm_demod m ${CMAKE_REQUIRED_LIBRARIES})

add_executable(cohpsk_mod cohpsk_mod.c)
target_link_libraries(cohpsk_mod ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(ofdm_get_test_bits ofdm_get_test_bits.c)
target_link_libraries(ofdm_get_test_bits ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(ofdm_gen_test_bits ofdm_gen_test_bits.c)
target_link_libraries(ofdm_gen_test_bits ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(ofdm_put_test_bits ofdm_put_test_bits.c)
target_link_libraries(ofdm_put_test_bits ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(ofdm_mod ofdm_mod.c)
target_link_libraries(ofdm_mod ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(ofdm_demod ofdm_demod.c octave.c)
target_link_libraries(ofdm_demod ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(fmfsk_mod fmfsk_mod.c)
target_link_libraries(fmfsk_mod ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(fmfsk_demod fmfsk_demod.c modem_probe.c octave.c)
target_link_libraries(fmfsk_demod ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(vhf_deframe_c2 vhf_deframe_c2.c)
target_link_libraries(vhf_deframe_c2  ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(vhf_frame_c2 vhf_frame_c2.c)
target_link_libraries(vhf_frame_c2  ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(cohpsk_demod cohpsk_demod.c octave.c)
target_link_libraries(cohpsk_demod ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(cohpsk_get_test_bits cohpsk_get_test_bits.c)
target_link_libraries(cohpsk_get_test_bits ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(cohpsk_put_test_bits cohpsk_put_test_bits.c octave.c)
target_link_libraries(cohpsk_put_test_bits ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(cohpsk_ch cohpsk_ch.c)
target_link_libraries(cohpsk_ch ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(ldpc_noise ldpc_noise.c)
target_link_libraries(ldpc_noise m)

add_executable(ldpc_enc ldpc_enc.c)
target_link_libraries(ldpc_enc ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(ldpc_dec ldpc_dec.c)
target_link_libraries(ldpc_dec ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(drs232 drs232.c)
target_link_libraries(drs232 ${CMAKE_REQUIRED_LIBRARIES})

add_executable(drs232_ldpc drs232_ldpc.c)
target_link_libraries(drs232_ldpc ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_definitions(-DHORUS_L2_RX -DINTERLEAVER -DSCRAMBLER -DRUN_TIME_TABLES)
add_executable(horus_demod horus_demod.c horus_api.c horus_l2.c golay23.c fsk.c kiss_fft.c)
target_link_libraries(horus_demod m ${CMAKE_REQUIRED_LIBRARIES})

install(TARGETS codec2 EXPORT codec2-config
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/codec2
)

#install(EXPORT codec2-config
#    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/codec2
#)

install(TARGETS
    c2enc
    c2dec
    fdmdv_get_test_bits
    fdmdv_mod fdmdv_demod
    fm_demod
    fsk_mod
    fdmdv_put_test_bits
    insert_errors
    drs232
    drs232_ldpc
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
