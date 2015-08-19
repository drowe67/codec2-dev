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
    target_link_libraries(generate_codebook ${CMAKE_REQUIRED_LIBRARIES})
    # genlspdtcb
    add_executable(genlspdtcb genlspdtcb.c)
    target_link_libraries(genlspdtcb ${CMAKE_REQUIRED_LIBRARIES})
    # Make native builds available for cross-compiling.
    export(TARGETS generate_codebook genlspdtcb
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
    cohpsk.c
    fifo.c
    fdmdv.c
    fm.c
    kiss_fft.c
    linreg.c
    interp.c
    lsp.c
    phase.c
    quantise.c
    pack.c
    codebook.c
    codebookd.c
    codebookjvm.c
    codebookmel.c
    codebooklspmelvq.c
    codebookge.c
    golay23.c
    freedv_api.c
    varicode.c
    modem_stats.c
)

set(CODEC2_PUBLIC_HEADERS
    golay23.h
    codec2.h
    codec2_fdmdv.h
    codec2_cohpsk.h
    codec2_fm.h
    codec2_fifo.h
    comp.h
    comp_prim.h
    modem_stats.h
    kiss_fft.h
    freedv_api.h
    varicode.h
    freedv_api_internal.h
)

#
# Setup the codec2 library
#
add_library(codec2 ${CODEC2_SRCS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(codec2 ${CMAKE_REQUIRED_LIBRARIES})
set_target_properties(codec2 PROPERTIES
    SOVERSION ${SOVERSION}
    PUBLIC_HEADER "${CODEC2_PUBLIC_HEADERS}"
)

add_executable(c2demo c2demo.c)
target_link_libraries(c2demo ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(c2enc c2enc.c)
target_link_libraries(c2enc ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(c2dec c2dec.c)
target_link_libraries(c2dec ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(c2sim c2sim.c ampexp.c phaseexp.c)
target_link_libraries(c2sim ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(fdmdv_get_test_bits fdmdv_get_test_bits.c fdmdv.c kiss_fft.c)
target_link_libraries(fdmdv_get_test_bits ${CMAKE_REQUIRED_LIBRARIES})

add_executable(fdmdv_mod fdmdv_mod.c fdmdv.c kiss_fft.c)
target_link_libraries(fdmdv_mod ${CMAKE_REQUIRED_LIBRARIES})

add_executable(fdmdv_demod fdmdv_demod.c fdmdv.c kiss_fft.c octave.c modem_stats.c)
target_link_libraries(fdmdv_demod ${CMAKE_REQUIRED_LIBRARIES})

add_executable(fdmdv_put_test_bits fdmdv_put_test_bits.c fdmdv.c kiss_fft.c)
target_link_libraries(fdmdv_put_test_bits ${CMAKE_REQUIRED_LIBRARIES})

add_executable(fdmdv_channel fdmdv_channel.c)
target_link_libraries(fdmdv_channel ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(fdmdv_interleave fdmdv_interleave.c)
target_link_libraries(fdmdv_interleave ${CMAKE_REQUIRED_LIBRARIES})

add_executable(insert_errors insert_errors.c)
target_link_libraries(insert_errors ${CMAKE_REQUIRED_LIBRARIES})

add_executable(fec_enc fec_enc.c golay23.c)
target_link_libraries(fec_enc ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(fec_dec fec_dec.c golay23.c)
target_link_libraries(fec_dec ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(freedv_tx freedv_tx.c)
target_link_libraries(freedv_tx ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(freedv_rx freedv_rx.c)
target_link_libraries(freedv_rx ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(fm_demod fm_demod.c fm.c)
target_link_libraries(fm_demod ${CMAKE_REQUIRED_LIBRARIES})

add_executable(cohpsk_mod cohpsk_mod.c)
target_link_libraries(cohpsk_mod ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(cohpsk_demod cohpsk_demod.c octave.c)
target_link_libraries(cohpsk_demod ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(cohpsk_get_test_bits cohpsk_get_test_bits.c)
target_link_libraries(cohpsk_get_test_bits ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(cohpsk_put_test_bits cohpsk_put_test_bits.c octave.c)
target_link_libraries(cohpsk_put_test_bits ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(cohpsk_ch cohpsk_ch.c)
target_link_libraries(cohpsk_ch ${CMAKE_REQUIRED_LIBRARIES} codec2)

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
    c2demo
    c2demo
    c2enc
    c2dec
    c2sim
    fdmdv_get_test_bits
    fdmdv_mod fdmdv_demod
    fm_demod
    fdmdv_put_test_bits
    fdmdv_interleave
    insert_errors
    fec_enc
    fec_dec
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
