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

# lspd VQ quantisers

set(CODEBOOKSVQ
    ${D}/lsp1.txt
    ${D}/lsp2.txt
    ${D}/lsp3.txt
    ${D}/lsp4.txt
    ${CMAKE_SOURCE_DIR}/unittest/lsp45678910.txt
)

set(CODEBOOKSJND
    ${D}/lsp1.txt
    ${D}/lsp2.txt
    ${D}/lsp3.txt
    ${D}/lsp4.txt
    ${CMAKE_SOURCE_DIR}/unittest/lspjnd5-10.txt
)

set(CODEBOOKSDT
    ${D}/lspdt1.txt
    ${D}/lspdt2.txt
    ${D}/lspdt3.txt
    ${D}/lspdt4.txt
    ${D}/lspdt5.txt
    ${D}/lspdt6.txt
    ${D}/lspdt7.txt
    ${D}/lspdt8.txt
    ${D}/lspdt9.txt
    ${D}/lspdt10.txt
)

set(CODEBOOKSJVM
    ${D}/lspjvm1.txt
    ${D}/lspjvm2.txt
    ${D}/lspjvm3.txt
)

set(CODEBOOKSVQANSSI
    ${D}/lspvqanssi1.txt
    ${D}/lspvqanssi2.txt
    ${D}/lspvqanssi3.txt
    ${D}/lspvqanssi4.txt
)

set(CODEBOOKSGE ${D}/gecb.txt)

# when crosscompiling import the executable targets from a file
if(CMAKE_CROSSCOMPILING)
    SET(IMPORT_EXECUTABLES "IMPORTFILE-NOTFOUND" CACHE FILEPATH "Point it to the export file from a native build")
    INCLUDE(${IMPORT_EXECUTABLES})
endif(CMAKE_CROSSCOMPILING)

# Build code generator binaries. These do not get installed.
if(NOT CMAKE_CROSSCOMPILING)
    # generate_codebook
    add_executable(generate_codebook generate_codebook.c)
    target_link_libraries(generate_codebook ${CMAKE_REQUIRED_LIBRARIES})
    # genlspdtcb
    add_executable(genlspdtcb genlspdtcb.c)
    target_link_libraries(genlspdtcb ${CMAKE_REQUIRED_LIBRARIES})
    # Make native builds available for cross-compiling.
    export(TARGETS generate_codebook genlspdtcb
        FILE ${CMAKE_BINARY_DIR}/ImportExecutables.cmake)
endif(NOT CMAKE_CROSSCOMPILING)

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

# codebookdt.c
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/codebookdt.c
    COMMAND generate_codebook lsp_cbdt ${CODEBOOKSDT} > ${CMAKE_CURRENT_BINARY_DIR}/codebookdt.c
    DEPENDS generate_codebook ${CODEBOOKSDT}
)

# codebookvq.c
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/codebookvq.c
    COMMAND generate_codebook lsp_cbvq ${CODEBOOKSVQ} > ${CMAKE_CURRENT_BINARY_DIR}/codebookvq.c
    DEPENDS generate_codebook ${CODEBOOKSVQ}
)

# codebookjnd.c
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/codebookjnd.c
    COMMAND generate_codebook lsp_cbjnd ${CODEBOOKSJND} > ${CMAKE_CURRENT_BINARY_DIR}/codebookjnd.c
    DEPENDS generate_codebook ${CODEBOOKSJND}
)

# codebookjvm.c
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/codebookjvm.c
    COMMAND generate_codebook lsp_cbjvm ${CODEBOOKSJVM} > ${CMAKE_CURRENT_BINARY_DIR}/codebookjvm.c
    DEPENDS generate_codebook ${CODEBOOKSJVM}
)

# codebookvqanssi.c
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/codebookvqanssi.c
    COMMAND generate_codebook lsp_cbvqanssi ${CODEBOOKSVQANSSI} > ${CMAKE_CURRENT_BINARY_DIR}/codebookvqanssi.c
    DEPENDS generate_codebook ${CODEBOOKSVQANSSI}
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
    fifo.c
    fdmdv.c
    kiss_fft.c
    interp.c
    lsp.c
    phase.c
    quantise.c
    pack.c
    codebook.c
    codebookd.c
    codebookvq.c
    codebookjnd.c
    codebookjvm.c
    codebookvqanssi.c
    codebookdt.c
    codebookge.c
    golay23.c
)

add_library(codec2 ${CODEC2_SRCS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(codec2 ${CMAKE_REQUIRED_LIBRARIES})
set_target_properties(codec2 PROPERTIES SOVERSION ${SOVERSION})

set(CODEC2_PUBLIC_HEADERS
    codec2.h
    codec2_fdmdv.h
    codec2_fifo.h
    golay23.h
)

install(FILES ${CODEC2_PUBLIC_HEADERS}
    DESTINATION include/codec2
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

add_executable(fdmdv_demod fdmdv_demod.c fdmdv.c kiss_fft.c octave.c)
target_link_libraries(fdmdv_demod ${CMAKE_REQUIRED_LIBRARIES})

add_executable(fdmdv_put_test_bits fdmdv_put_test_bits.c fdmdv.c kiss_fft.c)
target_link_libraries(fdmdv_put_test_bits ${CMAKE_REQUIRED_LIBRARIES})

add_executable(fdmdv_interleave fdmdv_interleave.c)
target_link_libraries(fdmdv_interleave ${CMAKE_REQUIRED_LIBRARIES})

add_executable(insert_errors insert_errors.c)
target_link_libraries(insert_errors ${CMAKE_REQUIRED_LIBRARIES})

add_executable(fec_enc fec_enc.c golay23.c)
target_link_libraries(fec_enc ${CMAKE_REQUIRED_LIBRARIES} codec2)

add_executable(fec_dec fec_dec.c golay23.c)
target_link_libraries(fec_dec ${CMAKE_REQUIRED_LIBRARIES} codec2)

install(TARGETS
    codec2
    c2demo
    c2demo
    c2enc
    c2dec
    c2sim
    fdmdv_get_test_bits
    fdmdv_mod fdmdv_demod
    fdmdv_put_test_bits
    fdmdv_interleave
    insert_errors
    fec_enc
    fec_dec
    LIBRARY DESTINATION lib${LIB_SUFFIX}
    ARCHIVE DESTINATION lib${LIB_SUFFIX}
    RUNTIME DESTINATION bin
)
