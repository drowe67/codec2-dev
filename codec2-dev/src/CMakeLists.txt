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
)

add_library(codec2 ${CODEC2_SRCS})
target_link_libraries(codec2 ${CMAKE_REQUIRED_LIBRARIES})
set_target_properties(codec2 PROPERTIES SOVERSION ${SOVERSION})

set(CODEC2_PUBLIC_HEADERS
	codec2.h
	codec2_fdmdv.h
	codec2_fifo.h
	comp.h
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
