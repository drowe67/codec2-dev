add_definitions(-DFLOATING_POINT -DVAR_ARRAYS)
include_directories(../src)

add_executable(genres genres.c ../src/lpc.c)
target_link_libraries(genres codec2)

add_executable(genlsp genlsp.c ../src/lpc.c ../src/lsp.c)
target_link_libraries(genlsp codec2)

add_executable(extract extract.c)
target_link_libraries(extract codec2)

add_executable(vqtrain vqtrain.c)
target_link_libraries(vqtrain codec2)

add_executable(vqtrainjnd vqtrainjnd.c)
target_link_libraries(vqtrainjnd codec2)

add_executable(vqtrainph vqtrainph.c)
target_link_libraries(vqtrainph codec2)

add_executable(vqtrainsp vqtrainsp.c)
target_link_libraries(vqtrainsp codec2)

add_executable(genphdata genphdata.c)
target_link_libraries(genphdata codec2)

add_executable(genampdata genampdata.c)
target_link_libraries(genampdata codec2)

add_executable(polar2rect polar2rect.c)
target_link_libraries(polar2rect codec2)

add_executable(vq_train_jvm vq_train_jvm.c)
target_link_libraries(vq_train_jvm codec2)

set(CODEBOOKS ../src/codebook.c ../src/codebookd.c ../src/codebookvq.c ../src/codebookjnd.c ../src/codebookdt.c ../src/codebookjvm.c ../src/codebookvqanssi.c ../src/codebookge.c)

add_executable(tnlp tnlp.c ../src/sine.c ../src/nlp.c ../src/kiss_fft.c ../src/dump.c)
target_link_libraries(tnlp codec2)

add_executable(tinterp tinterp.c ../src/sine.c ../src/kiss_fft.c ../src/interp.c ../src/lpc.c ../src/lsp.c ../src/quantise.c ${CODEBOOKS} ../src/dump.c)
target_link_libraries(tinterp codec2)

add_executable(tquant tquant.c ../src/quantise.c ../src/lpc.c ../src/lsp.c ../src/dump.c ../src/kiss_fft.c ${CODEBOOKS})
target_link_libraries(tquant codec2)

add_executable(scalarlsptest scalarlsptest.c ../src/quantise.c ../src/lpc.c ../src/lsp.c ../src/dump.c ../src/kiss_fft.c ${CODEBOOKS})
target_link_libraries(scalarlsptest codec2)

add_executable(tfdmdv tfdmdv.c ../src/fdmdv.c ../src/kiss_fft.c ../src/octave.c)
target_link_libraries(tfdmdv codec2)

add_executable(tcohpsk tcohpsk.c ../src/cohpsk.c ../src/octave.c)
target_link_libraries(tcohpsk codec2)

add_executable(test_cohpsk_ch test_cohpsk_ch.c ../src/cohpsk.c ../src/octave.c)
target_link_libraries(test_cohpsk_ch codec2)

add_executable(t16_8 t16_8.c ../src/fdmdv.c ../src/kiss_fft.c)
target_link_libraries(t16_8 codec2)

add_executable(t16_8_short t16_8_short.c ../src/fdmdv.c ../src/kiss_fft.c)
target_link_libraries(t16_8_short codec2)

add_executable(tfsk tfsk.c ../src/kiss_fft.c ../src/kiss_fftr.c ../src/octave.c ../src/modem_probe.c)
target_link_libraries(tfsk m)

add_executable(tfmfsk tfmfsk.c ../src/octave.c ../src/modem_probe.c)
target_link_libraries(tfmfsk m)

add_executable(tdeframer tdeframer.c)
target_link_libraries(tdeframer m codec2)

add_executable(tfreedv_data_channel tfreedv_data_channel.c)
target_link_libraries(tfreedv_data_channel codec2)

#add_executable(t48_8 t48_8.c ../src/fdmdv.c ../src/kiss_fft.c)
#target_link_libraries(t48_8 codec2)

add_executable(lspsync lspsync.c ../src/quantise.c ../src/lpc.c ../src/lsp.c ../src/dump.c ../src/kiss_fft.c ../src/codec2.c ../src/sine.c ../src/nlp.c ../src/postfilter.c ../src/phase.c ../src/interp.c ../src/pack.c ${CODEBOOKS})
target_link_libraries(lspsync codec2)

add_executable(create_interleaver create_interleaver.c)
target_link_libraries(create_interleaver codec2)

add_executable(tlspsens tlspsens.c ../src/quantise.c ../src/lpc.c ../src/lsp.c ../src/dump.c ../src/kiss_fft.c ../src/codec2.c ../src/sine.c ../src/nlp.c ../src/pack.c ../src/interp.c ../src/postfilter.c ../src/phase.c ${CODEBOOKS})
target_link_libraries(tlspsens codec2)

add_executable(tprede tprede.c ../src/lpc.c)
target_link_libraries(tprede codec2)

add_executable(pre pre.c ../src/lpc.c)
target_link_libraries(pre codec2)

add_executable(de de.c ../src/lpc.c)
target_link_libraries(de codec2)

add_executable(tfifo tfifo.c ../src/fifo.c)
target_link_libraries(tfifo codec2 ${CMAKE_THREAD_LIBS_INIT})

add_executable(speexnoisesup speexnoisesup.c)
target_link_libraries(speexnoisesup ${SPEEXDSP_LIBRARY})
set_target_properties(speexnoisesup
    PROPERTIES INCLUDE_DIRECTORIES ${SPEEXDSP_INCLUDE_DIR}
)

add_executable(fdmdv_mem fdmdv_mem.c)

add_executable(raw2h raw2h.c)
target_link_libraries(raw2h codec2)

add_definitions(-D__UNITTEST__)
add_executable(c2validate c2validate.c)
target_link_libraries(c2validate codec2)

add_executable(tnewamp1 tnewamp1.c ../src/quantise.c ../src/newamp1.c ../src/mbest.c ../src/kiss_fft.c ../src/sine.c ../src/nlp.c ../src/dump.c ../src/octave.c ${CODEBOOKS})
target_link_libraries(tnewamp1 codec2)

add_executable(tsrc tsrc.c)
target_link_libraries(tsrc samplerate)

add_executable(tlininterp tlininterp.c)
add_executable(tdec tdec.c)
