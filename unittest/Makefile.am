AM_CFLAGS = -I../src -fPIC -g  -DFLOATING_POINT -DVAR_ARRAYS -O2 -Wall
AUTOMAKE_OPTS = gnu
NAME = libcodec2
AM_CPPFLAGS = $(AM_CFLAGS)

noinst_PROGRAMS = genres genlsp extract vqtrain vqtrainjnd tnlp tinterp tquant vq_train_jvm scalarlsptest tfdmdv t48_8 lspsync create_interleaver tlspsens tphinterp vqtrainph

genres_SOURCES = genres.c ../src/lpc.c
genres_LDADD = $(lib_LTLIBRARIES) 
genres_LDFLAGS = $(LIBS)

genlsp_SOURCES = genlsp.c ../src/lpc.c ../src/lsp.c
genlsp_LDADD = $(lib_LTLIBRARIES) 
genlsp_LDFLAGS = $(LIBS)

extract_SOURCES = extract.c
extract_LDADD = $(lib_LTLIBRARIES) 
extract_LDFLAGS = $(LIBS)

vqtrain_SOURCES = vqtrain.c
vqtrain_LDADD = $(lib_LTLIBRARIES) 
vqtrain_LDFLAGS = $(LIBS)

vqtrainjnd_SOURCES = vqtrainjnd.c
vqtrainjnd_LDADD = $(lib_LTLIBRARIES) 
vqtrainjnd_LDFLAGS = $(LIBS)

vqtrainph_SOURCES = vqtrainph.c
vqtrainph_LDADD = $(lib_LTLIBRARIES) 
vqtrainph_LDFLAGS = $(LIBS)

vq_train_jvm_SOURCES = vq_train_jvm.c
vq_train_jvm_LDADD = $(lib_LTLIBRARIES) 
vq_train_jvm_LDFLAGS = $(LIBS)

CODEBOOKS = ../src/codebook.c ../src/codebookd.c ../src/codebookvq.c ../src/codebookjnd.c ../src/codebookdt.c ../src/codebookjvm.c ../src/codebookge.c

tnlp_SOURCES = tnlp.c ../src/sine.c ../src/nlp.c ../src/kiss_fft.c ../src/dump.c
tnlp_LDADD = $(lib_LTLIBRARIES) 
tnlp_LDFLAGS = $(LIBS)

tinterp_SOURCES = tinterp.c ../src/sine.c ../src/kiss_fft.c ../src/interp.c ../src/lpc.c ../src/lsp.c ../src/quantise.c $(CODEBOOKS) ../src/dump.c
tinterp_LDADD = $(lib_LTLIBRARIES) 
tinterp_LDFLAGS = $(LIBS)

tquant_SOURCES = tquant.c ../src/quantise.c ../src/lpc.c ../src/lsp.c ../src/dump.c ../src/kiss_fft.c $(CODEBOOKS)
tquant_LDADD = $(lib_LTLIBRARIES) 
tquant_LDFLAGS = $(LIBS)

scalarlsptest_SOURCES = scalarlsptest.c ../src/quantise.c ../src/lpc.c ../src/lsp.c ../src/dump.c ../src/kiss_fft.c $(CODEBOOKS)
scalarlsptest_LDADD = $(lib_LTLIBRARIES) 
scalarlsptest_LDFLAGS = $(LIBS)

tfdmdv_SOURCES = tfdmdv.c ../src/fdmdv.c ../src/kiss_fft.c ../src/octave.c
tfdmdv_LDADD = $(lib_LTLIBRARIES) 
tfdmdv_LDFLAGS = $(LIBS)

t48_8_SOURCES = t48_8.c ../src/fdmdv.c ../src/kiss_fft.c
t48_8_LDADD = $(lib_LTLIBRARIES) 
t48_8_LDFLAGS = $(LIBS)

lspsync_SOURCES = lspsync.c ../src/quantise.c ../src/lpc.c ../src/lsp.c ../src/dump.c ../src/kiss_fft.c \
../src/codec2.c ../src/sine.c ../src/nlp.c ../src/postfilter.c ../src/phase.c ../src/interp.c ../src/pack.c $(CODEBOOKS)
lspsync_LDADD = $(lib_LTLIBRARIES) 
lspsync_LDFLAGS = $(LIBS)

create_interleaver_SOURCES = create_interleaver.c 
create_interleaver_LDADD = $(lib_LTLIBRARIES) 
create_interleaver_LDFLAGS = $(LIBS)

tlspsens_SOURCES = tlspsens.c ../src/quantise.c ../src/lpc.c ../src/lsp.c ../src/dump.c ../src/kiss_fft.c ../src/codec2.c ../src/sine.c ../src/nlp.c ../src/pack.c ../src/interp.c ../src/postfilter.c ../src/phase.c $(CODEBOOKS)
tlspsens_LDADD = $(lib_LTLIBRARIES) 
tlspsens_LDFLAGS = $(LIBS)

tphinterp_SOURCES = tphinterp.c ../src/nlp.c ../src/sine.c ../src/kiss_fft.c ../src/dump.c
tphinterp_LDADD = $(lib_LTLIBRARIES) 
tphinterp_LDFLAGS = $(LIBS)

