AM_CFLAGS = -I../src -fPIC -g  -DFLOATING_POINT -DVAR_ARRAYS -O2
AUTOMAKE_OPTS = gnu
NAME = libcodec2
AM_CPPFLAGS = $(AM_CFLAGS)

bin_PROGRAMS = genres genlsp extract vqtrain tnlp tinterp tquant tcodec2


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

tnlp_SOURCES = tnlp.c ../src/sine.c ../src/nlp.c ../src/fft.c ../src/kiss_fft.c ../src/dump.c
tnlp_LDADD = $(lib_LTLIBRARIES) 
tnlp_LDFLAGS = $(LIBS)

tinterp_SOURCES = tinterp.c ../src/sine.c ../src/fft.c ../src/kiss_fft.c ../src/interp.c ../src/lpc.c ../src/lsp.c ../src/quantise.c ../src/codebook.c ../src/codebookd.c ../src/codebookdvq.c ../src/dump.c
tinterp_LDADD = $(lib_LTLIBRARIES) 
tinterp_LDFLAGS = $(LIBS)

tquant_SOURCES = tquant.c ../src/quantise.c ../src/lpc.c ../src/lsp.c ../src/dump.c ../src/fft.c ../src/kiss_fft.c ../src/codebook.c ../src/codebookd.c ../src/codebookdvq.c 
tquant_LDADD = $(lib_LTLIBRARIES) 
tquant_LDFLAGS = $(LIBS)

tcodec2_SOURCES = tcodec2.c ../src/quantise.c ../src/lpc.c ../src/lsp.c ../src/dump.c ../src/fft.c ../src/kiss_fft.c \
../src/codec2.c ../src/sine.c ../src/nlp.c ../src/postfilter.c ../src/phase.c ../src/interp.c \
../src/codebook.c ../src/codebookd.c ../src/codebookdvq.c ../src/pack.c
tcodec2_LDADD = $(lib_LTLIBRARIES) 
tcodec2_LDFLAGS = $(LIBS)

