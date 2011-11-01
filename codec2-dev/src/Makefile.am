AM_CFLAGS = -I../src -fPIC -Wall -O3
AUTOMAKE_OPTS = gnu
NAME = codec2
AM_CPPFLAGS = $(AM_CFLAGS)

D=codebook

# lsp quantisers

CODEBOOKS= \
  $D/lsp1.txt \
  $D/lsp2.txt \
  $D/lsp3.txt \
  $D/lsp4.txt \
  $D/lsp5.txt \
  $D/lsp6.txt \
  $D/lsp7.txt \
  $D/lsp8.txt \
  $D/lsp9.txt \
  $D/lsp10.txt

# lspd quantisers

CODEBOOKSD= \
  $D/dlsp1.txt \
  $D/dlsp2.txt \
  $D/dlsp3.txt \
  $D/dlsp4.txt \
  $(top_srcdir)/unittest/lsp45678910.txt

# lspd VQ quantisers

CODEBOOKSVQ= \
  $D/lsp1.txt \
  $D/lsp2.txt \
  $D/lsp3.txt \
  $D/lsp4.txt \
  $(top_srcdir)/unittest/lsp45678910.txt

CODEBOOKSRES= \
  $D/lspres_centre1.txt \
  $D/lspres_bw1.txt \
  $D/lsp3.txt \
  $D/lsp4.txt

noinst_PROGRAMS = generate_codebook

codebook.$(OBJEXT): codebook.c
codebookd.$(OBJEXT): codebookd.c
codebookvq.$(OBJEXT): codebookvq.c
codebookres.$(OBJEXT): codebookres.c

codebook.lo: codebook.c

codebook.c: generate_codebook $(CODEBOOKS)
	./generate_codebook lsp_cb $(CODEBOOKS) > codebook.c

codebookd.c: generate_codebook $(CODEBOOKSD)
	./generate_codebook lsp_cbd $(CODEBOOKSD) > codebookd.c

codebookvq.c: generate_codebook $(CODEBOOKSVQ)
	./generate_codebook lsp_cbvq $(CODEBOOKSVQ) > codebookvq.c

codebookres.c: generate_codebook $(CODEBOOKSRES)
	./generate_codebook lsp_cbres $(CODEBOOKSRES) > codebookres.c

clean-local:
	-rm -f codebook.c codebookd.c codebookdvq.c codebookres.c

lib_LTLIBRARIES	= libcodec2.la
libcodec2_la_SOURCES = dump.c \
lpc.c \
nlp.c \
postfilter.c \
sine.c \
codec2.c \
fft.c \
kiss_fft.c \
interp.c \
lsp.c \
phase.c \
quantise.c \
pack.c \
codebook.c \
codebookd.c \
codebookvq.c \
codebookres.c

libcodec2_la_CFLAGS = $(AM_CFLAGS)
libcodec2_la_LDFLAGS = $(LIBS)

library_includedir = $(prefix)
library_include_HEADERS	= codec2.h \
defines.h \
kiss_fft.h\
_kiss_fft_guts.h \
fft.h \
interp.h \
lsp.h \
phase.h \
quantise.h \
comp.h \
dump.h \
lpc.h \
nlp.h  \
postfilter.h \
sine.h

bin_PROGRAMS = c2demo c2enc c2dec c2sim

c2demo_SOURCES = c2demo.c
c2demo_LDADD = $(lib_LTLIBRARIES) 
c2demo_LDFLAGS = $(LIBS)

c2enc_SOURCES = c2enc.c
c2enc_LDADD = $(lib_LTLIBRARIES) 
c2enc_LDFLAGS = $(LIBS)

c2dec_SOURCES = c2dec.c
c2dec_LDADD = $(lib_LTLIBRARIES) 
c2dec_LDFLAGS = $(LIBS)

c2sim_SOURCES = c2sim.c
c2sim_LDADD = $(lib_LTLIBRARIES) 
c2sim_LDFLAGS = $(LIBS)
