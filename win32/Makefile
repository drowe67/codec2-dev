# Win32 Makefile derived from the fine work of Peter Lawrence:
#   http://code.google.com/p/codec2demo/
#
# David Rowe June 2012
#
# This makefile assumes the Linux version has been compiled first to generate
# the machine-generated codebook source files.  It is designed to be cross compiled on
# Linux.  
# 

# uncomment under Windows
#CC_PREFIX = 
# uncomment under Linux (to cross-compile for Windows)
CC_PREFIX = i586-mingw32msvc-

DLLNAME = codec2.dll
LIBNAME = codec2.lib

CFLAGS  = -O2 -Wall -D__CODEC2_WIN32__

HDRS = ../src/codec2.h ../src/codec2_fdmdv.h

OBJS = dump.o lpc.o nlp.o postfilter.o sine.o codec2.o fdmdv.o kiss_fft.o interp.o lsp.o phase.o quantise.o pack.o codebook.o codebookd.o codebookvq.o codebookjnd.o codebookjvm.o codebookdt.o codebookge.o

all: $(DLLNAME) c2demo.exe c2enc.exe c2dec.exe fdmdv_get_test_bits.exe fdmdv_mod.exe fdmdv_demod.exe fdmdv_put_test_bits.exe

$(DLLNAME): $(OBJS)
	$(CC_PREFIX)gcc -shared -o $(DLLNAME) $(OBJS) -Wl,--out-implib,$(LIBNAME)

%.o: ../src/%.c $(HDRS) Makefile
	$(CC_PREFIX)gcc $(CFLAGS) -D__CODEC2_BUILDING_DLL__ -c $< -o $@

c2demo.exe: ../src/c2demo.c  $(DLLNAME) $(LIBNAME)
	$(CC_PREFIX)gcc $(CFLAGS) $< -o $@ $(LIBNAME)

c2enc.exe: ../src/c2enc.c $(DLLNAME) $(LIBNAME)
	$(CC_PREFIX)gcc $(CFLAGS) $< -o $@ $(LIBNAME)

c2dec.exe: ../src/c2dec.c $(DLLNAME) $(LIBNAME)
	$(CC_PREFIX)gcc $(CFLAGS) $< -o $@ $(LIBNAME)

fdmdv_get_test_bits.exe: ../src/fdmdv_get_test_bits.c $(DLLNAME) $(LIBNAME)
	$(CC_PREFIX)gcc $(CFLAGS) $< -o $@ $(LIBNAME)

fdmdv_mod.exe: ../src/fdmdv_mod.c $(DLLNAME) $(LIBNAME)
	$(CC_PREFIX)gcc $(CFLAGS) $< -o $@ $(LIBNAME)

fdmdv_demod.exe: ../src/fdmdv_demod.c ../src/octave.c $(DLLNAME) $(LIBNAME)
	$(CC_PREFIX)gcc $(CFLAGS) $< ../src/octave.c -o $@ $(LIBNAME)

fdmdv_put_test_bits.exe: ../src/fdmdv_put_test_bits.c $(DLLNAME) $(LIBNAME)
	$(CC_PREFIX)gcc $(CFLAGS) $< -o $@ $(LIBNAME)

clean:
	rm -f $(DLLNAME) $(LIBNAME)
	rm -f *.o *.exe

