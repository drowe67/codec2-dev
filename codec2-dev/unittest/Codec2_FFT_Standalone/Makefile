
CFLAGS =  -Wall -g -O3 -Wno-strict-overflow -std=gnu11 -fPIC 

LDLIBS = -lm


OBJS=codec2.o sine.o codec2_fft.o kiss_fft.o

tst_codec2_fft_init: tst_codec2_fft_init.o $(OBJS) 


tst_codec2_fft_init.o: tst_codec2_fft_init.c codec2_fft.h 

sine.o: sine.c sine.h

