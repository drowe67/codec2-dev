# Requires FLTK 1.3 & Portaudio V19

FLTK_CFLAGS += $(shell fltk-config --ldstaticflags)
CFLAGS = -O3 -g
LIBS = ../src/.libs/libcodec2.a -lm -lrt -lportaudio -pthread
LC2POC_C = fl_fdmdv.cxx

all: fl_fdmdv

fl_fdmdv: Makefile $(LC2POC_C)
	g++ $(LC2POC_C) -I../src/ -o fl_fdmdv $(CFLAGS) $(FLTK_CFLAGS) $(LIBS)

clean:
	rm -f fl_fdmdv
