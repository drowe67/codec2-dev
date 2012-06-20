FLTK_CFLAGS += $(shell fltk-config --ldstaticflags)
CFLAGS = -O3 -g

LC2POC_C = fl_fdmdv.cxx

all: fl_fdmdv

fl_fdmdv: Makefile $(LC2POC_C)
	g++ $(LC2POC_C) -I../src/ -o fl_fdmdv $(CFLAGS) $(FLTK_CFLAGS) ../src/.libs/libcodec2.a

clean:
	rm -f fl_fdmdv
