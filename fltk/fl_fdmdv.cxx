/*
  fl_fdmdv.cxx
  Created 14 June 2012
  David Rowe

  Fltk based GUI program to prototype spectrum, waterfall, and other
  FDMDV GUI displays.
*/

#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Box.H>
#include <FL/fl_draw.H>
#include <FL/Fl_Group.H>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "fdmdv.h"

#define MIN_DB -40.0
#define MAX_DB   0.0
#define BETA     0.1
#define MIN_HZ     0
#define MAX_HZ  4000

class Spectrum;
char         *fin_name = NULL;
FILE         *fin = NULL;
struct FDMDV *fdmdv;
Fl_Group     *agroup;
Fl_Window    *window;
Spectrum     *aSpectrum;


class Spectrum: public Fl_Box {
protected:
    void draw() {
	float x_px_per_point = 0.0;
	float y_px_per_dB = 0.0;
	int   i, x1, y1, x2, y2, len;
	float mag1, mag2;
	char  label[20];
	float px_per_hz;

	Fl_Box::draw();
	fl_color(FL_BLACK);
	fl_rectf(x(),y(),w(),h());
	fl_color(FL_GREEN);
	fl_line_style(FL_SOLID);

	fl_push_clip(x(),y(),w(),h());
	//printf("%d %d\n", w(), h());
	x_px_per_point = (float)w()/FDMDV_NSPEC;
	y_px_per_dB = (float)h()/(MAX_DB - MIN_DB);

	// plot spectrum

	for(i=0; i<FDMDV_NSPEC-1; i++) {
	    mag1 = av_mag[i];
	    mag2 = av_mag[i+1];

	    x1 = x() + i*x_px_per_point;
	    y1 = y() + -mag1*y_px_per_dB;
	    x2 = x() + (i+1)*x_px_per_point;
	    y2 = y() + -mag2*y_px_per_dB;
	    fl_line(x1,y1,x2,y2);   
	}

	// y axis graticule

	fl_line_style(FL_DOT);
	for(i=MIN_DB; i<MAX_DB; i+=10) {
	    x1 = x();
	    y1 = y() + -i*y_px_per_dB;
	    x2 = x() + w();
	    y2 = y1;
	    //printf("%d %d %d %d\n", x1, y1, x2, y2);
	    fl_line(x1,y1,x2,y2);   
	    sprintf(label, "%d", i);
	    fl_draw(label, x1, y1);
	}

	// x axis graticule

	px_per_hz = (float)w()/(MAX_HZ-MIN_HZ);
	fl_line_style(FL_DOT);
	for(i=500; i<MAX_HZ; i+=500) {
	    x1 = x() + i*px_per_hz;
	    y1 = y();
	    x2 = x1;
	    y2 = y() + h();
	    //printf("i=%d %d %d %d %d\n", i, x1, y1, x2, y2);
	    fl_line(x1,y1,x2,y2);   
	    sprintf(label, "%d", i);
	    fl_draw(label, x1, y2);
	}

	fl_pop_clip();
    }

public:
    float  av_mag[FDMDV_NSPEC];

    Spectrum(int x, int y, int h, int w): Fl_Box(x,y,h,w, "Spectrum")
    {
	int i;

	for(i=0; i<FDMDV_NSPEC; i++)
	    av_mag[i] = -40.0;

	align(FL_ALIGN_TOP);
	labelsize(10);
    };

    // update average of each spectrum point
    
    void new_data(float mag_dB[]) {
	int i;

	for(i=0; i<FDMDV_NSPEC; i++)
	    av_mag[i] = (1.0 - BETA)*av_mag[i] + BETA*mag_dB[i];
    }

};

void idle(void*) {
    int   nin = FDMDV_NOM_SAMPLES_PER_FRAME;
    short rx_fdm_scaled[FDMDV_MAX_SAMPLES_PER_FRAME];
    float rx_fdm[FDMDV_MAX_SAMPLES_PER_FRAME];
    float rx_spec[FDMDV_NSPEC];
    int   i;

    if (fread(rx_fdm_scaled, sizeof(short), nin, fin) == nin) {
	for(i=0; i<nin; i++)
	    rx_fdm[i] = (float)rx_fdm_scaled[i]/FDMDV_SCALE;
	fdmdv_get_rx_spectrum(fdmdv, rx_spec, rx_fdm, nin);
	aSpectrum->new_data(rx_spec);
	aSpectrum->redraw();
	usleep(20000);
    }
}

int arg_callback(int argc, char **argv, int &i) {
    if (argv[i][1] == 'i') {
	if ((i+1) >= argc) 
	    return 0;
	fin_name = argv[i+1];
	i += 2;
	return 2;
    }
    return 0;
}

int main(int argc, char **argv) {
    int ret;
    int i = 1;

    Fl::args(argc,argv,i,arg_callback);
    
    if (argc != 3) {
	printf("usage: %s -i inputFdmdvRawFile\n", argv[0]);
	exit(0);
    }

    fin = fopen(fin_name,"rb");
    if (fin == NULL) {
	fprintf(stderr, "Error opening input fdmdv raw file %s\n", argv[1]);
	exit(1);
    }
    
    window = new Fl_Window(800, 600, "fl_fmdv");
    window->size_range(400,200);
    aSpectrum = new Spectrum(20, 20, 800-40, 300);
    window->add_resizable(*aSpectrum);
    fdmdv = fdmdv_create();

    Fl::add_idle(idle);

    window->end();

    window->show(argc, argv);
    ret = Fl::run();

    fdmdv_destroy(fdmdv);
    fclose(fin);

    return ret;
}
