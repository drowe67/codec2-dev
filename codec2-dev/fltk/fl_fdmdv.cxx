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

#define MIN_DB             -40.0
#define MAX_DB               0.0
#define BETA                 0.1
#define MIN_HZ               0
#define MAX_HZ            4000
#define WATERFALL_SECS_Y     5    // number of seconds respresented by y axis of waterfall
#define DT                   0.1  // time between samples 
#define FS                8000

class Spectrum;
class Waterfall;

char         *fin_name = NULL;
FILE         *fin = NULL;
struct FDMDV *fdmdv;
Fl_Group     *agroup;
Fl_Window    *window;
Spectrum     *aSpectrum;
Waterfall    *aWaterfall;

float  av_mag[FDMDV_NSPEC]; // shared between a few classes

float  Ts = 0.0;

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
    Spectrum(int x, int y, int h, int w): Fl_Box(x,y,h,w, "Spectrum")
    {
	align(FL_ALIGN_TOP);
	labelsize(10);
    };

};


// question: how to map block sizes to pixels?
// need to handle: FDMDV_NSPEC across x, even when FDMDV_NSPEC < w(),
//                 could use axis as min of w(), FDMDV_NSPEC?
// set height of blocks based on update rate, number of seconds to display?
// do we take a snapshot every x mseconds of plot average?
// maybe just take a first pass
// draw graticule or axis

/*

  Notes:

  The height h() pixels represents WATERFALL_SECS_Y of data.  Every DT
  seconds we get a vector of FDMDV_NSPEC spectrum samples which we use
  to update the last row.  The height of each row is dy pixels, which
  maps to DT seconds.  We call each dy high rectangle of pixels a
  block.

*/

class Waterfall: public Fl_Box {
protected:

    uchar *pixel_buf;

    void draw() {
	float  spec_index_per_px, intensity_per_dB;
	int    px_per_sec;
	int    index, dy, dy_blocks, bytes_in_row_of_blocks, b;
	int    px, py, intensity;
	uchar *last_row, *pdest, *psrc;

	Fl_Box::draw();

	// determine dy, the height of one "block"

	px_per_sec = (float)h()/WATERFALL_SECS_Y;
	dy = DT*px_per_sec;

	// number of dy high blocks in spectrogram

	dy_blocks = h()/dy;

	// shift previous bit map
					       
	bytes_in_row_of_blocks = dy*w()*sizeof(uchar);

	for(b=0; b<dy_blocks-1; b++) {
	    pdest = pixel_buf + b*w()*dy;
	    psrc  = pixel_buf + (b+1)*w()*dy;
	    memcpy(pdest, psrc, bytes_in_row_of_blocks);
	}

	// create a new row of blocks at bottom

	spec_index_per_px = (float)FDMDV_NSPEC/(float)w();
	intensity_per_dB = (float)256/(MAX_DB - MIN_DB);
	last_row = pixel_buf + dy*(dy_blocks - 1)*w()*sizeof(uchar);

	for(px=0; px<w(); px++) {
	    index = px*spec_index_per_px;
	    intensity = intensity_per_dB * (av_mag[index] - MIN_DB);
	    if (intensity > 255) intensity = 255;
	    if (intensity < 0) intensity = 0;

	    for(py=0; py<dy; py++)
		last_row[px+py*w()] = intensity;
	}

	// update bit map

	fl_draw_image_mono(pixel_buf, x(), y(), w(), h(), 1, w());
    }

public:

    Waterfall(int x, int y, int h, int w): Fl_Box(x,y,h,w, "Waterfall")
    {
	int buf_sz, i;

	align(FL_ALIGN_TOP);
	labelsize(10);

	buf_sz = h*w;
	pixel_buf = new uchar[buf_sz];
	for(i=0; i<buf_sz; i++)
	    pixel_buf[i] = 0;
    };

    ~Waterfall() {
	delete pixel_buf;
    }
};


// update average of each spectrum point
    
void new_data(float mag_dB[]) {
    int i;

    for(i=0; i<FDMDV_NSPEC; i++)
	av_mag[i] = (1.0 - BETA)*av_mag[i] + BETA*mag_dB[i];
}

void idle(void*) {
    int   nin = FDMDV_NOM_SAMPLES_PER_FRAME;
    short rx_fdm_scaled[FDMDV_MAX_SAMPLES_PER_FRAME];
    float rx_fdm[FDMDV_MAX_SAMPLES_PER_FRAME];
    float rx_spec[FDMDV_NSPEC];
    int   i;

    if (fread(rx_fdm_scaled, sizeof(short), nin, fin) == nin) {
	Ts += (float)nin/FS;
	//printf("Ts %f\n", Ts);
	for(i=0; i<nin; i++)
	    rx_fdm[i] = (float)rx_fdm_scaled[i]/FDMDV_SCALE;
	fdmdv_get_rx_spectrum(fdmdv, rx_spec, rx_fdm, nin);
	new_data(rx_spec);
	if (Ts >= DT) {
	    Ts -= DT;
	    aSpectrum->redraw();
	    aWaterfall->redraw();
	}
    }
    usleep(20000);
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
    int i;

    i = 1;
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
    
    for(i=0; i<FDMDV_NSPEC; i++)
	av_mag[i] = -40.0;

    // reccommended to prevent dithering and stopped display being
    // covered by black flickering squares

    Fl::visual(FL_RGB);

    window = new Fl_Window(800, 20+300+20+20+300+20, "fl_fmdv");
    window->size_range(400,200);
    aSpectrum = new Spectrum(20, 20, 800-40, 300);
    aWaterfall = new Waterfall(20, 20+300+20+20, 800-40, 300);
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
