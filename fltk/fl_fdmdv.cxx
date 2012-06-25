/*
  fl_fdmdv.cxx
  Created 14 June 2012
  David Rowe

  Fltk 1.3 based GUI program to prototype spectrum, waterfall, and other
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
#define DT                   0.02  // time between samples 
#define FS                8000

#define SCATTER_MEM       (FDMDV_NSYM)*50
#define SCATTER_X_MAX        3.0
#define SCATTER_Y_MAX        3.0

#define W                  1200
#define W3                 (W/3)
#define H                  600
#define H2                 (H/2)
#define SP                  20


class Spectrum;
class Waterfall;
class Scatter;

char         *fin_name = NULL;
FILE         *fin = NULL;
struct FDMDV *fdmdv;
Fl_Group     *agroup;
Fl_Window    *window;
Spectrum     *aSpectrum;
Waterfall    *aWaterfall;
Scatter      *aScatter;

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
    Spectrum(int x, int y, int w, int h): Fl_Box(x, y, w, h, "Spectrum")
    {
	align(FL_ALIGN_TOP);
	labelsize(10);
    };

};


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

    int       prev_w, prev_h;
    unsigned *pixel_buf;
    unsigned  heatmap_lut[256];
    int       greyscale;

    void new_pixel_buf(int w, int h) {
	int buf_sz, i;

	prev_w = w; prev_h = h;
	buf_sz = h*w;
	pixel_buf = new unsigned[buf_sz];
	for(i=0; i<buf_sz; i++)
	    pixel_buf[i] = 0;
    }
    
    // map val to a rgb colour
    // from http://eddiema.ca/2011/01/21/c-sharp-heatmaps/

    unsigned heatmap(float val, float min, float max) {
	unsigned r = 0;
	unsigned g = 0;
	unsigned b = 0;

	val = (val - min) / (max - min);
	if(val <= 0.2) {
	    b = (unsigned)((val / 0.2) * 255);
	} else if(val >  0.2 &&  val <= 0.7) {
	    b = (unsigned)((1.0 - ((val - 0.2) / 0.5)) * 255);
	}
	if(val >= 0.2 &&  val <= 0.6) {
	    g = (unsigned)(((val - 0.2) / 0.4) * 255);
	} else if(val >  0.6 &&  val <= 0.9) {
	    g = (unsigned)((1.0 - ((val - 0.6) / 0.3)) * 255);
	}
	if(val >= 0.5) {
	    r = (unsigned)(((val - 0.5) / 0.5) * 255);
	}
    
	//printf("%f %x %x %x\n", val, r, g, b);

	return  (b << 16) + (g << 8) + r;
    }

    void draw() {
	float  spec_index_per_px, intensity_per_dB;
	int    px_per_sec;
	int    index, dy, dy_blocks, bytes_in_row_of_blocks, b;
	int    px, py, intensity;
	unsigned *last_row, *pdest, *psrc;

	/* detect resizing of window */

	if ((h() != prev_h) || (w() != prev_w)) {
	    delete pixel_buf;
	    new_pixel_buf(w(), h());
	}

	Fl_Box::draw();

	// determine dy, the height of one "block"

	px_per_sec = (float)h()/WATERFALL_SECS_Y;
	dy = DT*px_per_sec;

	// number of dy high blocks in spectrogram

	dy_blocks = h()/dy;

	// shift previous bit map
					       
	bytes_in_row_of_blocks = dy*w()*sizeof(unsigned);

	for(b=0; b<dy_blocks-1; b++) {
	    pdest = pixel_buf + b*w()*dy;
	    psrc  = pixel_buf + (b+1)*w()*dy;
	    memcpy(pdest, psrc, bytes_in_row_of_blocks);
	}

	// create a new row of blocks at bottom

	spec_index_per_px = (float)FDMDV_NSPEC/(float)w();
	intensity_per_dB = (float)256/(MAX_DB - MIN_DB);
	last_row = pixel_buf + dy*(dy_blocks - 1)*w();

	for(px=0; px<w(); px++) {
	    index = px*spec_index_per_px;
	    intensity = intensity_per_dB * (av_mag[index] - MIN_DB);
	    if (intensity > 255) intensity = 255;
	    if (intensity < 0) intensity = 0;

	    if (greyscale) {
		for(py=0; py<dy; py++)
		    last_row[px+py*w()] = intensity<<8;
	    }
	    else {
		for(py=0; py<dy; py++)
		    last_row[px+py*w()] = heatmap_lut[intensity];
	    }
	}

	// update bit map

	fl_draw_image((uchar*)pixel_buf, x(), y(), w(), h(), 4, 0);

    }

public:

    Waterfall(int x, int y, int w, int h): Fl_Box(x, y, w, h, "Waterfall")
    {
	float f;
	int   i;

	for(i=0; i<255; i++) {
	    heatmap_lut[i] = heatmap((float)i, 0.0, 255.0);
	}
	greyscale = 0;

	align(FL_ALIGN_TOP);
	labelsize(10);
	new_pixel_buf(w,h);
    };

    ~Waterfall() {
	delete pixel_buf;
    }
};


class Scatter: public Fl_Box {
protected:
    int  first;
    COMP mem[SCATTER_MEM];
    COMP new_samples[FDMDV_NSYM];
    int  prev_w, prev_h;

    void draw() {
	float x_scale;
	float y_scale;
	int   i, j, x1, y1;

	Fl_Box::draw();

	/* detect resizing of window */

	if ((h() != prev_h) || (w() != prev_w)) {
	    fl_color(FL_BLACK);
	    fl_rectf(x(),y(),w(),h());
	    prev_h = h(); prev_w = w();
	}

	fl_push_clip(x(),y(),w(),h());

	x_scale = w()/SCATTER_X_MAX;
	y_scale = h()/SCATTER_Y_MAX;

	// erase last samples

	fl_color(FL_BLACK);
	for(i=0; i<FDMDV_NSYM; i++) {
	    x1 = x_scale * mem[i].real + x() + w()/2;
	    y1 = y_scale * mem[i].imag + y() + h()/2;
	    fl_point(x1, y1);
	    mem[i] = mem[i+FDMDV_NSYM];
	}

	// shift memory

	for(i=FDMDV_NSYM; i<SCATTER_MEM-FDMDV_NSYM; i++) {
	    mem[i] = mem[i+FDMDV_NSYM];
	}

	// draw new samples

	fl_color(FL_GREEN);
	for(i=SCATTER_MEM-FDMDV_NSYM, j=0; i<SCATTER_MEM; i++,j++) {
	    x1 = x_scale * new_samples[j].real + x() + w()/2;
	    y1 = y_scale * new_samples[j].imag + y() + h()/2;
	    fl_point(x1, y1);
	    mem[i] = new_samples[j];
	}
	fl_pop_clip();
    }

public:
    Scatter(int x, int y, int w, int h): Fl_Box(x, y, w, h, "Scatter")
    {
	int i;

	align(FL_ALIGN_TOP);
	labelsize(10);

	for(i=0; i<SCATTER_MEM; i++) {
	    mem[i].real = 0.0;
	    mem[i].imag = 0.0;
	}

	prev_w = 0; prev_h = 0;
    };

    void add_new_samples(COMP samples[]) {
	int i;

	for(i=0; i<FDMDV_NSYM; i++)
	    new_samples[i] = samples[i];
    }

};

// update average of each spectrum point
    
void new_data(float mag_dB[]) {
    int i;

    for(i=0; i<FDMDV_NSPEC; i++)
	av_mag[i] = (1.0 - BETA)*av_mag[i] + BETA*mag_dB[i];
}

// simulates real time operation by reading a raw file then pausing

void idle(void*) {
    struct FDMDV_STATS stats;
    int   rx_bits[FDMDV_BITS_PER_FRAME];
    int   sync_bit;
    int   nin, nin_prev;
    short rx_fdm_scaled[FDMDV_MAX_SAMPLES_PER_FRAME];
    float rx_fdm[FDMDV_MAX_SAMPLES_PER_FRAME];
    float rx_spec[FDMDV_NSPEC];
    int   i;

    nin = FDMDV_NOM_SAMPLES_PER_FRAME;

    if (fread(rx_fdm_scaled, sizeof(short), nin, fin) == nin) {
	Ts += (float)nin/FS;
	
	// demod per frame processing

	for(i=0; i<nin; i++)
	    rx_fdm[i] = (float)rx_fdm_scaled[i]/FDMDV_SCALE;
	nin_prev = nin;
	fdmdv_demod(fdmdv, rx_bits, &sync_bit, rx_fdm, &nin);

	// get stats and spectrum and update

	fdmdv_get_demod_stats(fdmdv, &stats);
	fdmdv_get_rx_spectrum(fdmdv, rx_spec, rx_fdm, nin_prev);
	new_data(rx_spec);
	aScatter->add_new_samples(stats.rx_symbols);

	// update plots every DT

	if (Ts >= DT) {
	    Ts -= DT;
	    aSpectrum->redraw();
	    aWaterfall->redraw();
	    aScatter->redraw();
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

    window = new Fl_Window(W, SP+H2+SP+SP+H2+SP, "fl_fmdv");
    window->size_range(100, 100);
    window->resizable();
    aSpectrum = new Spectrum(SP, SP, 2*W3-2*SP, H2);
    aWaterfall = new Waterfall(SP, SP+H2+SP+SP, 2*W3-2*SP, H2);
    aScatter = new Scatter(2*W3, SP, W3, H2);
    fdmdv = fdmdv_create();

    Fl::add_idle(idle);

    window->end();

    window->show(argc, argv);
    ret = Fl::run();

    fdmdv_destroy(fdmdv);
    fclose(fin);

    return ret;
}
