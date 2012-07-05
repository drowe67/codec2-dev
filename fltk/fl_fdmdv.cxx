/*
  fl_fdmdv.cxx
  Created 14 June 2012
  David Rowe

  Fltk 1.3 based GUI program to prototype spectrum, waterfall, and other
  FDMDV GUI displays.
*/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Box.H>
#include <FL/fl_draw.H>
#include <FL/Fl_Group.H>
#include <FL/names.h>

#include "portaudio.h"

#include "fdmdv.h"

#define MIN_DB             -40.0 
#define MAX_DB               0.0
#define BETA                 0.1  // constant for time averageing spectrum data
#define MIN_HZ               0
#define MAX_HZ            4000
#define WATERFALL_SECS_Y     5    // number of seconds respresented by y axis of waterfall
#define DT                   0.02 // time between samples 
#define FS                8000    // FDMDV modem sample rate

#define SCATTER_MEM       (FDMDV_NSYM)*50
#define SCATTER_X_MAX        3.0
#define SCATTER_Y_MAX        3.0

// main window params

#define W                  1200
#define W3                 (W/3)
#define H                  600
#define H2                 (H/2)
#define SP                  20

#define SOUND_CARD_FS      48000

// forward class declarations

class Spectrum;
class Waterfall;
class Scatter;
class Scalar;

// Globals --------------------------------------

char         *fin_name = NULL;
char         *sound_dev_name = NULL;
FILE         *fin = NULL;
struct FDMDV *fdmdv;
float         av_mag[FDMDV_NSPEC]; // shared between a few classes

// GUI variables --------------------------------

Fl_Group     *agroup;
Fl_Window    *window;
Fl_Window    *zoomSpectrumWindow = NULL;
Fl_Window    *zoomWaterfallWindow = NULL;
Spectrum     *aSpectrum;
Spectrum     *aZoomedSpectrum;
Waterfall    *aWaterfall;
Waterfall    *aZoomedWaterfall;
Scatter      *aScatter;
Scalar       *aTimingEst;
Scalar       *aFreqEst;
Scalar       *aSNR;
int          zoom_spectrum = 0;

// Main processing loop states ------------------

float  Ts = 0.0;
int    nbuf = 0;
short  rx_fdm_scaled[7*FDMDV_NOM_SAMPLES_PER_FRAME];
int    nin = FDMDV_NOM_SAMPLES_PER_FRAME;

// Portaudio states -----------------------------

PaStreamParameters inputParameters;
PaStream *stream = NULL;
PaError err;

// Class for each window type  ------------------

class Spectrum: public Fl_Box {
protected:
    int handle(int event) {

	//  detect a left mouse down if inside the spectrum window

	if ((event == FL_NO_EVENT) && (Fl::event_button() == 1)) {
	    if ((Fl::event_x() > x()) && (Fl::event_x() < (x() + w())) &&
		(Fl::event_y() > y()) && (Fl::event_y() < (y() + h()))) {

		// show zoomed spectrum window

		zoomSpectrumWindow->show();
	    }
	    
	}
	return 0;
    }

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
    
    int handle(int event) {

	//  detect a left mouse down if inside the window

	if ((event == FL_NO_EVENT) && (Fl::event_button() == 1)) {
	    if ((Fl::event_x() > x()) && (Fl::event_x() < (x() + w())) &&
		(Fl::event_y() > y()) && (Fl::event_y() < (y() + h()))) {

		// show zoomed spectrum window

		zoomWaterfallWindow->show();
	    }
	    
	}
	return 0;
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
    COMP mem[SCATTER_MEM];
    COMP new_samples[FDMDV_NSYM];
    int  prev_w, prev_h, prev_x, prev_y;

    void draw() {
	float x_scale;
	float y_scale;
	int   i, j, x1, y1;

	Fl_Box::draw();

	/* detect resizing of window */

	if ((h() != prev_h) || (w() != prev_w) || (x() != prev_x) || (y() != prev_y)) {
	    fl_color(FL_BLACK);
	    fl_rectf(x(),y(),w(),h());
	    prev_h = h(); prev_w = w(); prev_x = x(); prev_y = y();
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

	prev_w = 0; prev_h = 0; prev_x = 0; prev_y = 0;
    };

    void add_new_samples(COMP samples[]) {
	int i;

	for(i=0; i<FDMDV_NSYM; i++)
	    new_samples[i] = samples[i];
    }

};


// general purpose way of plotting scalar values that are 
// updated once per frame

class Scalar: public Fl_Box {
protected:
    int    x_max, y_max;
    float *mem;              /* array of x_max samples */
    float  new_sample;
    int    index, step;
    int    prev_w, prev_h, prev_x, prev_y;

    int clip(int y1) {
	if (y1 > (h()/2 - 10))
	    y1 = h()/2 - 10;       
	if (y1 < -(h()/2 - 10))
	    y1 = -(h()/2 - 10);       
	return y1;
    }

    void draw() {
	float x_scale;
	float y_scale;
	int   i, j, x1, y1, x2, y2;
	char  label[100];

	Fl_Box::draw();

	/* detect resizing of window */

	if ((h() != prev_h) || (w() != prev_w) || (x() != prev_x) || (y() != prev_y)) {
	    fl_color(FL_BLACK);
	    fl_rectf(x(),y(),w(),h());
	    prev_h = h(); prev_w = w(); prev_x = x(); prev_y = y();
	}

	fl_push_clip(x(),y(),w(),h());

	x_scale = (float)w()/x_max;
	y_scale = (float)h()/(2.0*y_max);

	// erase last sample

	fl_color(FL_BLACK);
	x1 = x_scale * index + x();
	y1 = y_scale * mem[index];
	y1 = clip(y1);
	y1 = y() + h()/2 - y1;
	fl_point(x1, y1);

	// draw new sample

	fl_color(FL_GREEN);
	x1 = x_scale * index + x();
	y1 = y_scale * new_sample;
	y1 = clip(y1);
	y1 = y() + h()/2 - y1;
	fl_point(x1, y1);
	mem[index] = new_sample;

	index++;
	if (index >=  x_max)
	    index = 0;

	// y axis graticule

	step = 10;

	while ((2.0*y_max/step) > 10)
	    step *= 2.0;
	while ((2.0*y_max/step) < 4)
	    step /= 2.0;

	fl_color(FL_DARK_GREEN);
	fl_line_style(FL_DOT);
	for(i=-y_max; i<y_max; i+=step) {
	    x1 = x();
	    y1 = y() + h()/2 - i*y_scale;
	    x2 = x() + w();
	    y2 = y1;
	    fl_line(x1,y1,x2,y2);   
	}

	// y axis graticule labels

	fl_color(FL_GREEN);
	fl_line_style(FL_SOLID);
	for(i=-y_max; i<y_max; i+=step) {
	    x1 = x();
	    y1 = y() + h()/2 - i*y_scale;
	    sprintf(label, "%d", i);
	    fl_draw(label, x1, y1);
	}

	fl_pop_clip();
    }

public:
    Scalar(int x, int y, int w, int h, int x_max_, int y_max_, const char name[]): Fl_Box(x, y, w, h, name)
    {
	int i;

	align(FL_ALIGN_TOP);
	labelsize(10);

	x_max = x_max_; y_max = y_max_;
	mem = new float[x_max];
	for(i=0; i<x_max; i++) {
	    mem[i] = 0.0;
	}

	prev_w = 0; prev_h = 0; prev_x = 0; prev_y = 0;
	index = 0;
    };

    ~Scalar() {
	delete mem;
    }

    void add_new_sample(float sample) {
	new_sample = sample;
    }

};


// update average of each spectrum point
    
void new_data(float mag_dB[]) {
    int i;

    for(i=0; i<FDMDV_NSPEC; i++)
	av_mag[i] = (1.0 - BETA)*av_mag[i] + BETA*mag_dB[i];
}

/*
  The sample source could be a sound card or file.  The sample source
  supplies a fixed number of samples with each call.  However
  fdmdv_demod requires a variable number of samples for each call.
  This function will buffer as appropriate and call fdmdv_demod with
  the correct number of samples.

  The processing sequence is:

  collect demod input samples from sound card 1 A/D
  while we have enough samples:
    demod samples into bits
    decode bits into speech samples
    output a buffer of speech samples to sound card 2 D/A

  Note that sound card 1 and sound card 2 will have slightly different
  sample rates, as their sample clocks are not syncronised.  We
  effectively lock the system to the demod A/D (sound card 1) sample
  rate. This ensures the demod gets a continuous sequence of samples,
  maintaining sync. Sample underflow or overflow will instead occur on
  the sound card 2 D/A.  This is acceptable as a buffer of lost or
  extra speech samples is unlikely to be noticed.

  The situation is actually a little more complex than that.  Through
  the demod timing estimation the buffers supplied to sound card D/A 2
  are effectively clocked at the remote modulator sound card D/A clock
  rate.  We slip/gain buffers supplied to sound card 2 to compensate.

  idle() is the FLTK function that gets continusouly called when FLTK
  is not doing GUI work.  In a non-GUI program this could just be a
  while(1) loop.
*/

void idle(void*) {
    struct FDMDV_STATS stats;
    int   rx_bits[FDMDV_BITS_PER_FRAME];
    int   sync_bit;
    float rx_fdm[FDMDV_MAX_SAMPLES_PER_FRAME];
    float rx_spec[FDMDV_NSPEC];
    int   i, j, nin_prev, ret;

    if (fin_name != NULL) {
	ret = fread(&rx_fdm_scaled[nbuf], 
		    sizeof(short), 
		    FDMDV_NOM_SAMPLES_PER_FRAME, 
		    fin);

	// simulate time delay from real world A/D input

	usleep(20000);
    }

    if (sound_dev_name != NULL) {    
	err = Pa_ReadStream(stream, 
			   &rx_fdm_scaled[nbuf], 
			   FDMDV_NOM_SAMPLES_PER_FRAME);
	if (err & paInputOverflow) { 
	    fprintf( stderr, "Input Overflow.\n" );
	}
    }

    nbuf += FDMDV_NOM_SAMPLES_PER_FRAME;
    assert(nbuf <= (2*FDMDV_NOM_SAMPLES_PER_FRAME));

    // this will run the demod 0, 1 (nominal) or 2 time

    while(nbuf >= nin) {

	Ts += (float)nin/FS;
	
	// demod per frame processing

	for(i=0; i<nin; i++)
	    rx_fdm[i] = (float)rx_fdm_scaled[i]/FDMDV_SCALE;
	nin_prev = nin;
	//fdmdv_demod(fdmdv, rx_bits, &sync_bit, rx_fdm, &nin);
	nbuf -= nin_prev;
	assert(nbuf >= 0);

	// get spectrum, stats and update windows

	fdmdv_get_rx_spectrum(fdmdv, rx_spec, rx_fdm, nin_prev);
	fdmdv_get_demod_stats(fdmdv, &stats);
	new_data(rx_spec);
	aScatter->add_new_samples(stats.rx_symbols);
	aTimingEst->add_new_sample(stats.rx_timing);
	aFreqEst->add_new_sample(stats.foff);
	aSNR->add_new_sample(stats.snr_est);

	// update plots every DT

	if (Ts >= DT) {
	    Ts -= DT;
	    if (!zoomSpectrumWindow->shown() && !zoomWaterfallWindow->shown()) {
		aSpectrum->redraw();
		aWaterfall->redraw();
		aScatter->redraw();
		aTimingEst->redraw();
		aFreqEst->redraw();
		aSNR->redraw();
	    }
	    if (zoomSpectrumWindow->shown())		
		aZoomedSpectrum->redraw();		
	    if (zoomWaterfallWindow->shown())		
		aZoomedWaterfall->redraw();		
	}

	// shift buffer

	for(i=0; i<nbuf; i++)
	    rx_fdm_scaled[i] = rx_fdm_scaled[i+nin_prev];

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
    if (argv[i][1] == 's') {
	if ((i+1) >= argc) 
	    return 0;
	sound_dev_name = argv[i+1];
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

    if (argc == 1) {
	printf("usage: %s [-i inputFdmdvRawFile] [-s inputSoundDevice]\n", argv[0]);
	exit(0);
    }

    if (fin_name != NULL) {
	fin = fopen(fin_name,"rb");
	if (fin == NULL) {
	    fprintf(stderr, "Error opening input fdmdv raw file %s\n", fin_name);
	    exit(1);
	}
    }

    if (sound_dev_name != NULL) {
	err = Pa_Initialize();
	if( err != paNoError ) goto pa_error;
	inputParameters.device = Pa_GetDefaultInputDevice(); /* default input device */
	printf( "Input device # %d.\n", inputParameters.device );
	printf( "Input LL: %g s\n", Pa_GetDeviceInfo( inputParameters.device )->defaultLowInputLatency );
	printf( "Input HL: %g s\n", Pa_GetDeviceInfo( inputParameters.device )->defaultHighInputLatency );
	inputParameters.channelCount = 1;
	inputParameters.sampleFormat = paInt16;
	inputParameters.suggestedLatency = Pa_GetDeviceInfo( inputParameters.device )->defaultHighInputLatency ;
	inputParameters.hostApiSpecificStreamInfo = NULL;

	err = Pa_OpenStream(
              &stream,
              &inputParameters,
              NULL,                         // no output
              /*SOUND_CARD_FS*/8000,                
              FDMDV_NOM_SAMPLES_PER_FRAME,
              paClipOff,      
              NULL,                         // no callback, use blocking API
              NULL );                       // no callback, so no callback userData
	if( err != paNoError ) goto pa_error;
	err = Pa_StartStream( stream );
	if( err != paNoError ) goto pa_error;
    }
    for(i=0; i<FDMDV_NSPEC; i++)
	av_mag[i] = -40.0;

    // recommended to prevent dithering and stopped display being
    // covered by black flickering squares

    Fl::visual(FL_RGB);

    // set up main window

    window = new Fl_Window(W, SP+H2+SP+SP+H2+SP, "fl_fmdv");
    //window->size_range(100, 100);
    //window->resizable();
    aSpectrum = new Spectrum(SP, SP, W3-2*SP, H2);
    aWaterfall = new Waterfall(SP, SP+H2+SP+SP, W3-2*SP, H2);
    aScatter = new Scatter(W3+SP, SP, W3-2*SP, H2);
    aTimingEst = new Scalar(W3+SP, SP+H2+SP+SP, W3-2*SP, H2, 100, 80, "Timing Est");
    aFreqEst = new Scalar(2*W3+SP, SP, W3-2*SP, H2, 100, 100, "Frequency Est");
    aSNR = new Scalar(2*W3+SP, SP+H2+SP+SP, W3-2*SP, H2, 100, 10, "SNR");
    fdmdv = fdmdv_create();

    Fl::add_idle(idle);

    window->end();

    // set up zoomed spectrum window

    zoomSpectrumWindow = new Fl_Window(W, H, "Spectrum");
    aZoomedSpectrum = new Spectrum(SP, SP, W-2*SP, H-2*SP);
    zoomSpectrumWindow->end();

    // set up zoomed waterfall window

    zoomWaterfallWindow = new Fl_Window(W, H, "Waterfall");
    aZoomedWaterfall = new Waterfall(SP, SP, W-2*SP, H-2*SP);
    zoomWaterfallWindow->end();

    // show the main window and start running

    window->show(argc, argv);
    Fl::run();

    if (sound_dev_name != NULL) {
	err = Pa_StopStream( stream );
	if( err != paNoError ) goto pa_error;
	Pa_CloseStream( stream );
	Pa_Terminate();
    }

    fdmdv_destroy(fdmdv);
    if (fin_name != NULL)
	fclose(fin);

    return ret;

    // Portaudio error handling

pa_error:
    if( stream ) {
       Pa_AbortStream( stream );
       Pa_CloseStream( stream );
    }
    Pa_Terminate();
    fprintf( stderr, "An error occured while using the portaudio stream\n" );
    fprintf( stderr, "Error number: %d\n", err );
    fprintf( stderr, "Error message: %s\n", Pa_GetErrorText( err ) );
    return -1;
}
