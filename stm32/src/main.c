#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "stm32f4xx_conf.h"
#include "stm32f4xx.h"
#include "gdb_stdio.h"
#include "codec2.h"
#include "dump.h"
#include "sine.h"
#include "machdep.h"

#ifdef __EMBEDDED__
#define printf gdb_stdio_printf
#define fopen gdb_stdio_fopen
#define fclose gdb_stdio_fclose
#define fread gdb_stdio_fread
#define fwrite gdb_stdio_fwrite
#endif

static void c2demo(int mode, char inputfile[], char outputfile[])
{
    struct CODEC2 *codec2;
    short         *inbuf, *outbuf;
    unsigned char *bits;
    int            nsam, nbit;
    FILE          *fin, *fout;
    int            frame;
    TIMER_VAR(enc_start, dec_start);

    codec2 = codec2_create(mode);
    nsam = codec2_samples_per_frame(codec2);
    outbuf = (short*)malloc(nsam*sizeof(short));
    inbuf = (short*)malloc(nsam*sizeof(short));
    nbit = codec2_bits_per_frame(codec2);
    bits = (unsigned char*)malloc(nbit*sizeof(char));

    fin = fopen(inputfile, "rb");
    if (fin == NULL) {
        printf("Error opening input file: %s\n\nTerminating....\n",inputfile);
        exit(1);
    }

    fout = fopen(outputfile, "wb");
    if (fout == NULL) {
        printf("Error opening output file: %s\n\nTerminating....\n",outputfile);
        exit(1);
    }

    #ifdef DUMP
    dump_on("stm32f4");
    #endif
    frame = 0;

    while (fread(inbuf, sizeof(short), nsam, fin) == nsam) {
        TIMER_SAMPLE(enc_start);
        codec2_encode(codec2, bits, inbuf);
        TIMER_SAMPLE_AND_LOG(dec_start, enc_start, "  enc");     
	codec2_decode(codec2, outbuf, bits);
        TIMER_SAMPLE_AND_LOG2(dec_start, "  dec");     
        TIMER_SAMPLE_AND_LOG2(enc_start, "  enc & dec");     
        fwrite((char*)outbuf, sizeof(short), nsam, fout);
        printf("frame: %d\n", ++frame);
        machdep_timer_print_logged_samples();
    }

    #ifdef DUMP
    dump_off("sm32f4");
    #endif

    fclose(fin);
    fclose(fout);
    free(inbuf);
    free(outbuf);
    free(bits);
    codec2_destroy(codec2);
}

#define SPEED_TEST_SAMPLES 24000

static void c2speedtest(int mode, char inputfile[])
{
    struct CODEC2 *codec2;
    short         *inbuf, *outbuf, *pinbuf;
    unsigned char *bits;
    int            nsam, nbit, nframes;
    FILE          *fin;
    int            f, nread;

    codec2 = codec2_create(mode);
    nsam = codec2_samples_per_frame(codec2);
    nframes = SPEED_TEST_SAMPLES/nsam;
    outbuf = (short*)malloc(nsam*sizeof(short));
    inbuf = (short*)malloc(SPEED_TEST_SAMPLES*sizeof(short));
    nbit = codec2_bits_per_frame(codec2);
    bits = (unsigned char*)malloc(nbit*sizeof(char));

    fin = fopen(inputfile, "rb");
    if (fin == NULL) {
        printf("Error opening input file: %s\nTerminating....\n",inputfile);
        exit(1);
    }

    nread = fread(inbuf, sizeof(short), SPEED_TEST_SAMPLES, fin);
    if (nread != SPEED_TEST_SAMPLES) {
        printf("error reading %s, %d samples reqd, %d read\n", 
               inputfile, SPEED_TEST_SAMPLES, nread);
    }
    fclose(fin);
    
    pinbuf = inbuf;
    for(f=0; f<nframes; f++) {
	GPIOD->ODR = (1 << 13);
        codec2_encode(codec2, bits, pinbuf);
        pinbuf += nsam;
	GPIOD->ODR &= ~(1 << 13);
	codec2_decode(codec2, outbuf, bits);
    }

    free(inbuf);
    free(outbuf);
    free(bits);
    codec2_destroy(codec2);
}

void gpio_init() {
    RCC->AHB1ENR |= RCC_AHB1ENR_GPIODEN; // enable the clock to GPIOD 
    GPIOD->MODER = (1 << 26);            // set pin 13 to be general 
                                         // purpose output
}

int main(int argc, char *argv[]) {
    SystemInit();
    gpio_init();
    machdep_timer_init ();
 
    printf("Starting c2demo\n");

    /* File I/O test for profiling or (with #define DUMP)
       dumping states for optimisation and tiuning */

    c2demo(CODEC2_MODE_1600, "stm_in.raw", "stm_out.raw");

    printf("Starting c2 speed test\n");
    
    /* Another test of execution speed. Look at PD13 with a
       oscilliscope.  On time is enc, off is dec */

    c2speedtest(CODEC2_MODE_1600, "stm_in.raw");

    printf("Finished\n");

    return 0;
}

