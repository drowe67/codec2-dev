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
    unsigned int   enc_start, dec_start;

    codec2 = codec2_create(mode);
    nsam = codec2_samples_per_frame(codec2);
    outbuf = (short*)malloc(nsam*sizeof(short));
    inbuf = (short*)malloc(nsam*sizeof(short));
    nbit = codec2_bits_per_frame(codec2);
    bits = (unsigned char*)malloc(nbit*sizeof(char));

    fin = fopen(inputfile, "rb");
    if (fin == NULL) {
        printf("Error opening input file: %s\n",inputfile);
        exit(1);
    }

    fout = fopen(outputfile, "wb");
    if (fout == NULL) {
        printf("Error opening output file: %s\n",outputfile);
        exit(1);
    }

    #ifdef DUMP
    dump_on("stm32f4");
    #endif
    frame = 0;

    while (fread(inbuf, sizeof(short), nsam, fin) == nsam) {
        enc_start = machdep_timer_sample();
	codec2_encode(codec2, bits, inbuf);
        dec_start = machdep_timer_sample_and_log(enc_start, "  enc");     
	codec2_decode(codec2, outbuf, bits);
        machdep_timer_sample_and_log(dec_start, "  dec");     
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

int main(void) {
    SystemInit();
    printf("Starting\n");
    machdep_timer_init ();

    c2demo(CODEC2_MODE_1600, "hts1a.raw", "hts1a_out.raw");
    printf("Finished\n");

    return 0;
}


/*
 * Dummy function to avoid compiler error
 */
void _init() {

}

