/*
   tdec.c
   David Rowe
   Jan 2017

   Trivial non filtered decimator for high ration sample rate conversion.

   build: gcc tdec.c -o tdec -Wall -O2

*/

#include <assert.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#define SIGNED_16BIT   0
#define SIGNED_8BIT    1

void display_help(void) {
    fprintf(stderr, "\nusage: tdec inputRawFile OutputRawFile DecimationRatio [-c]\n");
    fprintf(stderr, "\nUse - for stdin/stdout\n\n");
    fprintf(stderr, "-c complex signed 16 bit input and output\n");
    fprintf(stderr, "-d complex signed 8 bit input, complex signed 16 bit output\n\n");
}

int main(int argc, char *argv[]) {
    FILE       *fin, *fout;
    short       dec;

    if (argc < 3) {
	display_help();
	exit(1);
    }

    if (strcmp(argv[1], "-") == 0) 
        fin = stdin;
    else
        fin = fopen(argv[1], "rb");
    assert(fin != NULL);

    if (strcmp(argv[2], "-") == 0) 
        fout = stdout;
    else
        fout = fopen(argv[2], "wb");
    assert(fout != NULL);

    dec = atoi(argv[3]);

    int channels = 1;
    int opt;
    int format = SIGNED_16BIT;
    while ((opt = getopt(argc, argv, "cd")) != -1) {
        switch (opt) {
        case 'c': channels = 2; break;
        case 'd': channels = 2; format = SIGNED_8BIT; break;
        default:
            display_help();
            exit(1);
        }
    }

    if (format == SIGNED_16BIT) {
        short buf[dec*channels];
        while(fread(buf, sizeof(short)*channels, dec, fin) == dec) {
            fwrite(buf, sizeof(short), channels, fout);
        }
    }
    else {
        uint8_t inbuf[dec*channels];
        short   outbuf[channels];
        short   sam, i;
        
        while(fread(inbuf, sizeof(uint8_t)*channels, dec, fin) == dec) {
            for (i=0; i<channels; i++) {
                sam = (short)inbuf[i];
                sam <<= 8;
                outbuf[i] = sam;
            }
            fwrite(outbuf, sizeof(short), channels, fout);
        }

    }

    fclose(fout);
    fclose(fin);

    return 0;
}
