/*
   tlininterp.c
   David Rowe
   Jan 2017

   Linear interpolator, CPU effecient way of getting large oversample
   ratios, if the input is already limited to a fraction of the
   input sampling rate.

   build: gcc tdec.c -o tdec -Wall -O2

*/

#include <assert.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

void display_help(void) {
    fprintf(stderr, "\nusage: tdec inputRawFile OutputRawFile DecimationRatio [-c]\n");
    fprintf(stderr, "\nUse - for stdin/stdout\n\n");
    fprintf(stderr, "-c complex (two channel) resampling\n\n");
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
    while ((opt = getopt(argc, argv, "c")) != -1) {
        switch (opt) {
        case 'c': channels = 2; break;
        default:
            display_help();
            exit(1);
        }
    }

    short buf[dec*channels];
    fprintf(stderr, "dec: %d\n", dec);
    while(fread(buf, sizeof(short)*channels, dec, fin) == dec) {
        fwrite(buf, sizeof(short), channels, fout);
    }

    fclose(fout);
    fclose(fin);

    return 0;
}
