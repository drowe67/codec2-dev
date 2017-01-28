/*
   tlininterp.c
   David Rowe
   Jan 2017

   Linear interpolator, CPU effecient way of getting large oversample
   ratios, if the input is already limited to a fraction of the
   input sampling rate.

   build: gcc tlininterp.c -o tlininterp -Wall -O2

*/

#include <assert.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#define N    10000                   /* processing buffer size */

void display_help(void) {
    fprintf(stderr, "\nusage: tlininterp inputRawFile OutputRawFile OverSampleRatio [-c]\n");
    fprintf(stderr, "\nUse - for stdin/stdout\n\n");
    fprintf(stderr, "-c complex (two channel) resampling\n\n");
}

int main(int argc, char *argv[]) {
    FILE       *fin, *fout;
    short       left, right, out;
    float       oversample, t;

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

    oversample = atof(argv[3]);

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

    left = 0;
    t = 0.0;
    while(fread(&right, sizeof(short), 1, fin) == 1) {
        while (t < 1.0) {
            out = (1.0 - t)*left + t*right;
            fwrite(&out, sizeof(short), 1, fout);
            t += 1.0/oversample;
        }
        t -= 1.0;
        left = right;
    }

    fclose(fout);
    fclose(fin);

    return 0;
}
