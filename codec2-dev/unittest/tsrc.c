/*
   tsrc.c
   David Rowe
   Sat Nov 3 2012

   Unit test for libresample code.

   build: gcc tsrc.c -o tsrc -lm -lsamplerate

  */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <samplerate.h>

#define N8   160                     /* processing buffer size at 8 kHz       */
#define FS   8000                    /* nominal input sample rate             */
#define N48  ((int)N8*(48000/8000))  /* buf size assuming 48k max sample rate */

int main(int argc, char *argv[]) {
    FILE       *f8k, *fout;
    short       in8k_short[N8];
    float       in8k[N8];
    float       out[N48];
    short       out_short[N48];
    SRC_STATE  *src;
    SRC_DATA    data;
    int         error;

    if (argc != 4) {
	printf("usage %s inputRawFile OutputRawFile OutSampleRatio\n", argv[0]);
	exit(0);
    }

    if (strcmp(argv[1], "-") == 0) 
        f8k = stdin;
    else
        f8k = fopen(argv[1], "rb");
    assert(f8k != NULL);

    if (strcmp(argv[2], "-") == 0) 
        fout = stdout;
    else
        fout = fopen(argv[2], "wb");
    assert(fout != NULL);

    src = src_new(SRC_SINC_FASTEST, 1, &error);
    assert(src != NULL);

    data.data_in = in8k;
    data.data_out = out;
    data.input_frames = N8;
    data.output_frames = N48;
    data.end_of_input = 0;
    data.src_ratio = atof(argv[3]);

    while(fread(in8k_short, sizeof(short), N8, f8k) == N8) {
	src_short_to_float_array(in8k_short, in8k, N8);
	src_process(src, &data);
	//fprintf(stderr, "%d %d\n", (int)data.output_frames , (int)data.output_frames_gen);
	assert(data.output_frames_gen <= N48);
	src_float_to_short_array(out, out_short, data.output_frames_gen);
	fwrite(out_short, sizeof(short), data.output_frames_gen, fout);
        if (fout == stdout) fflush(stdout);
        if (f8k == stdin) fflush(stdin);
   }

    fclose(fout);
    fclose(f8k);

    return 0;
}
