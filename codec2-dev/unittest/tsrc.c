/*
   tsrc.c
   David Rowe
   Sat Nov 3 2012

   Unit test for libresample code.

   build: gcc tsrc.c -o tsrc -lm -lsamplerate

  */

#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <samplerate.h>

#define N    10000                   /* processing buffer size */

int main(int argc, char *argv[]) {
    FILE       *fin, *fout;
    short       in_short[N], out_short[N];
    float       in[N], out[N];
    SRC_STATE  *src;
    SRC_DATA    data;
    int         error, nin, nremaining, i;

    if (argc != 4) {
	printf("usage %s inputRawFile OutputRawFile OutSampleRatio\n", argv[0]);
	exit(0);
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

    src = src_new(SRC_SINC_FASTEST, 1, &error);
    //src = src_new(SRC_LINEAR, 1, &error);
    assert(src != NULL);

    data.data_in = in;
    data.data_out = out;
    data.input_frames = N;
    data.output_frames = N;
    data.end_of_input = 0;
    data.src_ratio = atof(argv[3]);

    int total_in = 0;
    int total_out = 0;

    nin = N;
    nremaining = 0;
    while(fread(&in_short[nremaining], sizeof(short), nin, fin) == nin) {
	src_short_to_float_array(in_short, in, N);
	error = src_process(src, &data);
        assert(error == 0);
	src_float_to_short_array(out, out_short, data.output_frames_gen);

	fwrite(out_short, sizeof(short), data.output_frames_gen, fout);
        if (fout == stdout) fflush(stdout);

        nremaining = N - data.input_frames_used;
        nin = data.input_frames_used;
	//fprintf(stderr, "input frames: %d output_frames %d nremaining: %d\n", 
        //        (int)data.input_frames_used, (int)data.output_frames_gen, nremaining);
        for(i=0; i<nremaining; i++)
            in_short[i] = in_short[i+nin];

        total_in  += data.input_frames_used;
        total_out += data.output_frames_gen;
    }

    //fprintf(stderr, "total_in: %d total_out: %d\n", total_in, total_out);

    fclose(fout);
    fclose(fin);

    return 0;
}
