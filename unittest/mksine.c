/*
  mksine.c
  David Rowe 
  10 Nov 2010

  Creates a file of sine wave samples.
*/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define TWO_PI     6.283185307
#define FS         8000.0
#define AMP        10000.0

int main(int argc, char *argv[]) {
    FILE *f;
    int   i,n;
    float freq, length;
    short *buf;

    if (argc != 4) {
	printf("usage: %s outputFile frequencyHz lengthSecs\n", argv[0]);
	exit(1);
    }

    f = fopen(argv[1] ,"wb");
    freq = atof(argv[2]);
    length = atof(argv[3]);
    n = length*FS;
    buf = (short*)malloc(sizeof(short)*n);
    assert(buf != NULL);

    for(i=0; i<n; i++)
	buf[i] = AMP*cos(freq*i*(TWO_PI/FS));

    fwrite(buf, sizeof(short), n, f);

    return 0;
}
