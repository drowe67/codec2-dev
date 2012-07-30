/*
  genphdata.c

  Generates test phase data for trainvqph testing.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>

typedef struct {
    float real;
    float imag;
} COMP;

#define N 100000
#define K 2
#define M 8
#define PI         3.141592654	/* mathematical constant                */
#define TWO_PI     6.283185307	/* mathematical constant                */

int main(void) {
    FILE *f=fopen("testph.flt", "wb");
    int   i;
    float angle;
    COMP  c;

    #ifdef TEST1
    for(i=0; i<M*K; i++) {
	c.real = cos(i*TWO_PI/(M*K));
	c.imag = sin(i*TWO_PI/(M*K));
	fwrite(&c, sizeof(COMP), 1, f);
    }
    #endif

    for(i=0; i<N; i++) {
	angle = PI*(1.0 - 2.0*rand()/RAND_MAX);
	c.real = cos(angle);
	c.imag = sin(angle);
	fwrite(&c, sizeof(COMP), 1, f);
    }

    return 0;
}
