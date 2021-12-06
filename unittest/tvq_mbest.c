/*
  tvq_mbest.c
  David Rowe Dec 2019

  Generate some test vectors to exercise misc/vq_mbest.c
*/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

void write_float_file(char fn[], float *values, int n) {
    FILE *f=fopen(fn,"wb");
    assert(f != NULL);
    assert(fwrite(values, sizeof(float), n, f) == n);
    fclose(f);
}

int main(void) {
    float target[] = {1.0,1.0};
    write_float_file("target.f32", target, 2);
    float vq1[] = {0.9,0.9,  /* this will be a better match on first stage */
		   0.8,0.8}; /* but after second stage should choose this  */
    write_float_file("vq1.f32", vq1, 4);
    float vq2[] = {0.3,0.3,
		   0.2,0.2}; /* 0.8+0.2 == 1.0 so best 2nd stage entry     */
    write_float_file("vq2.f32", vq2, 4);
    return 0;
}
