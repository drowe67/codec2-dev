#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "codec2.h"
#include "codec2_internal.h"
#include "defines.h"

static const float expect_w[] = {
  0.004293, 0.004301, 0.004309, 0.004315,
  0.004320, 0.004323, 0.004326, 0.004328,
  0.004328, 0.004328, 0.004326, 0.004323,
  0.004320, 0.004315, 0.004309, 0.004301};


static const float expect_W[] = {
 -0.002176,  0.002195,  0.004429, -0.008645,
 -0.012196,  0.065359,  0.262390,  0.495616, 
  0.601647,  0.495616,  0.262390,  0.065359, 
 -0.012196, -0.008645,  0.004429,  0.002195};


int float_cmp(float a, float b) {
    if ( fabsf(a - b) < 1e-6f ) return 1;
    else return 0;
    }

int main(int argc, char *argv[]) {

    struct CODEC2 *codec2;
    int i, j;

    ////////
    codec2 = codec2_create(CODEC2_MODE_700C);

    j = (codec2->c2const.m_pitch / 2) - 8;
    for (i=0; i<16; i++) {
        printf("w[%d] = %f", j+i, 
                (double)codec2->w[j+i]);
        if (!float_cmp(codec2->w[j+i], expect_w[i])) {
            printf(" Error, expected %f", (double)expect_w[i]);
            }
        printf("\n");
        }

    printf("\n");

    j = (FFT_ENC / 2) - 8;
    for (i=0; i<16; i++) {
        printf("W[%d] = %f", j+i, 
                (double)codec2->W[j+i].real);
        if (!float_cmp(codec2->W[j+i].real, expect_W[i])) {
            printf(" Error, expected %f", (double)expect_W[i]);
            }
        printf("\n");
        }

    printf("\nEnd of Test\n");
    fclose(stdout);
    fclose(stderr);

    return(0);
}

/* vi:set ts=4 et sts=4: */
