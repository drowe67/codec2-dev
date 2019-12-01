/*
  timpulse.c
  David Rowe Dec 2019

  Generate an impulse train from a sum of sinusoids.  Test program for
  phaseNN project.
*/

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define FS 8000

int main(void) {
    float Wo = M_PI*200.0/FS;
    int L = M_PI/Wo;
    short buf[FS] = {0};
    
    for(int i=0; i<FS; i++) {
	for(int m=1; m<L; m++)
	    buf[i] += 1000*cos(m*Wo*i);
    }
    fwrite(buf, sizeof(short), FS, stdout);
    
    return 0;
}
