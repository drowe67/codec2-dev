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
#define F0 200.0

int main(int argc, char *argv[]) {
    float Wo = 2.0*M_PI*F0/FS;
    int L = M_PI/Wo;
    short buf[FS] = {0};
    float n0 = 0;
    
    if (argc != 2) {
	printf("usage: %s n0\n", argv[0]);
	exit(1);
    }
    else
	n0 = atof(argv[1]);
    
    for(int i=0; i<FS; i++) {
	for(int m=1; m<L; m++)
	    buf[i] += 1000*cos(m*Wo*(i+n0));
    }
    fwrite(buf, sizeof(short), FS, stdout);
    
    return 0;
}
