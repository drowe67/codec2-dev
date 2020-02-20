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

int main(int argc, char *argv[]) {
    short buf[FS] = {0};
    float f0, n0;
    
    if (argc < 2) {
	printf("usage: %s f0 n0 [filter]\n", argv[0]);
	exit(1);
    }
    else {
	f0 = atof(argv[1]);
	n0 = atof(argv[2]);
    }
    float Wo = 2.0*M_PI*f0/FS;
    int L = M_PI/Wo;
   
    for(int i=0; i<FS; i++) {
	for(int m=1; m<L; m++)
	    buf[i] += (1000.0/L)*cos(m*Wo*(i+n0));
    }

    if (argc == 4) {
        /* filter with optional 2nd order system */
	float alpha = 0.25*M_PI, gamma=0.95;
	float a[2] = {-2.0*gamma*cos(alpha), gamma*gamma};
	float mem[2] = {0};
	for(int i=0; i<FS; i++) {
	    float x = (float)buf[i];
	    float y = (x - mem[0]*a[0] - mem[1]*a[1]);
	    mem[1] = mem[0]; mem[0] = y;
	    //printf("x: %f y: %f\n", x,y);
	    buf[i] = (short)y;
	}
    }
    
    fwrite(buf, sizeof(short), FS, stdout);
    
    return 0;
}
