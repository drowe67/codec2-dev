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
#include <getopt.h>

#define FS 8000

int main(int argc, char *argv[]) {
    short buf[FS] = {0};
    float f0 = 60.0;
    float n0 = 0.0;
    int   Nsecs = 1;
    int   randf0 = 0;
    int   filter = 0;
    
    int o = 0;
    int opt_idx = 0;
    while( o != -1 ) {
        static struct option long_opts[] = {
            {"help",   no_argument,       0, 'h'},
            {"n0",     required_argument, 0, 'n'},
            {"f0",     required_argument, 0, 'f'},
            {"secs",   required_argument, 0, 's'},
            {"randf0", no_argument, 0, 'r'},
            {"filter", no_argument, 0, 'i'},
            {0, 0, 0, 0}
        };
        
        o = getopt_long(argc,argv,"hn:f:s:r",long_opts,&opt_idx);
        
        switch(o) {
        case 'n':
	    n0 = atof(optarg);
            break;
        case 'f':
            f0 = atof(optarg);
	    break;
        case 's':
            Nsecs = atoi(optarg);
	    break;
        case 'r':
            randf0 = 1;
	    break;
        case 'i':
            filter = 1;
	    break;
        case '?':
        case 'h':
	    fprintf(stderr, "usage: %s [--f0 f0Hz] [--n0 samples] [--secs Nsecs]\n"
	                    "[--randf0]  choose a random F0 every second\n\n", argv[0]);
	    exit(1);      
	break;
        }
    }

    int t = 0;
    
    /* optionally filter with 2nd order system */
    float alpha = 0.25*M_PI, gamma=0.99;
    float a[2] = {-2.0*gamma*cos(alpha), gamma*gamma};
    float mem[2] = {0};
    
    for (int j=0; j<Nsecs; j++) {
	if (randf0) {
	    float pitch_period = FS/400.0 + (FS/50.0 - FS/400.0)*rand()/RAND_MAX;
	    f0 = (float)FS/pitch_period;
	    fprintf(stderr, "P: %f f0: %f\n", pitch_period, f0);
	}
	float Wo = 2.0*M_PI*f0/FS;
	int L = M_PI/Wo;
	for(int i=0; i<FS; i++) {
	    buf[i] = 0;
	    for(int m=1; m<L; m++)
		buf[i] += (1000.0/L)*cos(m*Wo*(t + n0));
	    t++;
	}
	if (filter) {
	    for(int i=0; i<FS; i++) {
		float x = (float)buf[i];
		float y = (x - mem[0]*a[0] - mem[1]*a[1]);
		mem[1] = mem[0]; mem[0] = y;
		buf[i] = (short)y;
	    }
	}
    
	fwrite(buf, sizeof(short), FS, stdout);
    }
    
    return 0;
}
