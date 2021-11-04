/*
  est_n0.c
  David Rowe Dec 2019

  Estimate the position of n0, the impulse that excites each frame in
  the source/filter model.  This defines the linear component of the
  phase spectra.  Uses a file of Codec 2 model parameters on stdin as
  input.
*/

#include <assert.h>
#include <complex.h>
#include <getopt.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "defines.h"

int main(int argc, char *argv[]) {
    MODEL model;
    FILE *fdisp;
    
    int o = 0; int opt_idx = 0;
    int extract = 0, add = 0;
    while (o != -1) {
       static struct option long_opts[] = {
            {"remove", no_argument, 0, 'r'},
            {"add",    required_argument, 0, 'a'},
            {0, 0, 0, 0}
        };

        o = getopt_long(argc,argv,"hra:",long_opts,&opt_idx);

        switch (o) {
        case 'a':
	    add = 1; fdisp = fopen(optarg, "rb");
	    if (fdisp == NULL) {
		fprintf(stderr, "Error opening disp model file: %s\n", optarg);
		exit(1);
	    }
	    break;
        case 'r': extract = 1; break;
        case 'h':
            fprintf(stderr, "usage: %s [-r] [-a disp.model]\n", argv[0]);
	    fprintf(stderr, "-r remove linear phase term and ouput model records\n");
	    fprintf(stderr, "-a est linear phase from model file on stdin, then add disp.model\n");
            exit(1);
        }
    }

    int wr = 0;
    while(fread(&model, sizeof(MODEL), 1, stdin)) {
	float Wo = model.Wo; int L = model.L;
	float P = 2.0*M_PI/Wo;
	float best_error = 1E32; float best_n0=0.0;

	/* note weighting MSE by log10(Am) works much better than
           Am*Am, the latter tends to fit a linear phase model between
           the two highest amplitude harmonics */

	for(float n0=0; n0<=P; n0+=0.25) {
	    float error = 0.0;
	    for(int m=1; m<=L; m++) {
		complex diff = cexp(I*model.phi[m]) - cexp(I*n0*m*Wo);
		error += log10(model.A[m])*cabs(diff)*cabs(diff);
	    }
	    if (error < best_error) {
		best_error = error;
		best_n0 = n0;
	    }
	}
	if (extract) {
	    for(int m=1; m<=L; m++) {
		complex diff = cexp(I*model.phi[m])*cexp(-I*best_n0*m*Wo);		
		model.phi[m] = carg(diff);
	    }
	    assert(fwrite(&model, sizeof(MODEL), 1, stdout));
	    wr++;
	}
	else if (add) {
	    MODEL disp;
	    assert(fread(&disp, sizeof(MODEL), 1, fdisp) == 1);
	    for(int m=1; m<=L; m++) {
		complex combined = cexp(I*disp.phi[m])*cexp(I*best_n0*m*Wo);		
		model.phi[m] = carg(combined);
	    }
	    assert(fwrite(&model, sizeof(MODEL), 1, stdout));
	    wr++;
	}
	else
	    printf("%f\n", best_n0);
    }
    fprintf(stderr, "wr: %d\n", wr);
    return 0;
}
