/*
  est_n0.c
  David Rowe Dec 2019

  Estimate the position of n0, the impulse that excites each frame in
  the source/filter model.  Uses a file of Codec 2 model parameters on
  stdin as input.
*/

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "defines.h"

int main(void) {
    MODEL model;

    while(fread(&model, sizeof(MODEL), 1, stdin)) {
	float Wo = model.Wo; int L = model.L;
	float P = 2.0*M_PI/Wo;
	float best_error = 1E32; float best_n0=0.0;
	for(float n0=0; n0<=P; n0+=0.25) {
	    float error = 0.0;
	    for(int m=1; m<=L; m++) {
		complex diff = cexp(I*model.phi[m]) - cexp(I*n0*m*Wo);
		error += model.A[m]*model.A[m]*cabs(diff)*cabs(diff);
	    }
	    if (error < best_error) {
		best_error = error;
		best_n0 = n0;
	    }
	}
	printf("%f\n", best_n0);
    }
    return 0;
}
