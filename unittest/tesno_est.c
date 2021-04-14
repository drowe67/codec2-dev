/*---------------------------------------------------------------------------*\

  FILE........: tesno_est.c
  AUTHORS.....: David Rowe
  DATE CREATED: Mar 2021

  Test for C port of Es/No estimator.
  
\*---------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ofdm_internal.h"

int main(int argc, char *argv[])
{
    FILE *fin = fopen(argv[1],"rb"); assert(fin != NULL);
    size_t nsym = atoi(argv[2]); assert(nsym >= 0);
    complex float rx_sym[nsym];
    size_t nread = fread(rx_sym, sizeof(complex float), nsym, fin);
    assert(nread == nsym);
    fclose(fin);
    
    float sig_var, noise_var;
    esno_est_calc(&sig_var, &noise_var, rx_sym, nsym);
    printf("%f\n",10.0*log10(sig_var/noise_var));

    return 0;
}
