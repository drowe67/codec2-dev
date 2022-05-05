/*
  hpf.c
  David Rowe
  May 2022

  Streaming high pass filter tool.
*/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "lpc.h"
#include "bpf.h"
#include "bpfb.h"

#define N_SAMP 80

int main(void) {
    short buf[N_SAMP];
    float buf_float[N_SAMP];
    float bpf_buf[BPF_N+N_SAMP];
    int   i;
    
    for(i=0; i<BPF_N; i++)
        bpf_buf[i] = 0.0;

    while(fread(buf, sizeof(short), N_SAMP, stdin)) {
        for(i=0; i<N_SAMP; i++)
            bpf_buf[BPF_N+i] = buf[i];
        inverse_filter(&bpf_buf[BPF_N], hpf200, N_SAMP, buf_float, HPF200_N);
        for(i=0; i<N_SAMP; i++)
            buf[i] = buf_float[i];
        int nwrite = fwrite(buf,sizeof(short), N_SAMP, stdout);
        assert(nwrite == N_SAMP);
        for(i=0; i<BPF_N; i++)
            bpf_buf[i] = bpf_buf[N_SAMP+i];
    }

    return 0;
}
