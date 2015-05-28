/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: cohpsk_ch.c
  AUTHOR......: David Rowe  
  DATE CREATED: May 2015
                                                                             
  Channel impairment program for testing command line versions of
  cohpsk modem.
     
  TODO: 
    [ ] measure and prints pwrs to check, prints warning
    [ ] SNR in 3000Hz input
    [ ] example operation with sox for sample rate change
    [ ] way to calibrate for different input pwrs
    [ ] HT to do real->complex
        [ ] check no BER hit just through HT
        [ ] unit test HT
    [ ] clipping detect

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2015 David Rowe

  All rights reserved.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License version 2, as
  published by the Free Software Foundation.  This program is
  distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this program; if not, see <http://www.gnu.org/licenses/>.
*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include "codec2_cohpsk.h"
#include "comp_prim.h"
#include "../unittest/noise_samples.h"
#include "ht_coeff.h"
#include "codec2_fdmdv.h"

#define BUF_N        160
#define HF_DELAY_MS  2.0

/* This file gets generated using the function write_noise_file in tcohpsk.m.  You have to run
   tcohpsk first (any variant) to load the function into Octave, e.g.:

  octave:17> tcohpsk
  octave:18> write_noise_file("../raw/fading_samples.float", 7500, 7500*60)
*/

#define FADING_FILE_NAME "../../raw/fading_samples.float"

int main(int argc, char *argv[])
{
    FILE          *fin, *ffading, *fout;
    float          EsNodB, foff_hz;
    int            fading_en, nhfdelay;

    short          buf[BUF_N];
    float          htbuf[HT_N+BUF_N];
    COMP           ch_in[BUF_N];
    COMP           ch_fdm[BUF_N];
                                           
    COMP           phase_ch;
    int            noise_r, noise_end;
    float          EsNo, variance;
    COMP           scaled_noise;
    float          hf_gain;
    COMP          *ch_fdm_delay, aspread, aspread_2ms, delayed, direct;
    float          tx_pwr, rx_pwr, noise_pwr;
    int            frames, i, j, k, ret;
    float          sam;

    if (argc == 6) {
        if (strcmp(argv[1], "-")  == 0) fin = stdin;
        else if ( (fin = fopen(argv[1],"rb")) == NULL ) {
            fprintf(stderr, "Error opening input modem raw file: %s: %s.\n",
                    argv[1], strerror(errno));
            exit(1);
        }

        if (strcmp(argv[2], "-") == 0) fout = stdout;
        else if ( (fout = fopen(argv[2],"wb")) == NULL ) {
            fprintf(stderr, "Error opening output modem raw file: %s: %s.\n",
                    argv[2], strerror(errno));
            exit(1);
        }

        EsNodB = atof(argv[3]);
        foff_hz = atof(argv[4]);
        fading_en = atoi(argv[5]);
    }
    else {
        fprintf(stderr, "usage: %s InputRealModemRawFileFs7500Hz OutputRealModemRawFileFs7500Hz SNR3000Hz FoffHz FadingEn\n", argv[0]);
        exit(1);
    }
    fprintf(stderr, "EsNodB: %4.2f foff: %4.2f Hz fading: %d\n", EsNodB, foff_hz, fading_en);
    
    phase_ch.real = 1.0; phase_ch.imag = 0.0; 
    noise_r = 0; 
    noise_end = sizeof(noise)/sizeof(COMP);

    /*  each carrier has power = 2, total power 2Nc, total symbol rate
        NcRs, noise BW B=Fs Es/No = (C/Rs)/(N/B), N = var =
        2NcFs/NcRs(Es/No) = 2Fs/Rs(Es/No) */

    EsNo = pow(10.0, EsNodB/10.0);
    variance = 2.0*COHPSK_FS/(COHPSK_RS*EsNo);
    
    tx_pwr = rx_pwr = noise_pwr = 0.0;

    /* init HF fading model */

    if (fading_en) {
        ffading = fopen(FADING_FILE_NAME, "rb");
        if (ffading == NULL) {
            printf("Can't find fading file: %s\n", FADING_FILE_NAME);
            exit(1);
        }
        nhfdelay = floor(HF_DELAY_MS*COHPSK_FS/1000);
        ch_fdm_delay = (COMP*)malloc((nhfdelay+COHPSK_SAMPLES_PER_FRAME)*sizeof(COMP));
        assert(ch_fdm_delay != NULL);
        for(i=0; i<nhfdelay+COHPSK_SAMPLES_PER_FRAME; i++) {
            ch_fdm_delay[i].real = 0.0;
            ch_fdm_delay[i].imag = 0.0;
        }

        /* first values in file are HF gains */

        for (i=0; i<4; i++)
            ret = fread(&hf_gain, sizeof(float), 1, ffading);
        fprintf(stderr, "hf_gain: %f\n", hf_gain);
    }

    for(i=0; i<HT_N; i++) {
        htbuf[i] = 0.0;
    }

    /* --------------------------------------------------------*\
	                          Main Loop
    \*---------------------------------------------------------*/
    
    while(fread(buf, sizeof(short), BUF_N, fin) == BUF_N) {
	frames++;

        /* Hilbert Transform to produce complex signal so we can do
           single sided freq shifts.  Allows us to use real signal I/O
           which is handy */

        for(i=0, j=HT_N; i<BUF_N; i++,j++) {

            /* 
               Hilbert Transform to produce complex signal so we can do
               single sided freq shifts.  Essential filters out negative
               freqencies.
            */

            htbuf[j] = (float)buf[i]/FDMDV_SCALE;

            /* FIR filter with HT to get imag, just delay to get real */

            ch_in[i].real = 0.0;
            ch_in[i].imag = 0.0;
            for(k=0; k<HT_N; k++) {
                ch_in[i].real += htbuf[j-k]*ht_coeff[k].real;
                ch_in[i].imag += htbuf[j-k]*ht_coeff[k].imag;
            }
            //printf("%d %f %f\n", i, ch_in[i].real, ch_in[i].imag);
        }
        assert(j <= (BUF_N+HT_N));
        
        /* update HT memory */

        for(i=0; i<HT_N; i++)
           htbuf[i] = htbuf[i+BUF_N];

        for(i=0; i<BUF_N; i++) {
            //printf("%d %f %f\n", i, ch_in[i].real, ch_in[i].imag);
            tx_pwr += pow(ch_in[i].real, 2.0) + pow(ch_in[i].imag, 2.0);
        }

        /* +3dB factor is beacuse we ouput a real signal, this has half
           the power of the complex version but as the noise reflects
           across from the -ve side the same noise power. */

        for(i=0; i<BUF_N; i++) {
            ch_in[i] = fcmult(sqrt(2.0), ch_in[i]);
        }

	/* --------------------------------------------------------*\
	                          Channel
	\*---------------------------------------------------------*/

        fdmdv_freq_shift(ch_fdm, ch_in, foff_hz, &phase_ch, BUF_N);

        /* optional HF fading -------------------------------------*/

        if (fading_en) {

            /* update delayed signal buffer */

            for(i=0; i<nhfdelay; i++)
                ch_fdm_delay[i] = ch_fdm_delay[i+BUF_N];
            for(j=0; j<BUF_N; i++, j++)
                ch_fdm_delay[i] = ch_fdm[j];

            /* combine direct and delayed paths, both multiplied by
               "spreading" (doppler) functions */

            for(i=0; i<BUF_N; i++) {
                ret = fread(&aspread, sizeof(COMP), 1, ffading);
                assert(ret == 1);
                ret = fread(&aspread_2ms, sizeof(COMP), 1, ffading);
                assert(ret == 1);
                //printf("%f %f %f %f\n", aspread.real, aspread.imag, aspread_2ms.real, aspread_2ms.imag);
                
                direct    = cmult(aspread, ch_fdm[i]);
                delayed   = cmult(aspread_2ms, ch_fdm_delay[i]);
                ch_fdm[i] = fcmult(hf_gain, cadd(direct, delayed));
            }
        }

        /* we only output the real signal, which is half the power. */

        for(i=0; i<BUF_N; i++) {
            rx_pwr += pow(ch_fdm[i].real, 2.0);
        }

        /* AWGN noise ------------------------------------------*/

        for(i=0; i<BUF_N; i++) {
            scaled_noise = fcmult(sqrt(variance), noise[noise_r]);
            ch_fdm[i] = cadd(ch_fdm[i], scaled_noise);
            noise_pwr += pow(noise[noise_r].real, 2.0) + pow(noise[noise_r].imag, 2.0);
            noise_r++;
            if (noise_r > noise_end) {
                noise_r = 0;
                //fprintf(stderr, "  [%d] noise wrap\n", f);            
            }
               
        }

	/* scale and save to disk as shorts */

	for(i=0; i<BUF_N; i++) {
            sam = FDMDV_SCALE * ch_fdm[i].real;
            if (fabs(sam) > 32767.0)
                fprintf(stderr,"clipping: %f\n", sam);
	    buf[i] = sam;
        }

 	fwrite(buf, sizeof(short), BUF_N, fout);

	/* if this is in a pipeline, we probably don't want the usual
	   buffering to occur */

        if (fout == stdout) fflush(stdout);
        if (fin == stdin) fflush(stdin);         
    }

    fclose(fin);
    fclose(fout);

    fprintf(stderr, "tx var: %f noise var: %f rx var: %f\n", 
           tx_pwr/(frames*BUF_N), 
           noise_pwr/(frames*BUF_N),
           rx_pwr/(frames*BUF_N) 
           );

    return 0;
}

