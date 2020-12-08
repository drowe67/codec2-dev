/*---------------------------------------------------------------------------*\

  FILE........: cohpsk_ch.c
  AUTHOR......: David Rowe
  DATE CREATED: May 2015

  Channel impairment program for testing command line versions of
  cohpsk (and other) modems.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2015 David Rowe

  All rights reserved.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License version 2.1, as
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

#include "codec2_fdmdv.h"
#include "codec2_cohpsk.h"
#include "comp_prim.h"
#include "noise_samples.h"
#include "ht_coeff.h"
#include "ssbfilt_coeff.h"

#include "debug_alloc.h"

#define BUF_N                 160
#define FAST_FADING_DELAY_MS  2.0
#define SLOW_FADING_DELAY_MS  0.5
#define PAPR_TARGET           7.0

/*
  Use Octave to generate the fading channel models:

     octave:24> pkg load signal
     octave:25> cohpsk_ch_fading("../raw/faster_fading_samples.float", 8000, 2.0, 8000*60)
     octave:26> cohpsk_ch_fading("../raw/fast_fading_samples.float", 8000, 1.0, 8000*60)
     octave:27> cohpsk_ch_fading("../raw/slow_fading_samples.float", 8000, 0.1, 8000*60)

   Note: for Fs=8000Hz operation 7500 Hz is OK - these are just the two path fading complex numbers,
   a few % different in fading bandwidth won't matter.
*/

#define DEFAULT_RAW_DIR "../../raw"
#define SLOW_FADING_FILE_NAME   "slow_fading_samples.float"
#define FAST_FADING_FILE_NAME   "fast_fading_samples.float"
#define FASTER_FADING_FILE_NAME "faster_fading_samples.float"

int opt_exists(char *argv[], int argc, char opt[]) {
    int i;
    for (i=0; i<argc; i++) {
        if (strcmp(argv[i], opt) == 0) {
            return i;
        }
    }
    return 0;
}

int main(int argc, char *argv[])
{
    FILE          *fin, *ffading, *fout;
    char	   *raw_dir;
    float          NodB, foff_hz;
    int            fading_en, nhfdelay;

    short          buf[BUF_N];
    float          htbuf[HT_N+BUF_N];
    COMP           ch_in[BUF_N];
    COMP           ch_fdm[BUF_N];
    COMP           ssbfiltbuf[SSBFILT_N+BUF_N];
    COMP           ssbfiltout[BUF_N];

    COMP           phase_ch;
    int            noise_r, noise_end;
    float          No, variance;
    COMP           scaled_noise;
    float          hf_gain;
    COMP          *ch_fdm_delay = NULL, aspread, aspread_2ms, delayed, direct;
    float          tx_pwr, tx_pwr_fade, noise_pwr;
    int            frames, i, j, k, Fs, ret, clipped, ssbfilt_en, complex_out;
    float          sam, peak, inclip, papr, CNo, snr3k;

    if (argc > 3) {
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

        NodB = atof(argv[3]);
        Fs = COHPSK_FS; foff_hz = 0.0; fading_en = 0; inclip = 1.0; ssbfilt_en = 1; complex_out = 0;
        raw_dir = strdup(DEFAULT_RAW_DIR);

        for(int i=4; i<argc; i++) {
            if (!strcmp(argv[i],"--Fs")) { Fs = atoi(argv[i+1]); i++; }
            else if (!strcmp(argv[i], "-f")) { foff_hz = atoi(argv[i+1]); i++; }
            else if (!strcmp(argv[i], "--slow")) fading_en = 1;
            else if (!strcmp(argv[i], "--fast")) fading_en = 2;
            else if (!strcmp(argv[i], "--faster")) fading_en = 3;
            else if (!strcmp(argv[i], "--clip")) { inclip = atof(argv[i+1]); i++; }
            else if (!strcmp(argv[i], "--ssbfilt")) { ssbfilt_en = atof(argv[i+1]); i++; }
            else if (!strcmp(argv[i], "--complexout")) complex_out = 1;
            else if (!strcmp(argv[i], "--raw_dir")) {
                FREE(raw_dir); raw_dir = strdup(argv[i+1]); i++;
            } else {
                fprintf(stderr, "Unknown argument: %s\n", argv[i]);
                exit(1);
            }
        }
    }
    else {
        fprintf(stderr, "usage: %s InputRealModemRawFile OutputRealModemRawFile No(dB/Hz) [--Fs SampleRateHz]"
                        " [-f FoffHz] [--slow] [--fast] [--faster] [--clip 0to1] [--ssbfilt 0|1] [--raw_dir Path] [--complexout]\n", argv[0]);
        exit(1);
    }
    fprintf(stderr, "cohpsk_ch ----------------------------------------------------------------------------------\n");
    fprintf(stderr, "Fs: %d NodB: %4.2f foff: %4.2f Hz fading: %d inclip: %4.2f ssbfilt: %d complexout: %d\n",
            Fs, NodB, foff_hz, fading_en, inclip, ssbfilt_en, complex_out);
    fprintf(stderr, "cohpsk_ch ----------------------------------------------------------------------------------\n");

    phase_ch.real = 1.0; phase_ch.imag = 0.0;
    noise_r = 0;
    noise_end = sizeof(noise)/sizeof(COMP);

    /*  N = var = NoFs */

    No = pow(10.0, NodB/10.0);
    variance = Fs*No;

    tx_pwr = tx_pwr_fade = noise_pwr = 0.0;
    clipped = 0;
    peak = 0.0;

    /* init HF fading model */

    ffading = NULL;
    nhfdelay = 0;
    if (fading_en) {
        char fname[256];

        if (fading_en == 1) {
	          sprintf(fname, "%s/%s", raw_dir, SLOW_FADING_FILE_NAME);
            ffading = fopen(fname, "rb");
            if (ffading == NULL) {
            cant_load_fading_file:
                fprintf(stderr, "-----------------------------------------------------\n");
                fprintf(stderr, "cohpsk_ch ERROR: Can't find fading file: %s\n", fname);
                fprintf(stderr, "\nAdjust path --raw_dir or use GNU Octave to generate:\n\n");
                fprintf(stderr, "$ octave --no-gui\n");
                fprintf(stderr, "octave:24> pkg load signal\n");
                fprintf(stderr, "octave:25> cohpsk_ch_fading(\"../raw/faster_fading_samples.float\", 8000, 2.0, 8000*60)\n");
                fprintf(stderr, "octave:26> cohpsk_ch_fading(\"../raw/fast_fading_samples.float\", 8000, 1.0, 8000*60)\n");
                fprintf(stderr, "octave:27> cohpsk_ch_fading(\"../raw/slow_fading_samples.float\", 8000, 0.1, 8000*60)\n");
                fprintf(stderr, "-----------------------------------------------------\n");
                exit(1);
            }
            nhfdelay = floor(SLOW_FADING_DELAY_MS*Fs/1000);
        }

        if (fading_en == 2) {
	          sprintf(fname, "%s/%s", raw_dir, FAST_FADING_FILE_NAME);
            ffading = fopen(fname, "rb");
            if (ffading == NULL) goto cant_load_fading_file;
            nhfdelay = floor(FAST_FADING_DELAY_MS*Fs/1000);
        }

        if (fading_en == 3) {
	    sprintf(fname, "%s/%s", raw_dir, FASTER_FADING_FILE_NAME);
            ffading = fopen(fname, "rb");
            if (ffading == NULL) goto cant_load_fading_file;
            nhfdelay = floor(FAST_FADING_DELAY_MS*Fs/1000);
        }

        ch_fdm_delay = (COMP*)MALLOC((nhfdelay+COHPSK_NOM_SAMPLES_PER_FRAME)*sizeof(COMP));
        assert(ch_fdm_delay != NULL);
        for(i=0; i<nhfdelay+COHPSK_NOM_SAMPLES_PER_FRAME; i++) {
            ch_fdm_delay[i].real = 0.0;
            ch_fdm_delay[i].imag = 0.0;
        }

        /* first values in file are HF gains */

        for (i=0; i<4; i++)
            ret = fread(&hf_gain, sizeof(float), 1, ffading);
        //fprintf(stderr, "hf_gain: %f\n", hf_gain);
    }

    for(i=0; i<HT_N; i++) {
        htbuf[i] = 0.0;
    }
    for(i=0; i<SSBFILT_N; i++) {
        ssbfiltbuf[i].real = 0.0; ssbfiltbuf[i].imag = 0.0;
    }

    /* --------------------------------------------------------*\
	                          Main Loop
    \*---------------------------------------------------------*/
    frames = 0;
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

            sam = (float)buf[i];
            //printf("sam: %f ", sam);
            if (sam > inclip*32767.0)
                sam = inclip*32767.0;
            if (sam < -inclip*32767.0)
                sam = -inclip*32767.0;
            //printf("sam: %f\n", sam);
            htbuf[j] = sam/FDMDV_SCALE;

            if (fabs(htbuf[j]) > peak) {
                peak = fabs(htbuf[j]);
            }
            tx_pwr += pow(htbuf[j], 2.0);

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

	/* --------------------------------------------------------*\
	                          Channel
	\*---------------------------------------------------------*/

        fdmdv_freq_shift_coh(ch_fdm, ch_in, foff_hz, Fs, &phase_ch, BUF_N);

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
                if (ret == 0) {
                    fprintf(stderr, "Fading file finished - simulation stopping.  You may need more samples, see "
                            "top of src/cohpsk_ch.c for instructions.\n");
                    goto finish;
                    }
                ret = fread(&aspread_2ms, sizeof(COMP), 1, ffading);
                if (ret == 0) {
                    fprintf(stderr, "Fading file finished - simulation stopping.  You may need more samples, see "
                            "top of src/cohpsk_ch.c for instructions.\n");
                goto finish;
                }
                //printf("%f %f %f %f\n", aspread.real, aspread.imag, aspread_2ms.real, aspread_2ms.imag);

                direct    = cmult(aspread, ch_fdm[i]);
                delayed   = cmult(aspread_2ms, ch_fdm_delay[i]);
                ch_fdm[i] = fcmult(hf_gain, cadd(direct, delayed));
            }
        }

        /* Measure power after fading model to make sure average pwr
           is the same as AWGN channels. We only output the real
           signal, which is half the power. */

        for(i=0; i<BUF_N; i++) {
            tx_pwr_fade += pow(ch_fdm[i].real, 2.0);
        }

        /* AWGN noise ------------------------------------------*/

        for(i=0; i<BUF_N; i++) {
            scaled_noise = fcmult(sqrt(variance), noise[noise_r]);
            ch_fdm[i] = cadd(ch_fdm[i], scaled_noise);
            noise_pwr += pow(scaled_noise.real, 2.0);
            noise_r++;
            if (noise_r > noise_end) {
                noise_r = 0;
                //fprintf(stderr, "  [%d] noise wrap\n", f);
            }
        }

        /* FIR filter to simulate (a rather flat) SSB filter.  Might
           be useful to have an option for a filter with a few dB
           ripple too, to screw up the modem. This is mainly so analog
           SSB sounds realistic. */

        for(i=0, j=SSBFILT_N; i<BUF_N; i++,j++) {
            ssbfiltbuf[j] = ch_fdm[i];
            if (ssbfilt_en) {
                ssbfiltout[i].real = 0.0; ssbfiltout[i].imag = 0.0;
                for(k=0; k<SSBFILT_N; k++) {
                    ssbfiltout[i].real += ssbfiltbuf[j-k].real*ssbfilt_coeff[k];
                    ssbfiltout[i].imag += ssbfiltbuf[j-k].imag*ssbfilt_coeff[k];
                }
            }
            else {
                ssbfiltout[i] = ch_fdm[i];
            }
        }

        /* update SSB filter memory */
        for(i=0; i<SSBFILT_N; i++)
           ssbfiltbuf[i] = ssbfiltbuf[i+BUF_N];

	/* scale and save to disk as shorts */
        int nout = (complex_out+1)*BUF_N;
        short bufout[nout], *pout=bufout;
	for(i=0; i<BUF_N; i++) {
            sam = FDMDV_SCALE * ssbfiltout[i].real;
            if (sam >  32767.0) { clipped++; sam = 32767.0; }
            if (sam < -32767.0) { clipped++; sam = -32767.0; }
	    *pout++ = sam;
            if (complex_out) {
                sam = FDMDV_SCALE * ssbfiltout[i].imag;
                if (sam >  32767.0) { clipped++; sam = 32767.0; }
                if (sam < -32767.0) { clipped++; sam = -32767.0; }
                *pout++ = sam;
            }
        }

 	fwrite(bufout, sizeof(short), nout, fout);

	/* if this is in a pipeline, we probably don't want the usual
	   buffering to occur */

        if (fout == stdout) fflush(stdout);
        if (fin == stdin) fflush(stdin);
    }

 finish:
    fclose(fin);
    fclose(fout);

    fprintf(stderr, "cohpsk_ch -----------------------------------------------------------------------------\n");
    /*
    fprintf(stderr, "peak pwr: %7.2f av input pwr.: %7.2f av fading pwr: %7.2f noise pwr....: %7.2f\n",
            peak*peak,
            tx_pwr/(frames*BUF_N),
            tx_pwr_fade/(frames*BUF_N),
            noise_pwr/(frames*BUF_N)
           );
    */
    papr = 10*log10(peak*peak/(tx_pwr/(frames*BUF_N)));
    CNo = 10*log10(tx_pwr/(noise_pwr/(Fs/2))); // single sided spectrum magic IDFK!
    snr3k = CNo - 10*log10(3000);
    fprintf(stderr, "SNR3k(dB): %5.2f C/No: %4.1f PAPR: %4.1f \n", snr3k, CNo, papr);
    fprintf(stderr, "cohpsk_ch -----------------------------------------------------------------------------\n");

    if (ffading != NULL) fclose(ffading);
    if (ch_fdm_delay != NULL) FREE(ch_fdm_delay);
    return 0;
}
