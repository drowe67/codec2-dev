/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: fdmdv.c
  AUTHOR......: David Rowe
  DATE CREATED: April 14 2012
                                                                             
  Functions that implement the FDMDV modem.
                                                                             
\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2012 David Rowe

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

/*---------------------------------------------------------------------------*\
                                                                             
                               INCLUDES

\*---------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "fdmdv_internal.h"
#include "codec2_fdmdv.h"
#include "rn.h"
#include "rxdec_coeff.h"
#include "test_bits.h"
#include "pilot_coeff.h"
#include "kiss_fft.h"
#include "hanning.h"
#include "os.h"
#include "machdep.h"

static int sync_uw[] = {1,-1,1,-1,1,-1};

#ifdef __EMBEDDED__
#define printf gdb_stdio_printf
#endif

/*---------------------------------------------------------------------------*\
                                                                             
                               FUNCTIONS

\*---------------------------------------------------------------------------*/

static COMP cneg(COMP a)
{
    COMP res;

    res.real = -a.real;
    res.imag = -a.imag;

    return res;
}

static COMP cconj(COMP a)
{
    COMP res;

    res.real = a.real;
    res.imag = -a.imag;

    return res;
}

static COMP cmult(COMP a, COMP b)
{
    COMP res;

    res.real = a.real*b.real - a.imag*b.imag;
    res.imag = a.real*b.imag + a.imag*b.real;

    return res;
}

static COMP fcmult(float a, COMP b)
{
    COMP res;

    res.real = a*b.real;
    res.imag = a*b.imag;

    return res;
}

static COMP cadd(COMP a, COMP b)
{
    COMP res;

    res.real = a.real + b.real;
    res.imag = a.imag + b.imag;

    return res;
}

static float cabsolute(COMP a)
{
    return sqrtf(pow(a.real, 2.0) + pow(a.imag, 2.0));
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: fdmdv_create	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 16/4/2012 

  Create and initialise an instance of the modem.  Returns a pointer
  to the modem states or NULL on failure.  One set of states is
  sufficient for a full duplex modem.

\*---------------------------------------------------------------------------*/

struct FDMDV * fdmdv_create(int Nc)
{
    struct FDMDV *f;
    int           c, i, k;

    assert(NC == FDMDV_NC_MAX);  /* check public and private #defines match */
    assert(Nc <= NC);
    assert(FDMDV_NOM_SAMPLES_PER_FRAME == M);
    assert(FDMDV_MAX_SAMPLES_PER_FRAME == (M+M/P));

    f = (struct FDMDV*)malloc(sizeof(struct FDMDV));
    if (f == NULL)
	return NULL;
    
    f->Nc = Nc;

    f->ntest_bits = Nc*NB*4;
    f->current_test_bit = 0;
    f->rx_test_bits_mem = (int*)malloc(sizeof(int)*f->ntest_bits);
    assert(f->rx_test_bits_mem != NULL);
    for(i=0; i<f->ntest_bits; i++)
	f->rx_test_bits_mem[i] = 0;
    assert((sizeof(test_bits)/sizeof(int)) >= f->ntest_bits);

    f->old_qpsk_mapping = 0;

    f->tx_pilot_bit = 0;

    for(c=0; c<Nc+1; c++) {
	f->prev_tx_symbols[c].real = 1.0;
	f->prev_tx_symbols[c].imag = 0.0;
	f->prev_rx_symbols[c].real = 1.0;
	f->prev_rx_symbols[c].imag = 0.0;

	for(k=0; k<NSYM; k++) {
	    f->tx_filter_memory[c][k].real = 0.0;
	    f->tx_filter_memory[c][k].imag = 0.0;
	}

	/* Spread initial FDM carrier phase out as far as possible.
           This helped PAPR for a few dB.  We don't need to adjust rx
           phase as DQPSK takes care of that. */
	
	f->phase_tx[c].real = cosf(2.0*PI*c/(Nc+1));
 	f->phase_tx[c].imag = sinf(2.0*PI*c/(Nc+1));

	f->phase_rx[c].real = 1.0;
 	f->phase_rx[c].imag = 0.0;

	for(k=0; k<NT*P; k++) {
	    f->rx_filter_mem_timing[c][k].real = 0.0;
	    f->rx_filter_mem_timing[c][k].imag = 0.0;
	}
    }
    
    fdmdv_set_fsep(f, FSEP);
    f->freq[Nc].real = cosf(2.0*PI*0.0/FS);
    f->freq[Nc].imag = sinf(2.0*PI*0.0/FS);
    f->freq_pol[Nc]  = 2.0*PI*0.0/FS;
    
    f->fbb_rect.real     = cosf(2.0*PI*FDMDV_FCENTRE/FS);
    f->fbb_rect.imag     = sinf(2.0*PI*FDMDV_FCENTRE/FS);
    f->fbb_pol           = 2.0*PI*FDMDV_FCENTRE/FS;
    f->fbb_phase_tx.real = 1.0;
    f->fbb_phase_tx.imag = 0.0;
    f->fbb_phase_rx.real = 1.0;
    f->fbb_phase_rx.imag = 0.0;

    /* Generate DBPSK pilot Look Up Table (LUT) */

    generate_pilot_lut(f->pilot_lut, &f->freq[Nc]);

    /* freq Offset estimation states */

    f->fft_pilot_cfg = kiss_fft_alloc (MPILOTFFT, 0, NULL, NULL);
    assert(f->fft_pilot_cfg != NULL);

    for(i=0; i<NPILOTBASEBAND; i++) {
	f->pilot_baseband1[i].real = f->pilot_baseband2[i].real = 0.0;
	f->pilot_baseband1[i].imag = f->pilot_baseband2[i].imag = 0.0;
    }
    f->pilot_lut_index = 0;
    f->prev_pilot_lut_index = 3*M;
    
    for(i=0; i<NRXDEC-1+M; i++) {
        f->rxdec_lpf_mem[i].real = 0.0;
        f->rxdec_lpf_mem[i].imag = 0.0;
    }

    for(i=0; i<NPILOTLPF; i++) {
	f->pilot_lpf1[i].real = f->pilot_lpf2[i].real = 0.0;
	f->pilot_lpf1[i].imag = f->pilot_lpf2[i].imag = 0.0;
    }

    f->foff = 0.0;
    f->foff_phase_rect.real = 1.0;
    f->foff_phase_rect.imag = 0.0;

    for(i=0; i<NFILTER+M; i++) {
        f->rx_fdm_mem[i].real = 0.0; 
        f->rx_fdm_mem[i].imag = 0.0; 
    }

    f->fest_state = 0;
    f->sync = 0;
    f->timer = 0;
    for(i=0; i<NSYNC_MEM; i++)
        f->sync_mem[i] = 0;

    for(c=0; c<Nc+1; c++) {
	f->sig_est[c] = 0.0;
	f->noise_est[c] = 0.0;
    }

    for(i=0; i<2*FDMDV_NSPEC; i++)
	f->fft_buf[i] = 0.0;
    f->fft_cfg = kiss_fft_alloc (2*FDMDV_NSPEC, 0, NULL, NULL);
    assert(f->fft_cfg != NULL);
    
    f->sig_pwr_av = 0.0;

    return f;
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: fdmdv_destroy	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 16/4/2012

  Destroy an instance of the modem.

\*---------------------------------------------------------------------------*/

void fdmdv_destroy(struct FDMDV *fdmdv)
{
    assert(fdmdv != NULL);
    KISS_FFT_FREE(fdmdv->fft_pilot_cfg);
    KISS_FFT_FREE(fdmdv->fft_cfg);
    free(fdmdv->rx_test_bits_mem);
    free(fdmdv);
}


void fdmdv_use_old_qpsk_mapping(struct FDMDV *fdmdv) {
    fdmdv->old_qpsk_mapping = 1;  
}


int fdmdv_bits_per_frame(struct FDMDV *fdmdv)
{
    return (fdmdv->Nc * NB);
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: fdmdv_get_test_bits()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 16/4/2012

  Generate a frame of bits from a repeating sequence of random data.  OK so
  it's not very random if it repeats but it makes syncing at the demod easier
  for test purposes.

\*---------------------------------------------------------------------------*/

void fdmdv_get_test_bits(struct FDMDV *f, int tx_bits[])
{
    int i;
    int bits_per_frame = fdmdv_bits_per_frame(f);

    for(i=0; i<bits_per_frame; i++) {
	tx_bits[i] = test_bits[f->current_test_bit];
	f->current_test_bit++;
	if (f->current_test_bit > (f->ntest_bits-1))
	    f->current_test_bit = 0;
    }
 }

float fdmdv_get_fsep(struct FDMDV *f)
{
    return f->fsep;
}

void fdmdv_set_fsep(struct FDMDV *f, float fsep) {
    int   c;
    float carrier_freq;

    f->fsep = fsep;

    /* Set up frequency of each carrier */

    for(c=0; c<f->Nc/2; c++) {
	carrier_freq = (-f->Nc/2 + c)*f->fsep;
	f->freq[c].real = cosf(2.0*PI*carrier_freq/FS);
 	f->freq[c].imag = sinf(2.0*PI*carrier_freq/FS);
 	f->freq_pol[c]  = 2.0*PI*carrier_freq/FS;
    }

    for(c=f->Nc/2; c<f->Nc; c++) {
	carrier_freq = (-f->Nc/2 + c + 1)*f->fsep;
	f->freq[c].real = cosf(2.0*PI*carrier_freq/FS);
 	f->freq[c].imag = sinf(2.0*PI*carrier_freq/FS);
 	f->freq_pol[c]  = 2.0*PI*carrier_freq/FS;
    }
}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: bits_to_dqpsk_symbols()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 16/4/2012

  Maps bits to parallel DQPSK symbols. Generate Nc+1 QPSK symbols from
  vector of (1,Nc*Nb) input tx_bits.  The Nc+1 symbol is the +1 -1 +1
  .... BPSK sync carrier.

\*---------------------------------------------------------------------------*/

void bits_to_dqpsk_symbols(COMP tx_symbols[], int Nc, COMP prev_tx_symbols[], int tx_bits[], int *pilot_bit, int old_qpsk_mapping)
{
    int c, msb, lsb;
    COMP j = {0.0,1.0};

    /* Map tx_bits to to Nc DQPSK symbols.  Note legacy support for
       old (suboptimal) V0.91 FreeDV mapping */

    for(c=0; c<Nc; c++) {
	msb = tx_bits[2*c]; 
	lsb = tx_bits[2*c+1];
	if ((msb == 0) && (lsb == 0))
	    tx_symbols[c] = prev_tx_symbols[c];
	if ((msb == 0) && (lsb == 1))
            tx_symbols[c] = cmult(j, prev_tx_symbols[c]);
	if ((msb == 1) && (lsb == 0)) {
	    if (old_qpsk_mapping)
                tx_symbols[c] = cneg(prev_tx_symbols[c]);           
            else
                tx_symbols[c] = cmult(cneg(j),prev_tx_symbols[c]);
        }
	if ((msb == 1) && (lsb == 1)) {
	    if (old_qpsk_mapping)
                tx_symbols[c] = cmult(cneg(j),prev_tx_symbols[c]);  
            else
                tx_symbols[c] = cneg(prev_tx_symbols[c]);
        }
    }

    /* +1 -1 +1 -1 BPSK sync carrier, once filtered becomes (roughly)
       two spectral lines at +/- Rs/2 */
 
    if (*pilot_bit)
	tx_symbols[Nc] = cneg(prev_tx_symbols[Nc]);
    else
	tx_symbols[Nc] = prev_tx_symbols[Nc];

    if (*pilot_bit) 
	*pilot_bit = 0;
    else
	*pilot_bit = 1;
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: tx_filter()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 17/4/2012

  Given Nc*NB bits construct M samples (1 symbol) of Nc+1 filtered
  symbols streams.

\*---------------------------------------------------------------------------*/

void tx_filter(COMP tx_baseband[NC+1][M], int Nc, COMP tx_symbols[], COMP tx_filter_memory[NC+1][NSYM])
{
    int     c;
    int     i,j,k;
    float   acc;
    COMP    gain;

    gain.real = sqrtf(2.0)/2.0;
    gain.imag = 0.0;
    
    for(c=0; c<Nc+1; c++)
	tx_filter_memory[c][NSYM-1] = cmult(tx_symbols[c], gain);
    
    /* 
       tx filter each symbol, generate M filtered output samples for each symbol.
       Efficient polyphase filter techniques used as tx_filter_memory is sparse
    */

    for(i=0; i<M; i++) {
	for(c=0; c<Nc+1; c++) {

	    /* filter real sample of symbol for carrier c */

	    acc = 0.0;
	    for(j=0,k=M-i-1; j<NSYM; j++,k+=M)
		acc += M * tx_filter_memory[c][j].real * gt_alpha5_root[k];
	    tx_baseband[c][i].real = acc;	

	    /* filter imag sample of symbol for carrier c */

	    acc = 0.0;
	    for(j=0,k=M-i-1; j<NSYM; j++,k+=M)
		acc += M * tx_filter_memory[c][j].imag * gt_alpha5_root[k];
	    tx_baseband[c][i].imag = acc;

	}
    }

    /* shift memory, inserting zeros at end */

    for(i=0; i<NSYM-1; i++)
	for(c=0; c<Nc+1; c++)
	    tx_filter_memory[c][i] = tx_filter_memory[c][i+1];

    for(c=0; c<Nc+1; c++) {
	tx_filter_memory[c][NSYM-1].real = 0.0;
	tx_filter_memory[c][NSYM-1].imag = 0.0;
    }
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: fdm_upconvert()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 17/4/2012

  Construct FDM signal by frequency shifting each filtered symbol
  stream.  Returns complex signal so we can apply frequency offsets
  easily.

\*---------------------------------------------------------------------------*/

void fdm_upconvert(COMP tx_fdm[], int Nc, COMP tx_baseband[NC+1][M], COMP phase_tx[], COMP freq[],
                   COMP *fbb_phase, COMP fbb_rect)
{
    int   i,c;
    COMP  two = {2.0, 0.0};
    COMP  pilot;
    float mag;

    for(i=0; i<M; i++) {
	tx_fdm[i].real = 0.0;
	tx_fdm[i].imag = 0.0;
    }

    /* Nc/2 tones below centre freq */
  
    for (c=0; c<Nc/2; c++) 
	for (i=0; i<M; i++) {
	    phase_tx[c] = cmult(phase_tx[c], freq[c]);
	    tx_fdm[i] = cadd(tx_fdm[i], cmult(tx_baseband[c][i], phase_tx[c]));
	}

    /* Nc/2 tones above centre freq */

    for (c=Nc/2; c<Nc; c++) 
	for (i=0; i<M; i++) {
	    phase_tx[c] = cmult(phase_tx[c], freq[c]);
	    tx_fdm[i] = cadd(tx_fdm[i], cmult(tx_baseband[c][i], phase_tx[c]));
	}

    /* add centre pilot tone  */

    c = Nc;
    for (i=0; i<M; i++) {
	phase_tx[c] = cmult(phase_tx[c],  freq[c]);
	pilot = cmult(cmult(two, tx_baseband[c][i]), phase_tx[c]);
	tx_fdm[i] = cadd(tx_fdm[i], pilot);
    }

    /* shift whole thing up to carrier freq */

    for (i=0; i<M; i++) {
	*fbb_phase = cmult(*fbb_phase, fbb_rect);
	tx_fdm[i] = cmult(tx_fdm[i], *fbb_phase);
    }

    /*
      Scale such that total Carrier power C of real(tx_fdm) = Nc.  This
      excludes the power of the pilot tone.
      We return the complex (single sided) signal to make frequency
      shifting for the purpose of testing easier
    */

    for (i=0; i<M; i++) 
	tx_fdm[i] = cmult(two, tx_fdm[i]);

    /* normalise digital oscilators as the magnitude can drfift over time */

    for (c=0; c<Nc+1; c++) {
        mag = cabsolute(phase_tx[c]);
	phase_tx[c].real /= mag;	
	phase_tx[c].imag /= mag;	
    }

    mag = cabsolute(*fbb_phase);
    fbb_phase->real /= mag;	
    fbb_phase->imag /= mag;	
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: fdmdv_mod()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 26/4/2012

  FDMDV modulator, take a frame of FDMDV_BITS_PER_FRAME bits and
  generates a frame of FDMDV_SAMPLES_PER_FRAME modulated symbols.
  Sync bit is returned to aid alignment of your next frame.  

  The sync_bit value returned will be used for the _next_ frame.

  The output signal is complex to support single sided frequency
  shifting, for example when testing frequency offsets in channel
  simulation.

\*---------------------------------------------------------------------------*/

void fdmdv_mod(struct FDMDV *fdmdv, COMP tx_fdm[], int tx_bits[], int *sync_bit)
{
    COMP          tx_symbols[NC+1];
    COMP          tx_baseband[NC+1][M];
    TIMER_VAR(mod_start, tx_filter_start, fdm_upconvert_start);

    TIMER_SAMPLE(mod_start);
    bits_to_dqpsk_symbols(tx_symbols, fdmdv->Nc, fdmdv->prev_tx_symbols, tx_bits, &fdmdv->tx_pilot_bit, fdmdv->old_qpsk_mapping);
    memcpy(fdmdv->prev_tx_symbols, tx_symbols, sizeof(COMP)*(fdmdv->Nc+1));
    TIMER_SAMPLE_AND_LOG(tx_filter_start, mod_start, "    bits_to_dqpsk_symbols"); 
    tx_filter(tx_baseband, fdmdv->Nc, tx_symbols, fdmdv->tx_filter_memory);
    TIMER_SAMPLE_AND_LOG(fdm_upconvert_start, tx_filter_start, "    tx_filter"); 
    fdm_upconvert(tx_fdm, fdmdv->Nc, tx_baseband, fdmdv->phase_tx, fdmdv->freq, &fdmdv->fbb_phase_tx, fdmdv->fbb_rect);
    TIMER_SAMPLE_AND_LOG2(fdm_upconvert_start, "    fdm_upconvert"); 

    *sync_bit = fdmdv->tx_pilot_bit;
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: generate_pilot_fdm()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 19/4/2012

  Generate M samples of DBPSK pilot signal for Freq offset estimation.

\*---------------------------------------------------------------------------*/

void generate_pilot_fdm(COMP *pilot_fdm, int *bit, float *symbol, 
			float *filter_mem, COMP *phase, COMP *freq)
{
    int   i,j,k;
    float tx_baseband[M];

    /* +1 -1 +1 -1 DBPSK sync carrier, once filtered becomes (roughly)
       two spectral lines at +/- RS/2 */
 
    if (*bit)
	*symbol = -*symbol;
    else
	*symbol = *symbol;
    if (*bit) 
	*bit = 0;
    else
	*bit = 1;

    /* filter DPSK symbol to create M baseband samples */

    filter_mem[NFILTER-1] = (sqrtf(2)/2) * *symbol;
    for(i=0; i<M; i++) {
	tx_baseband[i] = 0.0; 
	for(j=M-1,k=M-i-1; j<NFILTER; j+=M,k+=M)
	    tx_baseband[i] += M * filter_mem[j] * gt_alpha5_root[k];
    }

    /* shift memory, inserting zeros at end */

    for(i=0; i<NFILTER-M; i++)
	filter_mem[i] = filter_mem[i+M];

    for(i=NFILTER-M; i<NFILTER; i++)
	filter_mem[i] = 0.0;

    /* upconvert */

    for(i=0; i<M; i++) {
	*phase = cmult(*phase, *freq);
	pilot_fdm[i].real = sqrtf(2)*2*tx_baseband[i] * phase->real;
	pilot_fdm[i].imag = sqrtf(2)*2*tx_baseband[i] * phase->imag;
    }
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: generate_pilot_lut()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 19/4/2012

  Generate a 4M sample vector of DBPSK pilot signal.  As the pilot signal
  is periodic in 4M samples we can then use this vector as a look up table
  for pilot signal generation in the demod.

\*---------------------------------------------------------------------------*/

void generate_pilot_lut(COMP pilot_lut[], COMP *pilot_freq)
{
    int   pilot_rx_bit = 0;
    float pilot_symbol = sqrtf(2.0);
    COMP  pilot_phase  = {1.0, 0.0};
    float pilot_filter_mem[NFILTER];
    COMP  pilot[M];
    int   i,f;

    for(i=0; i<NFILTER; i++)
	pilot_filter_mem[i] = 0.0;

    /* discard first 4 symbols as filter memory is filling, just keep
       last four symbols */

    for(f=0; f<8; f++) {
	generate_pilot_fdm(pilot, &pilot_rx_bit, &pilot_symbol, pilot_filter_mem, &pilot_phase, pilot_freq);
	if (f >= 4)
	    memcpy(&pilot_lut[M*(f-4)], pilot, M*sizeof(COMP));
    }

}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: lpf_peak_pick()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 20/4/2012

  LPF and peak pick part of freq est, put in a function as we call it twice.

\*---------------------------------------------------------------------------*/

void lpf_peak_pick(float *foff, float *max, COMP pilot_baseband[], 
		   COMP pilot_lpf[], kiss_fft_cfg fft_pilot_cfg, COMP S[], int nin)
{
    int   i,j,k;
    int   mpilot;
    COMP  s[MPILOTFFT];
    float mag, imax;
    int   ix;
    float r;

    /* LPF cutoff 200Hz, so we can handle max +/- 200 Hz freq offset */

    for(i=0; i<NPILOTLPF-nin; i++)
	pilot_lpf[i] = pilot_lpf[nin+i];
    for(i=NPILOTLPF-nin, j=NPILOTBASEBAND-nin; i<NPILOTLPF; i++,j++) {
	pilot_lpf[i].real = 0.0; pilot_lpf[i].imag = 0.0;
	for(k=0; k<NPILOTCOEFF; k++)
	    pilot_lpf[i] = cadd(pilot_lpf[i], fcmult(pilot_coeff[k], pilot_baseband[j-NPILOTCOEFF+1+k]));
    }

    /* decimate to improve DFT resolution, window and DFT */

    mpilot = FS/(2*200);  /* calc decimation rate given new sample rate is twice LPF freq */
    for(i=0; i<MPILOTFFT; i++) {
	s[i].real = 0.0; s[i].imag = 0.0;
    }
    for(i=0,j=0; i<NPILOTLPF; i+=mpilot,j++) {
	s[j] = fcmult(hanning[i], pilot_lpf[i]); 
    }

    kiss_fft(fft_pilot_cfg, (kiss_fft_cpx *)s, (kiss_fft_cpx *)S);

    /* peak pick and convert to Hz */

    imax = 0.0;
    ix = 0;
    for(i=0; i<MPILOTFFT; i++) {
	mag = S[i].real*S[i].real + S[i].imag*S[i].imag;
	if (mag > imax) {
	    imax = mag;
	    ix = i;
	}
    }
    r = 2.0*200.0/MPILOTFFT;     /* maps FFT bin to frequency in Hz */
  
    if (ix >= MPILOTFFT/2)
	*foff = (ix - MPILOTFFT)*r;
    else
	*foff = (ix)*r;
    *max = imax;

}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: rx_est_freq_offset()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 19/4/2012

  Estimate frequency offset of FDM signal using BPSK pilot.  Note that
  this algorithm is quite sensitive to pilot tone level wrt other
  carriers, so test variations to the pilot amplitude carefully.

\*---------------------------------------------------------------------------*/

float rx_est_freq_offset(struct FDMDV *f, COMP rx_fdm[], int nin)
{
    int  i,j;
    COMP pilot[M+M/P];
    COMP prev_pilot[M+M/P];
    float foff, foff1, foff2;
    float   max1, max2;

    assert(nin <= M+M/P);

    /* get pilot samples used for correlation/down conversion of rx signal */

    for (i=0; i<nin; i++) {
	pilot[i] = f->pilot_lut[f->pilot_lut_index];
	f->pilot_lut_index++;
	if (f->pilot_lut_index >= 4*M)
	    f->pilot_lut_index = 0;
	
	prev_pilot[i] = f->pilot_lut[f->prev_pilot_lut_index];
	f->prev_pilot_lut_index++;
	if (f->prev_pilot_lut_index >= 4*M)
	    f->prev_pilot_lut_index = 0;
    }

    /*
      Down convert latest M samples of pilot by multiplying by ideal
      BPSK pilot signal we have generated locally.  The peak of the
      resulting signal is sensitive to the time shift between the
      received and local version of the pilot, so we do it twice at
      different time shifts and choose the maximum.
    */

    for(i=0; i<NPILOTBASEBAND-nin; i++) {
	f->pilot_baseband1[i] = f->pilot_baseband1[i+nin];
	f->pilot_baseband2[i] = f->pilot_baseband2[i+nin];
    }

    for(i=0,j=NPILOTBASEBAND-nin; i<nin; i++,j++) {
       	f->pilot_baseband1[j] = cmult(rx_fdm[i], cconj(pilot[i]));
	f->pilot_baseband2[j] = cmult(rx_fdm[i], cconj(prev_pilot[i]));
    }

    lpf_peak_pick(&foff1, &max1, f->pilot_baseband1, f->pilot_lpf1, f->fft_pilot_cfg, f->S1, nin);
    lpf_peak_pick(&foff2, &max2, f->pilot_baseband2, f->pilot_lpf2, f->fft_pilot_cfg, f->S2, nin);

    if (max1 > max2)
	foff = foff1;
    else
	foff = foff2;
	
    return foff;
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: fdmdv_freq_shift()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 26/4/2012

  Frequency shift modem signal.  The use of complex input and output allows
  single sided frequency shifting (no images).

\*---------------------------------------------------------------------------*/

void fdmdv_freq_shift(COMP rx_fdm_fcorr[], COMP rx_fdm[], float foff, 
                      COMP *foff_phase_rect, int nin)
{
    COMP  foff_rect;
    float mag;
    int   i;

    foff_rect.real = cosf(2.0*PI*foff/FS);
    foff_rect.imag = sinf(2.0*PI*foff/FS);
    for(i=0; i<nin; i++) {
	*foff_phase_rect = cmult(*foff_phase_rect, foff_rect);
	rx_fdm_fcorr[i] = cmult(rx_fdm[i], *foff_phase_rect);
    }

    /* normalise digital oscilator as the magnitude can drfift over time */

    mag = cabsolute(*foff_phase_rect);
    foff_phase_rect->real /= mag;	 
    foff_phase_rect->imag /= mag;	 
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: fdm_downconvert     
  AUTHOR......: David Rowe			      
  DATE CREATED: 22/4/2012

  Frequency shift each modem carrier down to Nc+1 baseband signals.

\*---------------------------------------------------------------------------*/

void fdm_downconvert(COMP rx_baseband[NC+1][M+M/P], int Nc, COMP rx_fdm[], COMP phase_rx[], COMP freq[], int nin)
{
    int   i,c;
    float mag;

    /* maximum number of input samples to demod */

    assert(nin <= (M+M/P));

    /* Nc/2 tones below centre freq */
  
    for (c=0; c<Nc/2; c++) 
	for (i=0; i<nin; i++) {
	    phase_rx[c] = cmult(phase_rx[c], freq[c]);
	    rx_baseband[c][i] = cmult(rx_fdm[i], cconj(phase_rx[c]));
	}

    /* Nc/2 tones above centre freq */

    for (c=Nc/2; c<Nc; c++) 
	for (i=0; i<nin; i++) {
	    phase_rx[c] = cmult(phase_rx[c], freq[c]);
	    rx_baseband[c][i] = cmult(rx_fdm[i], cconj(phase_rx[c]));
	}

    /* centre pilot tone  */

    c = Nc;
    for (i=0; i<nin; i++) {
	phase_rx[c] = cmult(phase_rx[c],  freq[c]);
	rx_baseband[c][i] = cmult(rx_fdm[i], cconj(phase_rx[c]));
    }

    /* normalise digital oscilators as the magnitude can drift over time */

    for (c=0; c<Nc+1; c++) {
        mag = cabsolute(phase_rx[c]);
	phase_rx[c].real /= mag;	  
	phase_rx[c].imag /= mag;	  
    }
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: rx_filter()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 22/4/2012

  Receive filter each baseband signal at oversample rate P.  Filtering at
  rate P lowers CPU compared to rate M.

  Depending on the number of input samples to the demod nin, we
  produce P-1, P (usually), or P+1 filtered samples at rate P.  nin is
  occasionally adjusted to compensate for timing slips due to
  different tx and rx sample clocks.

\*---------------------------------------------------------------------------*/

void rx_filter(COMP rx_filt[NC+1][P+1], int Nc, COMP rx_baseband[NC+1][M+M/P], COMP rx_filter_memory[NC+1][NFILTER], int nin)
{
    int c, i,j,k,l;
    int n=M/P;

    /* rx filter each symbol, generate P filtered output samples for
       each symbol.  Note we keep filter memory at rate M, it's just
       the filter output at rate P */

    for(i=0, j=0; i<nin; i+=n,j++) {

	/* latest input sample */
	
	for(c=0; c<Nc+1; c++)
	    for(k=NFILTER-n,l=i; k<NFILTER; k++,l++)	
		rx_filter_memory[c][k] = rx_baseband[c][l];
	
	/* convolution (filtering) */

	for(c=0; c<Nc+1; c++) {
	    rx_filt[c][j].real = 0.0; rx_filt[c][j].imag = 0.0;
	    for(k=0; k<NFILTER; k++) 
		rx_filt[c][j] = cadd(rx_filt[c][j], fcmult(gt_alpha5_root[k], rx_filter_memory[c][k]));
	}

	/* make room for next input sample */
	
	for(c=0; c<Nc+1; c++)
	    for(k=0,l=n; k<NFILTER-n; k++,l++)	
		rx_filter_memory[c][k] = rx_filter_memory[c][l];
    }

    assert(j <= (P+1)); /* check for any over runs */
}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: rxdec_filter()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 31 July 2014

  +/- 1000Hz low pass filter, allows us to filter at rate Q to save CPU load.

\*---------------------------------------------------------------------------*/

void rxdec_filter(COMP rx_fdm_filter[], COMP rx_fdm[], COMP rxdec_lpf_mem[], int nin) {
    int i,j,k; 

    for(i=0; i<NRXDEC-1+M-nin; i++)
        rxdec_lpf_mem[i] = rxdec_lpf_mem[i+nin];
    for(i=0, j=NRXDEC-1+M-nin; i<nin; i++,j++)
        rxdec_lpf_mem[j] = rx_fdm[i];

    for(i=0; i<nin; i++) {
        rx_fdm_filter[i].real = 0.0;
        for(k=0; k<NRXDEC; k++)
            rx_fdm_filter[i].real += rxdec_lpf_mem[i+k].real * rxdec_coeff[k];
        rx_fdm_filter[i].imag = 0.0;
        for(k=0; k<NRXDEC; k++)
            rx_fdm_filter[i].imag += rxdec_lpf_mem[i+k].imag * rxdec_coeff[k];
    }
}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: fir_filter()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: July 2014

  Helper fir filter function.

\*---------------------------------------------------------------------------*/

static float fir_filter(float mem[], float coeff[], int dec_rate) {
    float acc = 0.0;
    int   m;

    for(m=0; m<NFILTER; m+=dec_rate) { 
        acc += coeff[m] * mem[2*m];
    }

    return dec_rate*acc;
}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: down_convert_and_rx_filter()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 30/6/2014

  Combined down convert and rx filter, more memory efficient but less
  intuitive design.  

  Depending on the number of input samples to the demod nin, we
  produce P-1, P (usually), or P+1 filtered samples at rate P.  nin is
  occasionally adjusted to compensate for timing slips due to
  different tx and rx sample clocks.

\*---------------------------------------------------------------------------*/

/*
   TODO: [ ] windback phase calculated once at init time
*/

void down_convert_and_rx_filter(COMP rx_filt[NC+1][P+1], int Nc, COMP rx_fdm[], 
                                COMP rx_fdm_mem[], COMP phase_rx[], COMP freq[], 
                                float freq_pol[], int nin, int dec_rate)
{
    int i,k,c,st,N;
    float windback_phase, mag;
    COMP  windback_phase_rect;
    COMP  rx_baseband[NFILTER+M];
    COMP  f_rect;

    //TIMER_VAR(windback_start,  downconvert_start, filter_start);

    /* update memory of rx_fdm */

    for(i=0; i<NFILTER+M-nin; i++)
        rx_fdm_mem[i] = rx_fdm_mem[i+nin];
    for(i=NFILTER+M-nin, k=0; i<NFILTER+M; i++, k++)
        rx_fdm_mem[i] = rx_fdm[k];

    for(c=0; c<Nc+1; c++) {

        /*
          now downconvert using current freq offset to get Nfilter+nin
          baseband samples.
     
                     Nfilter             nin
          |--------------------------|---------|
                                      | 
                                  phase_rx(c)
     
          This means winding phase(c) back from this point to ensure
          phase continuity.
        */

        //TIMER_SAMPLE(windback_start);
        windback_phase           = -freq_pol[c]*NFILTER;
        windback_phase_rect.real = cosf(windback_phase);
        windback_phase_rect.imag = sinf(windback_phase);
        phase_rx[c]              = cmult(phase_rx[c],windback_phase_rect);
        //TIMER_SAMPLE_AND_LOG(downconvert_start, windback_start, "        windback"); 
    
        /* down convert all samples in buffer */

        st  = NFILTER+M-1;    /* end of buffer                  */
        st -= nin-1;          /* first new sample               */
        st -= NFILTER;        /* first sample used in filtering */
        
        /* freq shift per dec_rate step is dec_rate times original shift */

        f_rect = freq[c];
        for(i=0; i<dec_rate-1; i++)
            f_rect = cmult(f_rect,freq[c]);

        for(i=st; i<NFILTER+M; i+=dec_rate) {
            phase_rx[c]    = cmult(phase_rx[c], f_rect);
            rx_baseband[i] = cmult(rx_fdm_mem[i],cconj(phase_rx[c]));
        }
        //TIMER_SAMPLE_AND_LOG(filter_start, downconvert_start, "        downconvert"); 
 
        /* now we can filter this carrier's P symbols */

        N=M/P;
        for(i=0, k=0; i<nin; i+=N, k++) {
           #ifdef ORIG 
           rx_filt[c][k].real = 0.0; rx_filt[c][k].imag = 0.0;
            
           for(m=0; m<NFILTER; m++) 
               rx_filt[c][k] = cadd(rx_filt[c][k], fcmult(gt_alpha5_root[m], rx_baseband[st+i+m]));
           #else
           rx_filt[c][k].real = fir_filter(&rx_baseband[st+i].real, (float*)gt_alpha5_root, dec_rate);
           rx_filt[c][k].imag = fir_filter(&rx_baseband[st+i].imag, (float*)gt_alpha5_root, dec_rate);
           #endif
        }
        //TIMER_SAMPLE_AND_LOG2(filter_start, "        filter"); 

        /* normalise digital oscilators as the magnitude can drift over time */

        mag = cabsolute(phase_rx[c]);
	phase_rx[c].real /= mag;	  
	phase_rx[c].imag /= mag;
	  
       //printf("phase_rx[%d] = %f %f\n", c, phase_rx[c].real, phase_rx[c].imag);
    }
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: rx_est_timing()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 23/4/2012

  Estimate optimum timing offset, re-filter receive symbols at optimum
  timing estimate.

\*---------------------------------------------------------------------------*/

float rx_est_timing(COMP rx_symbols[], 
                    int  Nc,
		    COMP rx_filt[NC+1][P+1], 
		    COMP rx_filter_mem_timing[NC+1][NT*P], 
		    float env[],
		    int nin)	 
{
    int   c,i,j;
    int   adjust;
    COMP  x, phase, freq;
    float rx_timing, fract, norm_rx_timing;
    int   low_sample, high_sample;

    /*
      nin  adjust 
      --------------------------------
      120  -1 (one less rate P sample)
      160   0 (nominal)
      200   1 (one more rate P sample)
    */

    adjust = P - nin*P/M;
    
    /* update buffer of NT rate P filtered symbols */
    
    for(c=0; c<Nc+1; c++) 
	for(i=0,j=P-adjust; i<(NT-1)*P+adjust; i++,j++)
	    rx_filter_mem_timing[c][i] = rx_filter_mem_timing[c][j];
    for(c=0; c<Nc+1; c++) 
	for(i=(NT-1)*P+adjust,j=0; i<NT*P; i++,j++)
	    rx_filter_mem_timing[c][i] = rx_filt[c][j];
	    
    /* sum envelopes of all carriers */

    for(i=0; i<NT*P; i++) {
	env[i] = 0.0;
	for(c=0; c<Nc+1; c++)
	    env[i] += cabsolute(rx_filter_mem_timing[c][i]);
    }

    /* The envelope has a frequency component at the symbol rate.  The
       phase of this frequency component indicates the timing.  So work
       out single DFT at frequency 2*pi/P */

    x.real = 0.0; x.imag = 0.0;
    freq.real = cosf(2*PI/P);
    freq.imag = sinf(2*PI/P);
    phase.real = 1.0;
    phase.imag = 0.0;

    for(i=0; i<NT*P; i++) {
	x = cadd(x, fcmult(env[i], phase));
	phase = cmult(phase, freq);
    }

    /* Map phase to estimated optimum timing instant at rate P.  The
       P/4 part was adjusted by experiment, I know not why.... */
    
    norm_rx_timing = atan2f(x.imag, x.real)/(2*PI);
    rx_timing      = norm_rx_timing*P + P/4;
    
    if (rx_timing > P)
	rx_timing -= P;
    if (rx_timing < -P)
	rx_timing += P;

    /* rx_filter_mem_timing contains Nt*P samples (Nt symbols at rate
       P), where Nt is odd.  Lets use linear interpolation to resample
       in the centre of the timing estimation window .*/

    rx_timing  += floorf(NT/2.0)*P;
    low_sample = floorf(rx_timing);
    fract = rx_timing - low_sample;
    high_sample = ceilf(rx_timing);

    //printf("rx_timing: %f low_sample: %d high_sample: %d fract: %f\n", rx_timing, low_sample, high_sample, fract);
    
    for(c=0; c<Nc+1; c++) {
        rx_symbols[c] = cadd(fcmult(1.0-fract, rx_filter_mem_timing[c][low_sample-1]), fcmult(fract, rx_filter_mem_timing[c][high_sample-1]));
        //rx_symbols[c] = rx_filter_mem_timing[c][high_sample];
    }
  	
    return norm_rx_timing*M;
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: qpsk_to_bits()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 24/4/2012

  Convert DQPSK symbols back to an array of bits, extracts sync bit
  from DBPSK pilot, and also uses pilot to estimate fine frequency
  error.

\*---------------------------------------------------------------------------*/

float qpsk_to_bits(int rx_bits[], int *sync_bit, int Nc, COMP phase_difference[], COMP prev_rx_symbols[], 
                   COMP rx_symbols[], int old_qpsk_mapping)
{
    int   c;
    COMP  pi_on_4;
    COMP  d;
    int   msb=0, lsb=0;
    float ferr, norm;

    pi_on_4.real = cosf(PI/4.0);
    pi_on_4.imag = sinf(PI/4.0);

    /* Extra 45 degree clockwise lets us use real and imag axis as
       decision boundaries. "norm" makes sure the phase subtraction
       from the previous symbol doesn't affect the amplitude, which
       leads to sensible scatter plots */

    for(c=0; c<Nc; c++) {
        norm = 1.0/(cabsolute(prev_rx_symbols[c])+1E-6);
	phase_difference[c] = cmult(cmult(rx_symbols[c], fcmult(norm,cconj(prev_rx_symbols[c]))), pi_on_4);
    }
				    
    /* map (Nc,1) DQPSK symbols back into an (1,Nc*Nb) array of bits */

    for (c=0; c<Nc; c++) {
      d = phase_difference[c];
      if ((d.real >= 0) && (d.imag >= 0)) {
          msb = 0; lsb = 0;
      }
      if ((d.real < 0) && (d.imag >= 0)) {
          msb = 0; lsb = 1;
      }
      if ((d.real < 0) && (d.imag < 0)) {
          if (old_qpsk_mapping) {
              msb = 1; lsb = 0;
          } else {
              msb = 1; lsb = 1;
          }
      }
      if ((d.real >= 0) && (d.imag < 0)) {
          if (old_qpsk_mapping) {
              msb = 1; lsb = 1;
          } else {
              msb = 1; lsb = 0;
          }
      }
      rx_bits[2*c] = msb;
      rx_bits[2*c+1] = lsb;
    }
 
    /* Extract DBPSK encoded Sync bit and fine freq offset estimate */

    norm = 1.0/(cabsolute(prev_rx_symbols[Nc])+1E-6);
    phase_difference[Nc] = cmult(rx_symbols[Nc], fcmult(norm, cconj(prev_rx_symbols[Nc])));
    if (phase_difference[Nc].real < 0) {
      *sync_bit = 1;
      ferr = phase_difference[Nc].imag*norm;    /* make f_err magnitude insensitive */
    }
    else {
      *sync_bit = 0;
      ferr = -phase_difference[Nc].imag*norm;
    }
    
    /* pilot carrier gets an extra pi/4 rotation to make it consistent
       with other carriers, as we need it for snr_update and scatter
       diagram */

    phase_difference[Nc] = cmult(phase_difference[Nc], pi_on_4);

    return ferr;
}
 
/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: snr_update()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 17 May 2012

  Given phase differences update estimates of signal and noise levels.

\*---------------------------------------------------------------------------*/

void snr_update(float sig_est[], float noise_est[], int Nc, COMP phase_difference[])
{
    float s[NC+1];
    COMP  refl_symbols[NC+1];
    float n[NC+1];
    COMP  pi_on_4;
    int   c;

    pi_on_4.real = cosf(PI/4.0);
    pi_on_4.imag = sinf(PI/4.0);

    /* mag of each symbol is distance from origin, this gives us a
       vector of mags, one for each carrier. */

    for(c=0; c<Nc+1; c++)
	s[c] = cabsolute(phase_difference[c]);

    /* signal mag estimate for each carrier is a smoothed version of
       instantaneous magntitude, this gives us a vector of smoothed
       mag estimates, one for each carrier. */

    for(c=0; c<Nc+1; c++)
	sig_est[c] = SNR_COEFF*sig_est[c] + (1.0 - SNR_COEFF)*s[c];

    /* noise mag estimate is distance of current symbol from average
       location of that symbol.  We reflect all symbols into the first
       quadrant for convenience. */
    
    for(c=0; c<Nc+1; c++) {
	refl_symbols[c].real = fabsf(phase_difference[c].real);
	refl_symbols[c].imag = fabsf(phase_difference[c].imag);    
	n[c] = cabsolute(cadd(fcmult(sig_est[c], pi_on_4), cneg(refl_symbols[c])));
    }
     
    /* noise mag estimate for each carrier is a smoothed version of
       instantaneous noise mag, this gives us a vector of smoothed
       noise power estimates, one for each carrier. */

    for(c=0; c<Nc+1; c++)
	noise_est[c] = SNR_COEFF*noise_est[c] + (1 - SNR_COEFF)*n[c];
}

// returns number of shorts in error_pattern[], one short per error

int fdmdv_error_pattern_size(struct FDMDV *f) {
    return f->ntest_bits;
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: fdmdv_put_test_bits()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 24/4/2012

  Accepts nbits from rx and attempts to sync with test_bits sequence.
  If sync OK measures bit errors.

\*---------------------------------------------------------------------------*/

void fdmdv_put_test_bits(struct FDMDV *f, int *sync, short error_pattern[],
			 int *bit_errors, int *ntest_bits, int rx_bits[])
{
    int   i,j;
    float ber;
    int   bits_per_frame = fdmdv_bits_per_frame(f);

    /* Append to our memory */

    for(i=0,j=bits_per_frame; i<f->ntest_bits-bits_per_frame; i++,j++)
	f->rx_test_bits_mem[i] = f->rx_test_bits_mem[j];
    for(i=f->ntest_bits-bits_per_frame,j=0; i<f->ntest_bits; i++,j++)
	f->rx_test_bits_mem[i] = rx_bits[j];
    
    /* see how many bit errors we get when checked against test sequence */
       
    *bit_errors = 0;
    for(i=0; i<f->ntest_bits; i++) {
        error_pattern[i] = test_bits[i] ^ f->rx_test_bits_mem[i];
	*bit_errors += error_pattern[i];
	//printf("%d %d %d %d\n", i, test_bits[i], f->rx_test_bits_mem[i], test_bits[i] ^ f->rx_test_bits_mem[i]);
    }

    /* if less than a thresh we are aligned and in sync with test sequence */

    ber = (float)*bit_errors/f->ntest_bits;
  
    *sync = 0;
    if (ber < 0.2)
	*sync = 1;
   
    *ntest_bits = f->ntest_bits;
    
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: freq_state(()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 24/4/2012

  Freq offset state machine.  Moves between coarse and fine states
  based on BPSK pilot sequence.  Freq offset estimator occasionally
  makes mistakes when used continuously.  So we use it until we have
  acquired the BPSK pilot, then switch to a more robust "fine"
  tracking algorithm.  If we lose sync we switch back to coarse mode
  for fast re-acquisition of large frequency offsets.

  The sync state is also useful for higher layers to determine when
  there is valid FDMDV data for decoding.  We want to reliably and
  quickly get into sync, stay in sync even on fading channels, and
  fall out of sync quickly if tx stops or it's a false sync.

  In multipath fading channels the BPSK sync carrier may be pushed
  down in the noise, despite other carriers being at full strength.
  We want to avoid loss of sync in these cases.

\*---------------------------------------------------------------------------*/

int freq_state(int *reliable_sync_bit, int sync_bit, int *state, int *timer, int *sync_mem)
{
    int next_state, sync, unique_word, i, corr;

    /* look for 6 symbols (120ms) 101010 of sync sequence */

    unique_word = 0;
    for(i=0; i<NSYNC_MEM-1; i++)
        sync_mem[i] = sync_mem[i+1];
    sync_mem[i] = 1 - 2*sync_bit;
    corr = 0;
    for(i=0; i<NSYNC_MEM; i++)
        corr += sync_mem[i]*sync_uw[i];
    if (abs(corr) == NSYNC_MEM)
        unique_word = 1;
    *reliable_sync_bit = (corr == NSYNC_MEM);

    /* iterate state machine */

    next_state = *state;
    switch(*state) {
    case 0:
	if (unique_word) {
	    next_state = 1;
            *timer = 0;
        }
	break;
    case 1:                   /* tentative sync state         */
	if (unique_word) {
            (*timer)++;
            if (*timer == 25) /* sync has been good for 500ms */
                next_state = 2;
        }
	else 
	    next_state = 0;  /* quickly fall out of sync     */
	break;
    case 2:                  /* good sync state */
	if (unique_word == 0) {
            *timer = 0;
	    next_state = 3;
        }
	break;
    case 3:                  /* tentative bad state, but could be a fade */
	if (unique_word)
	    next_state = 2;
	else  {
            (*timer)++;
            if (*timer == 50) /* wait for 1000ms in case sync comes back  */
                next_state = 0;
        }
	break;
    }

    //printf("state: %d next_state: %d uw: %d timer: %d\n", *state, next_state, unique_word, *timer);
    *state = next_state;
    if (*state)
	sync = 1;
    else
	sync = 0;
 
    return sync;
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: fdmdv_demod()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 26/4/2012

  FDMDV demodulator, take an array of FDMDV_SAMPLES_PER_FRAME
  modulated samples, returns an array of FDMDV_BITS_PER_FRAME bits,
  plus the sync bit.  

  The input signal is complex to support single sided frequency shifting
  before the demod input (e.g. fdmdv2 click to tune feature).

  The number of input samples nin will normally be M ==
  FDMDV_SAMPLES_PER_FRAME.  However to adjust for differences in
  transmit and receive sample clocks nin will occasionally be M-M/P,
  or M+M/P.

\*---------------------------------------------------------------------------*/

void fdmdv_demod(struct FDMDV *fdmdv, int rx_bits[], 
		 int *reliable_sync_bit, COMP rx_fdm[], int *nin)
{
    float         foff_coarse, foff_fine;
    COMP          rx_fdm_fcorr[M+M/P];
    COMP          rx_fdm_filter[M+M/P];
    COMP          rx_fdm_bb[M+M/P];
    COMP          rx_filt[NC+1][P+1];
    COMP          rx_symbols[NC+1];
    float         env[NT*P];
    int           sync_bit;
    TIMER_VAR(demod_start, fdmdv_freq_shift_start, down_convert_and_rx_filter_start);
    TIMER_VAR(rx_est_timing_start, qpsk_to_bits_start, snr_update_start, freq_state_start);

    /* shift down to complex baseband */

    fdmdv_freq_shift(rx_fdm_bb, rx_fdm, -FDMDV_FCENTRE, &fdmdv->fbb_phase_rx, *nin);

    /* freq offset estimation and correction */
   
    TIMER_SAMPLE(demod_start);
    foff_coarse = rx_est_freq_offset(fdmdv, rx_fdm_bb, *nin);
    TIMER_SAMPLE_AND_LOG(fdmdv_freq_shift_start, demod_start, "    rx_est_freq_offset"); 
    
    if (fdmdv->sync == 0)
	fdmdv->foff = foff_coarse;
    fdmdv_freq_shift(rx_fdm_fcorr, rx_fdm_bb, -fdmdv->foff, &fdmdv->foff_phase_rect, *nin);
    TIMER_SAMPLE_AND_LOG(down_convert_and_rx_filter_start, fdmdv_freq_shift_start, "    fdmdv_freq_shift"); 
 	
    /* baseband processing */

    rxdec_filter(rx_fdm_filter, rx_fdm_fcorr, fdmdv->rxdec_lpf_mem, *nin);
    down_convert_and_rx_filter(rx_filt, fdmdv->Nc, rx_fdm_filter, fdmdv->rx_fdm_mem, fdmdv->phase_rx, fdmdv->freq, 
                               fdmdv->freq_pol, *nin, M/Q);
    TIMER_SAMPLE_AND_LOG(rx_est_timing_start, down_convert_and_rx_filter_start, "    down_convert_and_rx_filter"); 
    fdmdv->rx_timing = rx_est_timing(rx_symbols, fdmdv->Nc, rx_filt, fdmdv->rx_filter_mem_timing, env, *nin);	 
    TIMER_SAMPLE_AND_LOG(qpsk_to_bits_start, rx_est_timing_start, "    rx_est_timing"); 
    
    /* Adjust number of input samples to keep timing within bounds */

    *nin = M;

    if (fdmdv->rx_timing > 2*M/P)
	*nin += M/P;
    
    if (fdmdv->rx_timing < 0)
	*nin -= M/P;
    
    foff_fine = qpsk_to_bits(rx_bits, &sync_bit, fdmdv->Nc, fdmdv->phase_difference, fdmdv->prev_rx_symbols, rx_symbols, 
                             fdmdv->old_qpsk_mapping);
    memcpy(fdmdv->prev_rx_symbols, rx_symbols, sizeof(COMP)*(fdmdv->Nc+1));
    TIMER_SAMPLE_AND_LOG(snr_update_start, qpsk_to_bits_start, "    qpsk_to_bits"); 
    snr_update(fdmdv->sig_est, fdmdv->noise_est, fdmdv->Nc, fdmdv->phase_difference);
    TIMER_SAMPLE_AND_LOG(freq_state_start, snr_update_start, "    snr_update"); 

    /* freq offset estimation state machine */

    fdmdv->sync = freq_state(reliable_sync_bit, sync_bit, &fdmdv->fest_state, &fdmdv->timer, fdmdv->sync_mem);
    TIMER_SAMPLE_AND_LOG2(freq_state_start, "    freq_state"); 
    fdmdv->foff  -= TRACK_COEFF*foff_fine;
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: calc_snr()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 17 May 2012

  Calculate current SNR estimate (3000Hz noise BW)

\*---------------------------------------------------------------------------*/

float calc_snr(int Nc, float sig_est[], float noise_est[])
{
    float S, SdB;
    float mean, N50, N50dB, N3000dB;
    float snr_dB;
    int   c;
   
    S = 0.0;
    for(c=0; c<Nc+1; c++)
	S += powf(sig_est[c], 2.0);
    SdB = 10.0*log10f(S+1E-12);
    
    /* Average noise mag across all carriers and square to get an
       average noise power.  This is an estimate of the noise power in
       Rs = 50Hz of BW (note for raised root cosine filters Rs is the
       noise BW of the filter) */

    mean = 0.0;
    for(c=0; c<Nc+1; c++)
	mean += noise_est[c];
    mean /= (Nc+1);
    N50 = powf(mean, 2.0);
    N50dB = 10.0*log10f(N50+1E-12);

    /* Now multiply by (3000 Hz)/(50 Hz) to find the total noise power
       in 3000 Hz */

    N3000dB = N50dB + 10.0*log10f(3000.0/RS);

    snr_dB = SdB - N3000dB;

    return snr_dB;
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: fdmdv_get_demod_stats()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 1 May 2012

  Fills stats structure with a bunch of demod information.

\*---------------------------------------------------------------------------*/

void fdmdv_get_demod_stats(struct FDMDV *fdmdv, struct FDMDV_STATS *fdmdv_stats)
{
    int   c;

    fdmdv_stats->Nc = fdmdv->Nc;
    fdmdv_stats->snr_est = calc_snr(fdmdv->Nc, fdmdv->sig_est, fdmdv->noise_est);
    fdmdv_stats->sync = fdmdv->sync;
    fdmdv_stats->foff = fdmdv->foff;
    fdmdv_stats->rx_timing = fdmdv->rx_timing;
    fdmdv_stats->clock_offset = 0.0; /* TODO - implement clock offset estimation */

    for(c=0; c<fdmdv->Nc+1; c++) {
	fdmdv_stats->rx_symbols[c] = fdmdv->phase_difference[c];
    }
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: fdmdv_8_to_16()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 9 May 2012

  Changes the sample rate of a signal from 8 to 16 kHz.  Support function for
  SM1000.

\*---------------------------------------------------------------------------*/

void fdmdv_8_to_16(float out16k[], float in8k[], int n)
{
    int i,j,k,l;

    /* make sure n is an integer multiple of the oversampling rate, ow
       this function breaks */

    assert((n % FDMDV_OS) == 0);

    for(i=0; i<n; i++) {
	for(j=0; j<FDMDV_OS; j++) {
	    out16k[i*FDMDV_OS+j] = 0.0;
	    for(k=0,l=0; k<FDMDV_OS_TAPS_16K; k+=FDMDV_OS,l++)
		out16k[i*FDMDV_OS+j] += fdmdv_os_filter[k+j]*in8k[i-l];
	    out16k[i*FDMDV_OS+j] *= FDMDV_OS;
	    
	}
    }	

    /* update filter memory */

    for(i=-(FDMDV_OS_TAPS_8K); i<0; i++)
	in8k[i] = in8k[i + n];

}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: fdmdv_16_to_8()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 9 May 2012

  Changes the sample rate of a signal from 16 to 8 kHz.
 
  n is the number of samples at the 8 kHz rate, there are FDMDV_OS*n
  samples at the 48 kHz rate.  As above however a memory of
  FDMDV_OS_TAPS samples is reqd for in16k[] (see t16_8.c unit test as example).

  Low pass filter the 16 kHz signal at 4 kHz using the same filter as
  the upsampler, then just output every FDMDV_OS-th filtered sample.

\*---------------------------------------------------------------------------*/

void fdmdv_16_to_8(float out8k[], float in16k[], int n)
{
    int i,j;

    for(i=0; i<n; i++) {
	out8k[i] = 0.0;
	for(j=0; j<FDMDV_OS_TAPS_16K; j++)
	    out8k[i] += fdmdv_os_filter[j]*in16k[i*FDMDV_OS-j];
    }

    /* update filter memory */

    for(i=-FDMDV_OS_TAPS_16K; i<0; i++)
	in16k[i] = in16k[i + n*FDMDV_OS];
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: fdmdv_get_rx_spectrum()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 9 June 2012

  Returns the FDMDV_NSPEC point magnitude spectrum of the rx signal in
  dB. The spectral samples are scaled so that 0dB is the peak, a good
  range for plotting is 0 to -40dB.

  Note only the real part of the complex input signal is used at
  present.  A complex variable is used for input for compatability
  with the other rx signal procesing.

  Successive calls can be used to build up a waterfall or spectrogram
  plot, by mapping the received levels to colours.

  The time-frequency resolution of the spectrum can be adjusted by varying
  FDMDV_NSPEC.  Note that a 2*FDMDV_NSPEC size FFT is reqd to get
  FDMDV_NSPEC output points. FDMDV_NSPEC must be a power of 2.

  See octave/tget_spec.m for a demo real time spectral display using
  Octave. This demo averages the output over time to get a smoother
  display:

     av = 0.9*av + 0.1*mag_dB

\*---------------------------------------------------------------------------*/

void fdmdv_get_rx_spectrum(struct FDMDV *f, float mag_spec_dB[], 
					       COMP rx_fdm[], int nin) 
{
    int   i,j;
    COMP  fft_in[2*FDMDV_NSPEC];
    COMP  fft_out[2*FDMDV_NSPEC];
    float full_scale_dB;

    /* update buffer of input samples */

    for(i=0; i<2*FDMDV_NSPEC-nin; i++)
	f->fft_buf[i] = f->fft_buf[i+nin];
    for(j=0; j<nin; j++,i++)
	f->fft_buf[i] = rx_fdm[j].real;
    assert(i == 2*FDMDV_NSPEC);

    /* window and FFT */

    for(i=0; i<2*FDMDV_NSPEC; i++) {
	fft_in[i].real = f->fft_buf[i] * (0.5 - 0.5*cosf((float)i*2.0*PI/(2*FDMDV_NSPEC)));
	fft_in[i].imag = 0.0;
    }

    kiss_fft(f->fft_cfg, (kiss_fft_cpx *)fft_in, (kiss_fft_cpx *)fft_out);

    /* FFT scales up a signal of level 1 FDMDV_NSPEC */

    full_scale_dB = 20*log10(FDMDV_NSPEC);

    /* scale and convert to dB */

    for(i=0; i<FDMDV_NSPEC; i++) {
	mag_spec_dB[i]  = 10.0*log10f(fft_out[i].real*fft_out[i].real + fft_out[i].imag*fft_out[i].imag + 1E-12);
	mag_spec_dB[i] -= full_scale_dB;
    }
}

/*---------------------------------------------------------------------------*\
                                                       
  Function used during development to test if magnitude of digital
  oscillators was drifting.  It was!

\*---------------------------------------------------------------------------*/

void fdmdv_dump_osc_mags(struct FDMDV *f) 
{
    int   i;

    fprintf(stderr, "phase_tx[]:\n");
    for(i=0; i<=f->Nc; i++)
	fprintf(stderr,"  %1.3f", (double)cabsolute(f->phase_tx[i]));
    fprintf(stderr,"\nfreq[]:\n");
    for(i=0; i<=f->Nc; i++)
	fprintf(stderr,"  %1.3f", (double)cabsolute(f->freq[i]));
    fprintf(stderr,"\nfoff_phase_rect: %1.3f", (double)cabsolute(f->foff_phase_rect));
    fprintf(stderr,"\nphase_rx[]:\n");
    for(i=0; i<=f->Nc; i++)
	fprintf(stderr,"  %1.3f", (double)cabsolute(f->phase_rx[i]));
    fprintf(stderr, "\n\n");
}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: randn()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 2 August 2014

  Simple approximation to normal (gaussian) random number generator
  with 0 mean and unit variance.

\*---------------------------------------------------------------------------*/

#define RANDN_IT 12  /* This magic number of iterations gives us a
                        unit variance.  I think beacuse var =
                        (b-a)^2/12 for one uniform random variable, so
                        for a sum of n random variables it's
                        n(b-a)^2/12, or for b=1, a = 0, n=12, we get
                        var = 12(1-0)^2/12 = 1 */ 

static float randn() {
    int   i;
    float rn = 0.0;

    for(i=0; i<RANDN_IT; i++)
        rn += (float)rand()/RAND_MAX;

    rn -= (float)RANDN_IT/2.0;
    return rn;
}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: fdmdv_simulate_channel()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 10 July 2014

  Simple channel simulation function to aid in testing.  Target SNR
  uses noise measured in a 3 kHz bandwidth.  

  TODO: Measured SNR is coming out a few dB higher than target_snr, this
  needs to be fixed.

\*---------------------------------------------------------------------------*/

void fdmdv_simulate_channel(struct FDMDV *f, COMP samples[], int nin, float target_snr) 
{
    float sig_pwr, target_snr_linear, noise_pwr, noise_pwr_1Hz, noise_pwr_4000Hz, noise_gain;
    int   i;

    /* estimate signal power */

    sig_pwr = 0.0;
    for(i=0; i<nin; i++)
        sig_pwr += samples[i].real*samples[i].real + samples[i].imag*samples[i].imag;
    
    sig_pwr /= nin;

    f->sig_pwr_av = 0.9*f->sig_pwr_av + 0.1*sig_pwr;

    /* det noise to meet target SNR */

    target_snr_linear = powf(10.0, target_snr/10.0);
    noise_pwr = f->sig_pwr_av/target_snr_linear;       /* noise pwr in a 3000 Hz BW     */
    noise_pwr_1Hz = noise_pwr/3000.0;                  /* noise pwr in a 1 Hz bandwidth */
    noise_pwr_4000Hz = noise_pwr_1Hz*4000.0;           /* noise pwr in a 4000 Hz BW, which 
                                                          due to fs=8000 Hz in our simulation noise BW */

    noise_gain = sqrtf(0.5*noise_pwr_4000Hz);          /* split noise pwr between real and imag sides  */
    
    for(i=0; i<nin; i++) {
        samples[i].real += noise_gain*randn();
        samples[i].imag += noise_gain*randn();
    }
    /*    
    fprintf(stderr, "sig_pwr: %f f->sig_pwr_av: %e target_snr_linear: %f noise_pwr_4000Hz: %e noise_gain: %e\n", 
            sig_pwr, f->sig_pwr_av, target_snr_linear, noise_pwr_4000Hz, noise_gain);
    */
}
