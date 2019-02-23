#ifndef __CODEC2_INTERNAL__
#define __CODEC2_INTERNAL__

#include "codec2_fft.h"
//#include "newamp1.h"
//#include "newamp2.h"

struct CODEC2 {
    int           mode;
    C2CONST       c2const;
    int           Fs;
    int           n_samp;
    int           m_pitch;
    codec2_fft_cfg  fft_fwd_cfg;           /* forward FFT config                        */
//    codec2_fftr_cfg fftr_fwd_cfg;          /* forward real FFT config                   */
    float        *w;	                   /* [m_pitch] time domain hamming window      */
    COMP          W[FFT_ENC];	           /* DFT of w[]                                */
    float        *Pn;	                   /* [2*n_samp] trapezoidal synthesis window   */
    float        *bpf_buf;                 /* buffer for band pass filter               */
    float        *Sn;                      /* [m_pitch] input speech                    */
    float         hpf_states[2];           /* high pass filter states                   */
    void         *nlp;                     /* pitch predictor states                    */
    int           gray;                    /* non-zero for gray encoding                */

//    codec2_fftr_cfg  fftr_inv_cfg;         /* inverse FFT config                        */
    float        *Sn_;	                   /* [2*n_samp] synthesised output speech      */
    float         ex_phase;                /* excitation model phase track              */
    float         bg_est;                  /* background noise estimate for post filter */
    float         prev_f0_enc;             /* previous frame's f0    estimate           */
    MODEL         prev_model_dec;          /* previous frame's model parameters         */
    float         prev_lsps_dec[LPC_ORD];  /* previous frame's LSPs                     */
    float         prev_e_dec;              /* previous frame's LPC energy               */

    int           lpc_pf;                  /* LPC post filter on                        */
    int           bass_boost;              /* LPC post filter bass boost                */
    float         beta;                    /* LPC post filter parameters                */
    float         gamma;

    float         xq_enc[2];               /* joint pitch and energy VQ states          */
    float         xq_dec[2];

    int           smoothing;               /* enable smoothing for channels with errors */
    float        *softdec;                 /* optional soft decn bits from demod        */

/*
    // newamp1 states

    float          rate_K_sample_freqs_kHz[NEWAMP1_K];
    float          prev_rate_K_vec_[NEWAMP1_K];
    float          Wo_left;
    int            voicing_left;
    codec2_fft_cfg phase_fft_fwd_cfg;
    codec2_fft_cfg phase_fft_inv_cfg;      
    
    //newamp2 states (also uses newamp1 states )
    float 			energy_prev ;
    float          n2_rate_K_sample_freqs_kHz[NEWAMP2_K];
    float          n2_prev_rate_K_vec_[NEWAMP2_K];
    float          n2_pwb_rate_K_sample_freqs_kHz[NEWAMP2_16K_K];
    float          n2_pwb_prev_rate_K_vec_[NEWAMP2_16K_K];

*/

    /* used to dump features for deep learning experiments */
    FILE *flspEWov;
};

#endif
