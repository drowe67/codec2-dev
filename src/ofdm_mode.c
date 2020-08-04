/*---------------------------------------------------------------------------*\

  FILE........: ofdm_mode.c
  AUTHORS.....: David Rowe & Steve Sampson
  DATE CREATED: July 2020

  Mode specific configuration for OFDM modem.

\*---------------------------------------------------------------------------*/

#include <assert.h>
#include "comp.h"
#include "ofdm_internal.h"
#include "ofdm_mode.h"

void ofdm_init_mode(char mode[], struct OFDM_CONFIG *config) {
    assert(mode != NULL);
    assert(config != NULL);
    
    /* Fill in default values - 700D */

    config->nc = 17;                            /* Number of carriers */
    config->np = 1;
    config->ns = 8;                             /* Number of Symbols per modem frame */
    config->ts = 0.018f;
    config->tcp = .002f;                        /* Cyclic Prefix duration */
    config->tx_centre = 1500.0f;                /* TX Carrier Frequency */
    config->rx_centre = 1500.0f;                /* RX Carrier Frequency */
    config->fs = 8000.0f;                       /* Sample rate */
    config->txtbits = 4;
    config->bps = 2;                            /* Bits per Symbol */
    config->nuwbits = 5 * config->bps;          /* default is 5 symbols of Unique Word bits */
    config->bad_uw_errors = 3;
    config->ftwindowwidth = 11;
    config->timing_mx_thresh = 0.30f;
    config->data_mode = 0;
    config->codename = "HRA_112_112";
    
    if (strcmp(mode,"700D") == 0) {   
    } else if (strcmp(mode,"2020") == 0) {
         config->ts = 0.0205;  config->nc = 31; config->codename = "HRAb_396_504";
    } else if (strcmp(mode,"qam16") == 0) {
        config->ns=5; config->np=5; config->tcp = 0.004; config->ts = 0.016; config->nc = 33;
        config->bps=4; config->txtbits = 0; config->nuwbits = 15*4; config->bad_uw_errors = 5;
        config->ftwindowwidth = 32; config->data_mode = 1;
    } else if (strcmp(mode,"datac1") == 0) {
        config->ns=5; config->np=18; config->tcp = 0.006; config->ts = 0.016; config-> nc = 18;
        config->txtbits = 0; config->nuwbits = 12; config->bad_uw_errors = 2;
        config->ftwindowwidth = 32; config->data_mode = 1; config->codename = "H2064_516_sparse";
    } else if (strcmp(mode,"datac2") == 0) {
        config->ns=5; config->np=36; config->tcp = 0.006; config->ts = 0.016; config->nc = 9;
        config->txtbits = 0; config->nuwbits = 12; config->bad_uw_errors = 1;
        config->ftwindowwidth = 32; config->data_mode = 1; config->codename = "H2064_516_sparse";
    } else if (strcmp(mode,"datac3") == 0) {
        config->ns=5; config->np=11; config->tcp = 0.006; config->ts = 0.016; config->nc = 9;
        config->txtbits = 0; config->nuwbits = 24; config->bad_uw_errors = 1; /* TODO 5 */
        config->ftwindowwidth = 32; config->timing_mx_thresh = 0.30; config->data_mode = 1;
        config->codename = "H_256_768_22";
        /* TODO custom UW */
        //uint8_t uw[] = {1,1,0,0, 1,0,1,0,  1,1,1,1, 0,0,0,0, 1,1,1,1, 0,0,0,0};
        //memcpy(config->tx_uw, uw, config->nuwbits);
    }
    else {
        assert(0);
    }
    config->rs=1.0f/config->ts;
}

