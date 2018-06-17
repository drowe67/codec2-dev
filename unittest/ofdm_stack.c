#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <unistd.h>

#include "ofdm_internal.h"
#include "codec2_ofdm.h"
#include "test_bits_ofdm.h"
#include "mpdecode_core.h"

#define DATA_BITSPERFRAME  (OFDM_BITSPERFRAME - (OFDM_NUWBITS + OFDM_NTXTBITS))
#define RX_OFFSET (OFDM_NUWBITS+OFDM_NTXTBITS)

#define MAX_ERRORS            32

struct OFDM    *ofdm;

// Top level data that is not really part of modem.
// Putting these here keeps them off the stack of the run_modem function.
int             tx_bits [DATA_BITSPERFRAME];        // Input to modulator
COMP            tx_rx   [OFDM_SAMPLESPERFRAME];     // Modulator to Demodulator
int             rx_bits [OFDM_BITSPERFRAME];        // Output of demodulator

// Forwards
void run_modem();
void dummy_code();

/////////////////////////////////////////////////////////////
///  MAIN()
int main(int argc, char *argv[]) {

    // Options
    int opt;
    int dummy = 0;  // flag to use dummy code
    int frames = 1; // how many frames
    int print = 0;  // flag to print all bits
    while ((opt = getopt(argc, argv, "df:p")) != -1) {
        switch (opt) {
            case 'd':
                dummy = 1;
                break;
            case 'f':
                frames = atoi(optarg);
                break;
            case 'p':
                print = 1;
                break;
            default:
                fprintf(stderr, "Usage: %s [-e] [-f <frames>] [-p]\n", argv[0]);
            }
        }

    int f,i; 

    for (f=0; f<frames; f++) {

        ////////
        // Prep inputs
        for(i=0; i<DATA_BITSPERFRAME; i++) {
            tx_bits[i] = payload_data_bits[(i % (sizeof(payload_data_bits)/sizeof(payload_data_bits[0])))];
        }

        ////////
        // Modem (or dummy)
        if (dummy) { dummy_code(); }
        else { run_modem(); }

        ////////
        // Compare results (or print)
        int errors = 0;
        if (print) {
            for(i=0; i<DATA_BITSPERFRAME; i++) {
                fprintf(stderr, "bit %3d: tx = %1d, rx = %1d",
                    i, tx_bits[i], rx_bits[i+RX_OFFSET]);
                if (tx_bits[i] != rx_bits[i+RX_OFFSET]) {
                    fprintf(stderr, " Error");
                    errors ++;
                }
                fprintf(stderr, "\n");
        }
        } else {
            for(i=0; i<DATA_BITSPERFRAME; i++) {
                if (tx_bits[i] != rx_bits[i+RX_OFFSET]) {
                    if (errors < MAX_ERRORS) {
                        fprintf(stderr, "Error in bit %3d: tx = %1d, rx = %1d\n",
                            i, tx_bits[i], rx_bits[i+RX_OFFSET]);
                    }
                    errors ++;
                }
            }
        }
        fprintf(stderr, "%d Errors\n", errors);

    } // for (f<frames

}   // end main()


//////////////////////////////////
void run_modem() {

    int             mod_bits[OFDM_SAMPLESPERFRAME];

    int i, j;

    ofdm = ofdm_create(NULL);

    ///////////
    // Mod
    ///////////

    for(i=0; i<OFDM_NUWBITS; i++) {
        mod_bits[i] = ofdm->tx_uw[i];
    }
    for(i=OFDM_NUWBITS; i<OFDM_NUWBITS+OFDM_NTXTBITS; i++) {
        mod_bits[i] = 0;
    }       

    for(j=0, i=OFDM_NUWBITS+OFDM_NTXTBITS; j<DATA_BITSPERFRAME; i++,j++) {
        mod_bits[i] = tx_bits[j];
    }
    for(j=0; j<DATA_BITSPERFRAME; i++,j++) {
        mod_bits[i] = tx_bits[j];
    }

    ofdm_mod(ofdm, (COMP*)tx_rx, mod_bits);

    ///////////
    // DeMod
    ///////////

    int  Nsam = OFDM_SAMPLESPERFRAME;
    int  prx = 0;
    int  nin =  OFDM_SAMPLESPERFRAME + 2*(OFDM_M+OFDM_NCP);

    int  lnew;
    COMP rxbuf_in[OFDM_MAX_SAMPLESPERFRAME];

    for (i=0; i<OFDM_SAMPLESPERFRAME ; i++,prx++) {
        ofdm->rxbuf[OFDM_RXBUF-nin+i] = tx_rx[prx].real + I*tx_rx[prx].imag;
    }
    for (i=OFDM_SAMPLESPERFRAME ; i<nin; i++) {
        ofdm->rxbuf[OFDM_RXBUF-nin+i] = 0 + I*0;
    }
    
    /* disable estimators for initial testing */
    ofdm_set_verbose(ofdm, false);
    ofdm_set_timing_enable(ofdm, true);
    ofdm_set_foff_est_enable(ofdm, true);
    ofdm_set_phase_est_enable(ofdm, true);

    ofdm->mean_amp = 1.0;

    nin = ofdm_get_nin(ofdm);

    /* Insert samples at end of buffer, set to zero if no samples
       available to disable phase estimation on future pilots on
       last frame of simulation. */

    if ((Nsam-prx) < nin) {
        lnew = Nsam-prx;
    } else {
        lnew = nin;
    }
    for(i=0; i<nin; i++) {
        rxbuf_in[i].real = 0.0;
        rxbuf_in[i].imag = 0.0;
    }

    if (lnew) {
        for(i=0; i<lnew; i++, prx++) {
            rxbuf_in[i] = tx_rx[prx];
        }
    }

    ofdm_demod(ofdm, rx_bits, rxbuf_in);
    

}   // end run_modem()


//////////////////////////////////
void dummy_code() {

    int i;

    for(i=0; i<DATA_BITSPERFRAME; i++) {
        rx_bits[i] = tx_bits[i];
    }

}   // end dummy_code()

/* vi:set ts=4 sts=4 et: */
