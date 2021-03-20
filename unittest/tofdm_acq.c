/*---------------------------------------------------------------------------*\

  FILE........: tofdm_acq.c
  AUTHORS.....: David Rowe
  DATE CREATED: Mar 2021

  Tests for the acquistion (sync) parts of the C version of the OFDM modem. 
  This program outputs a file of Octave vectors that are loaded and 
  automatically tested against the Octave version of the modem by the Octave 
  script tofdm_acq.m

\*---------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "ofdm_internal.h"
#include "codec2_ofdm.h"
#include "octave.h"

int main(int argc, char *argv[])
{
    struct OFDM *ofdm;
    struct OFDM_CONFIG ofdm_config;

    ofdm_init_mode("datac0", &ofdm_config);
    ofdm = ofdm_create(&ofdm_config);
    ofdm->data_mode = "burst";
    ofdm->verbose = 2;
    assert(ofdm != NULL);
    
    int nin_frame = ofdm_get_nin(ofdm);
    FILE *fin = fopen(argv[1],"rb"); assert(fin != NULL);
    short rx_scaled[ofdm_get_max_samples_per_frame(ofdm)];
    int f = 0;
    
    while (fread(rx_scaled, sizeof (short), nin_frame, fin) == nin_frame) {
        fprintf(stderr, "%3d ", f++);
        ofdm_sync_search_shorts(ofdm, rx_scaled, ofdm->amp_scale / 2.0f);
        // this is modified when a valid pre-amble is found. Force fixed nin to stay in acq state
        ofdm->nin = nin_frame; 
    }
    fclose(fin);
       
    /*---------------------------------------------------------*\
               Dump logs to Octave file for evaluation
                      by tofdm_acq.m Octave script
    \*---------------------------------------------------------*/

    FILE *fout = fopen("tofdm_acq_out.txt","wt");
    assert(fout != NULL);
    fprintf(fout, "# Created by tofdm_acq.c\n");
    octave_save_complex(fout, "tx_preamble_c", (COMP*)ofdm->tx_preamble, 1, ofdm->samplesperframe, ofdm->samplesperframe);
    fclose(fout);

    ofdm_destroy(ofdm);

    return 0;
}
