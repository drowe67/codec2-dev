/*---------------------------------------------------------------------------*\

  FILE........: ofdm_mem.c
  AUTHOR......: Don Reid
  DATE CREATED: 11 June 2018

  Prints out the memory used by the OFDM modem states.  Used to optimise
  memory use for thw STM32F4 port.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2018 Don Reid

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

#include "ofdm_internal.h"
#include "ofdm_bpf_coeff.h"

extern float pilot_coeff[];

int main(int argc, char *argv[])
{
    struct OFDM *ofdm;

    int used = 0;

    printf("struct OFDM.................: %ld\n", sizeof(struct OFDM));
    printf("config......................: %ld\n", sizeof(ofdm->config));
    used +=                                       sizeof(ofdm->config);
    printf("pilot_samples...............: %ld\n", sizeof(ofdm->pilot_samples));
    used +=                                       sizeof(ofdm->pilot_samples);
    printf("W...........................: %ld\n", sizeof(ofdm->W));
    used +=                                       sizeof(ofdm->W);
    printf("rxbuf.......................: %ld\n", sizeof(ofdm->rxbuf));
    used +=                                       sizeof(ofdm->rxbuf);
    printf("pilots......................: %ld\n", sizeof(ofdm->pilots));
    used +=                                       sizeof(ofdm->pilots);
    printf("rx_sym......................: %ld\n", sizeof(ofdm->rx_sym));
    used +=                                       sizeof(ofdm->rx_sym);
    printf("rx_np.......................: %ld\n", sizeof(ofdm->rx_np));
    used +=                                       sizeof(ofdm->rx_np);
    printf("rx_amp......................: %ld\n", sizeof(ofdm->rx_amp));
    used +=                                       sizeof(ofdm->rx_amp);
    printf("aphase_est_pilot_log........: %ld\n", sizeof(ofdm->aphase_est_pilot_log));
    used +=                                       sizeof(ofdm->aphase_est_pilot_log);
    printf("tx_uw.......................: %ld\n", sizeof(ofdm->tx_uw));
    used +=                                       sizeof(ofdm->tx_uw);
    printf("sync_state..................: %ld\n", sizeof(ofdm->sync_state));
    used +=                                       sizeof(ofdm->sync_state);
    printf("last_sync_state.............: %ld\n", sizeof(ofdm->last_sync_state));
    used +=                                       sizeof(ofdm->last_sync_state);
    printf("sync_state_interleaver......: %ld\n", sizeof(ofdm->sync_state_interleaver));
    used +=                                       sizeof(ofdm->sync_state_interleaver);
    printf("last_sync_state_interleaver.: %ld\n", sizeof(ofdm->last_sync_state_interleaver));
    used +=                                       sizeof(ofdm->last_sync_state_interleaver);

    // add in non-array sizes
    int single = 0;
    single +=  8 * sizeof(int);
    single += 13 * sizeof(float);
    single +=  1 * sizeof(complex float);
    single +=  1 * sizeof(float *);
    single +=  4 * sizeof(bool);
    printf("single values...............: %d\n",  single);
    used +=                                       single;

    printf("alignment...................: %ld\n", (sizeof(struct OFDM) - used));


    printf("--allocated separately--\n");
    printf("tx_bpf_buf..................: %ld\n", (sizeof(complex float)*(OFDM_BPF_N+OFDM_SAMPLESPERFRAME)));


    return 0;
}

