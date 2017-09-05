/*---------------------------------------------------------------------------*\

  FILE........: tdma.c
  AUTHOR......: Brady O'Brien
  DATE CREATED: 18 September 2016

  Skeletion of the TDMA FSK modem

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2017 Brady O'Brien

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


#include "fsk.h"
#include "freedv_vhf_framing.h"
#include "tdma.h"
#include <stdint.h>
#include <assert.h>
#include <stdio.h>

struct TDMA_MODEM * tdma_create(struct TDMA_MODE_SETTINGS mode){
    struct TDMA_MODEM * tdma;
    
    u32 Rs = mode.sym_rate;
    u32 Fs = mode.samp_rate;
    u32 slot_size = mode.slot_size;
    //u32 frame_size = mode.frame_size;
    u32 n_slots = mode.n_slots;
    u32 M = mode.fsk_m;
    u32 P = Fs/Rs;
    u32 Ts = Fs/Rs;
    
    assert( (Fs%Rs)==0 );
    assert( M==2 || M==4);

    /* allocate the modem */
    tdma = (struct TDMA_MODEM *) malloc(sizeof(struct TDMA_MODEM));
    /* TODO: Malloc NULL checks */

    /* Symbols over which pilot modem operates */
    u32 pilot_nsyms = slot_size/2;

    /* Set up pilot modem */
    struct FSK * pilot = fsk_create_hbr(Fs,Rs,P,M,Rs,Rs);
    fsk_set_nsym(pilot,pilot_nsyms);
    tdma->fsk_pilot = pilot;
    
    tdma->settings = mode;
    tdma->state = no_sync;
    tdma->sample_sync_offset = 0;

    /* Allocate buffer for incoming samples */
    /* TODO: We may only need a single slot's worth of samps -- look into this */
    COMP * samp_buffer = (COMP *) malloc(sizeof(COMP)*slot_size*Ts*n_slots);

    tdma->sample_buffer = samp_buffer;

    /* TODO: Allocate slot modems. Get pilot detection working first */

    return tdma;

}

void tdma_print_stuff(struct TDMA_MODEM * tdma){
    printf("symrate: %d\n",tdma->settings.sym_rate);
    printf("fsk_m: %d\n",tdma->settings.fsk_m);
    printf("samprate: %d\n",tdma->settings.samp_rate);
    printf("slotsize: %d\n",tdma->settings.slot_size);
    printf("framesize: %d\n",tdma->settings.frame_size);
    printf("nslots: %d\n",tdma->settings.n_slots);
    printf("frametype: %d\n",tdma->settings.frame_type);
    printf("sync_offset: %ld\n",tdma->sample_sync_offset);
}

void tdma_destroy(struct TDMA_MODEM * tdma){
    /* TODO: Free slot modems (need to create them first) */
    fsk_destroy(tdma->fsk_pilot);
    free(tdma->sample_buffer);
    free(tdma);
}

u32 tdma_get_N(struct TDMA_MODEM * tdma){
    u32 slot_size = tdma->settings.slot_size;
    u32 Fs = tdma->settings.samp_rate;
    u32 Rs = tdma->settings.sym_rate;
    return slot_size * (Fs/Rs);
}

void tdma_rx(struct TDMA_MODEM * tdma, COMP * samps){

}