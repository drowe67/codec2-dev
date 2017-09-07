/*---------------------------------------------------------------------------*\

  FILE........: tdma.h
  AUTHOR......: Brady O'Brien
  DATE CREATED: 18 September 2016

  Skeletion of the API for the TDMA FSK modem

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2016 Brady O'Brien

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

#ifndef __CODEC_2_TDMA_H
#define __CODEC_2_TDMA_H

#include "fsk.h"
#include "freedv_vhf_framing.h"
#include <stdint.h>
#include "comp_prim.h"


/* TODO: Replace these types with their full names */
/* I'm just feeling lazy right now */
typedef uint64_t u64;
typedef uint32_t u32;
typedef int32_t  i32;
typedef uint8_t  u8;
typedef float    f32;

//typedef void (*tdma_cb_rx_frame)()

/* The state for an individual slot */
enum slot_state {
    rx_no_sync,         /* Not synched */
    rx_sync,            /* Sunk */
    tx_client,          /* TX but timed from a different master */
    tx_master           /* TX in master mode */
};

/* The state of the entire TDMA modem */
enum tdma_state {
    no_sync,            /* No sync */
    pilot_sync,         /* Pilot modem has gotten sync, but slots haven't*/
    slot_sync,          /* One or more slots are sunk */
    master_sync,        /* This modem is the TDMA master */
};

/* TDMA frame type */
enum tdma_frame_type{
    frame_master,
    frame_client,
};

/* TDMA frame struct */
struct TDMA_FRAME {
    enum tdma_frame_type type;      /* Type of frame */
    int slot_idx;                   /* Index of slot from where frame was rx-ed */
    uint8_t frame_payload[];        /* Frame payload. TODO: figure out how to sling payloads around */
};

/* TDMA slot struct */

struct TDMA_SLOT {
    struct FSK * fsk;               /* The FSK modem for this slot */
    enum slot_state state;          /* Current local slot state */
    uint32_t slot_local_frame_offset;    /* Where the RX frame starts, in samples, from the perspective of the modem */
    struct TDMA_SLOT * next_slot;   /* Next slot in a linked list of slots */

};

/* Structure for tracking basic TDMA modem config */
struct TDMA_MODE_SETTINGS {
    uint32_t sym_rate;              /* Modem symbol rate */
    uint32_t fsk_m;                 /* Number of modem tones */
    uint32_t samp_rate;             /* Modem samplerate */
    uint32_t slot_size;             /* Number of symbols per slot, including quiet padding time */
    uint32_t frame_size;            /* Number of symbols per frame, not inclduing quiet padding */
    uint32_t n_slots;                /* Number of TDMA slots */
    uint32_t frame_type;             /* Frame type number for framer/deframer */
};

/* Declaration of basic 4800bps freedv tdma mode, defined in tdma.h */
//struct TDMA_MODE_SETTINGS FREEDV_4800T;

#define FREEDV_4800T {2400,4,48000,48,44,2,FREEDV_VHF_FRAME_AT}

/* TDMA modem */
struct TDMA_MODEM {
    struct FSK * fsk_pilot;         /* Pilot modem */
    enum tdma_state state;          /* Current state of modem */
    struct TDMA_SLOT * slots;       /* Linked list of slot structs */
    struct TDMA_MODE_SETTINGS settings; /* Basic TDMA config parameters */
    COMP * sample_buffer;          /* Buffer of incoming samples */
    size_t sample_sync_offset;      /* Offset into the sample buffer where slot 0 starts */
};

/* Allocate and setup a new TDMA modem */
struct TDMA_MODEM * tdma_create(struct TDMA_MODE_SETTINGS mode);

/* Tear down and free a TDMA modem */
void tdma_destroy(struct TDMA_MODEM * tdma);

/* Get the number of samples expected by RX for the next cycle */
u32 tdma_get_N(struct TDMA_MODEM * tdma);

/**
 Put 1 slot's worth of samples into the TDMA modem
 TODO: I'm still not entirely sure of what I want the semantics of this to look like
*/
void tdma_rx(struct TDMA_MODEM * tdma, COMP * samps,u64 timestamp);

/* Hideous debug function */
void tdma_print_stuff(struct TDMA_MODEM * tdma);

#endif
