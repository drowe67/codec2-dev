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

/*
 TODO: Lock sync point while transmitting in non-master mode and no master
 slot is detected
*/

/*
 TODO: Further testing and refinement with julia test suite
*/

#ifndef __CODEC_2_TDMA_H
#define __CODEC_2_TDMA_H

#include "fsk.h"
#include <stdint.h>
#include <stdbool.h>
#include "comp_prim.h"


#define TDMA_FRAME_A 3   /* 4800T Frame */

/* TODO: Replace these types with their full names */
/* I'm just feeling lazy right now */
typedef uint64_t u64;
typedef int64_t  i64;
typedef uint32_t u32;
typedef int32_t  i32;
typedef uint8_t  u8;
typedef float    f32;

typedef struct FSK fsk_t;
typedef struct TDMA_SLOT slot_t;
typedef struct TDMA_MODEM tdma_t;

/* The state for an individual slot */
enum slot_state {
    rx_no_sync = 0,         /* Not synched */
    rx_sync = 1,            /* Sunk */
    tx_client = 2,          /* TX but timed from a different master */
    tx_master = 3           /* TX in master mode */
};

/* The state of the entire TDMA modem */
enum tdma_state {
    no_sync = 0,            /* No sync */
    //pilot_sync = 1,         /* Pilot modem has gotten sync, but slots haven't*/
    slot_sync = 2,          /* One or more slots are sunk */
    master_sync = 3,        /* This modem is the TDMA master */
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
    fsk_t * fsk;                    /* The FSK modem for this slot */
    enum slot_state state;          /* Current local slot state */
    i32 slot_local_frame_offset;    /* Where the RX frame starts, in samples, from the perspective of the modem */
    u32 bad_uw_count;               /* How many bad UWs have we gotten since synchronized */
    i32 master_count;               /* How likely is this frame to be a synchronization master */
    slot_t * next_slot;             /* Next slot in a linked list of slots */
    bool single_tx;                 /* Are we TXing a single frame? */
};

/* Structure for tracking basic TDMA modem config */
struct TDMA_MODE_SETTINGS {
    u32 sym_rate;               /* Modem symbol rate */
    u32 fsk_m;                  /* Number of modem tones */
    u32 samp_rate;              /* Modem samplerate */
    u32 slot_size;              /* Number of symbols per slot, including quiet padding time */
    u32 frame_size;             /* Number of symbols per frame, not inclduing quiet padding */
    u32 n_slots;                /* Number of TDMA slots */
    u32 frame_type;             /* Frame type number for framer/deframer */
    u32 uw_len;                 /* Length of unique word */
    u32 pilot_sync_tol;         /* UW errors allowed for a valid pilot sync */
    u32 first_sync_tol;         /* UW errors allowed for a valid first frame sync */
    u32 frame_sync_tol;         /* UW errors allowed to maintain a frame sync */
    u32 frame_sync_baduw_tol;   /* How many bad UWs before calling a slot unsynced */
    i32 mastersat_max;          /* Maximum count for master detection counter */
    i32 mastersat_min;          /* Minimum count before frame considered 'master' */
    i32 loss_of_sync_frames;    /* How many bad frames before going from 'sync' to 'no_sync' for entire modem */
};

/* Declaration of basic 4800bps freedv tdma mode, defined in tdma.h */
//struct TDMA_MODE_SETTINGS FREEDV_4800T;

#define FREEDV_4800T {2400,4,48000,48,44,2,TDMA_FRAME_A,16,2,2,2,2,6,3,5};

/* Callback typedef that just returns the bits of the frame */
/* TODO: write this a bit better */
typedef void (*tdma_cb_rx_frame)(u8* frame_bits,u32 slot_i, slot_t * slot, tdma_t * tdma,u8 uw_type, void * cb_data);

/* Callback typedef when TDMA is ready to schedule a new frame */
/* Returns 1 if a frame is supplied, 0 if not */
/* If no frame supplied, slot is changed out of TX mode */
typedef int (*tdma_cb_tx_frame)(u8* frame_bits,u32 slot_i, slot_t * slot, tdma_t * tdma,u8 * uw_type, void * cb_data);

/* Callback to the radio front end to schedule a burst of TX samples */
typedef int (*tdma_cb_tx_burst)(tdma_t * tdma,COMP* samples, size_t n_samples,i64 timestamp,void * cb_data);

/* TDMA modem */
struct TDMA_MODEM {
    fsk_t * fsk_pilot;              /* Pilot modem */
    enum tdma_state state;          /* Current state of modem */
    slot_t * slots;                 /* Linked list of slot structs */
    struct TDMA_MODE_SETTINGS settings; /* Basic TDMA config parameters */
    COMP * sample_buffer;           /* Buffer of incoming samples */
    i32 sample_sync_offset;     /* Offset into the sample buffer where slot 0 starts */
    int64_t timestamp;             /* Timestamp of oldest sample in samp buffer */
    int64_t loop_delay;             /* Static offset applied to timestamp when scheduling 
                                        TX frames to account for delays in DSP and radio hardware */
    uint32_t tx_multislot_delay;    /* How many full slot periods in the future to delay TX burst scheduling */
    uint32_t slot_cur;              /* Current slot coming in */
    uint32_t sync_misses;           /* How many slots have been missed during this sync period */
    tdma_cb_rx_frame rx_callback;
    tdma_cb_tx_frame tx_callback;
    tdma_cb_tx_burst tx_burst_callback;
    void * rx_cb_data;
    void * tx_cb_data;
    void * tx_burst_cb_data;
    bool ignore_rx_on_tx;           /* Don't try and demod samples from a frame in a slot marked as TX */

    size_t master_bit_pos;          /* Where in the frame can we find the master indicator bit? */
    uint8_t uw_types;               /* How many different UWs does this framing format use? pulled from frame_type */
    uint8_t ** uw_list;             /* Pointer to list of valid UWs */
    

};


/* Allocate and setup a new TDMA modem */
tdma_t * tdma_create(struct TDMA_MODE_SETTINGS mode);

/* Tear down and free a TDMA modem */
void tdma_destroy(tdma_t * tdma);

/* Get the number of samples expected by RX for the next cycle */
u32 tdma_get_N(tdma_t * tdma);

/**
 Put 1 slot's worth of samples into the TDMA modem
 TODO: I'm still not entirely sure of what I want the semantics of this to look like
*/
void tdma_rx(tdma_t * tdma, COMP * samps,u64 timestamp);

/* Set the RX callback function */
void tdma_set_rx_cb(tdma_t * tdma,tdma_cb_rx_frame rx_callback,void * cb_data);

void tdma_set_tx_cb(tdma_t * tdma,tdma_cb_tx_frame tx_callback,void * cb_data);

void tdma_set_tx_burst_cb(tdma_t * tdma,tdma_cb_tx_burst tx_burst_callback, void * cb_data);

/* Set up TDMA to schedule the transmission of a single frame. The frame itself will be 
    passed in through the tx_frame callback
*/
void tdma_single_frame_tx(tdma_t * tdma, int slot_idx);

/* Start transmission of a bunch of frames on a particular slot
*/
void tdma_start_tx(tdma_t * tdma, int slot_idx);

/* Stop ongoing transmission of a bunch of frames for some slot
*/
void tdma_stop_tx(tdma_t * tdma, int slot_idx);

size_t tdma_nin(tdma_t * tdma);

size_t tdma_nout(tdma_t * tdma);

/* Convience function to look up a slot from it's index number */
slot_t * tdma_get_slot(tdma_t * tdma, u32 slot_idx);


#endif
