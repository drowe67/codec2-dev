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
    int slot_local_frame_offset;    /* Where the RX frame starts, in samples, from the perspective of the modem */
    struct TDMA_SLOT * next_slot;   /* Next slot in a linked list of slots */

};

/* TDMA modem */
struct TDMA_MODEM {
    struct FSK * fsk_pilot;         /* Pilot modem */
    enum tdma_state state;          /* Current state of modem */
    struct TDMA_SLOT * slots;       /* Linked list of slot structs */

    int total_slot_count;

};

#endif
