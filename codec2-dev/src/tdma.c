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
#include <stdlib.h>


static const uint8_t TDMA_UW_V[] =    {0,1,1,0,0,1,1,1,
                                       1,0,1,0,1,1,0,1};

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
    COMP * samp_buffer = NULL;
    
    size_t i;

    assert( (Fs%Rs)==0 );
    assert( M==2 || M==4);

    /* allocate the modem */
    tdma = (struct TDMA_MODEM *) malloc(sizeof(struct TDMA_MODEM));
    if(tdma == NULL) goto cleanup_bad_alloc;

    /* Symbols over which pilot modem operates */
    u32 pilot_nsyms = slot_size/2;

    /* Set up pilot modem */
    struct FSK * pilot = fsk_create_hbr(Fs,Rs,P,M,Rs,Rs);
    if(pilot == NULL) goto cleanup_bad_alloc;
    fsk_enable_burst_mode(pilot,pilot_nsyms);
    tdma->fsk_pilot = pilot;
    
    tdma->settings = mode;
    tdma->state = no_sync;
    tdma->sample_sync_offset = 960;
    tdma->slot_cur = 0;

    /* Allocate buffer for incoming samples */
    /* TODO: We may only need a single slot's worth of samps -- look into this */
    samp_buffer = (COMP *) malloc(sizeof(COMP)*slot_size*Ts*(n_slots+1));
    if(samp_buffer == NULL) goto cleanup_bad_alloc;

    tdma->sample_buffer = samp_buffer;
    for(i=0; i<slot_size*Ts*n_slots; i++){
        tdma->sample_buffer[i].real = 0;
        tdma->sample_buffer[i].imag = 0;
    }

    struct TDMA_SLOT * slot;
    struct TDMA_SLOT * last_slot;
    struct FSK * slot_fsk;
    last_slot = NULL;
    for(i=0; i<n_slots; i++){
        slot = (struct TDMA_SLOT *) malloc(sizeof(struct TDMA_SLOT));
        if(slot == NULL) goto cleanup_bad_alloc;
        slot->fsk = NULL;
        tdma->slots = slot;
        slot->next_slot = last_slot;
        slot->slot_local_frame_offset = 0;
        slot->state = rx_no_sync;
        slot->bad_uw_count = 0;
        slot_fsk = fsk_create_hbr(Fs,Rs,P,M,Rs,Rs);
        
        if(slot_fsk == NULL) goto cleanup_bad_alloc;

        /* NOTE/TODO: demod one extra symbol (probably of zeros) off the end of the frame */
        fsk_enable_burst_mode(slot_fsk, slot_size+1);
        
        slot->fsk = slot_fsk;
        last_slot = slot;
    }

    return tdma;

    /* Clean up after a failed malloc */
    /* TODO: Run with valgrind/asan, make sure I'm getting everything */
    cleanup_bad_alloc:
    fprintf(stderr,"Cleaning up\n");
    if(tdma == NULL) return NULL;

    struct TDMA_SLOT * cleanup_slot = tdma->slots;
    struct TDMA_SLOT * cleanup_slot_next;
    while(cleanup_slot != NULL){
        cleanup_slot_next = cleanup_slot->next_slot;
        if(cleanup_slot->fsk != NULL) fsk_destroy(cleanup_slot->fsk);
        if(cleanup_slot != NULL) free(cleanup_slot);
        cleanup_slot = cleanup_slot_next;
    }
    if(pilot != NULL) fsk_destroy(pilot);
    if(samp_buffer != NULL) free(samp_buffer);
    free(tdma);
    return NULL;
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

void tdma_destroy(tdma_t * tdma){
    /* TODO: Free slot modems (need to create them first) */
    fsk_destroy(tdma->fsk_pilot);
    free(tdma->sample_buffer);
    free(tdma);
}

u32 tdma_get_N(tdma_t * tdma){
    u32 slot_size = tdma->settings.slot_size;
    u32 Fs = tdma->settings.samp_rate;
    u32 Rs = tdma->settings.sym_rate;
    return slot_size * (Fs/Rs);
}

/* Convience function to look up a slot from it's index number */
static slot_t * tdma_get_slot(tdma_t * tdma, u32 slot_idx){
    /* Don't try and index beyond the end */
    if(slot_idx >= tdma->settings.n_slots) return NULL;

    size_t i;
    slot_t * cur = tdma->slots;
    for(i = 0; i < slot_idx; i++){
        /* Don't break */
        if(cur == NULL) return NULL;
        cur = cur->next_slot;
    }
    return cur;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

int tdma_demod_end_slot(tdma_t * tdma,u32 slot_idx, u8 * bit_buf){

    struct TDMA_MODE_SETTINGS mode = tdma->settings;
    u32 Rs = mode.sym_rate;
    u32 Fs = mode.samp_rate;
    u32 slot_size = mode.slot_size;
    u32 frame_size = mode.frame_size;
    u32 n_slots = mode.n_slots;
    u32 M = mode.fsk_m;
    u32 Ts = Fs/Rs;
    u32 bits_per_sym = M==2?1:2;
    u32 slot_samps = slot_size*Ts;
    u32 buffer_samps = slot_samps*n_slots;
    size_t nbits = (slot_size+1)*bits_per_sym;
    u32 frame_bits = frame_size*bits_per_sym;

    /* Samples that belong to this frame */
    COMP frame_samps[(slot_size+1)*Ts];
    COMP * sample_buffer = tdma->sample_buffer;
    slot_t * slot = tdma_get_slot(tdma,slot_idx);
    struct FSK * fsk = slot->fsk;

    /* Zero out tail end of bit buffer so we can get last symbol out of demod */
    /* TODO: This is a hack. Look into better burst mode support in FSK */
    size_t i;
    for(i = slot_samps; i< (slot_size+1)*Ts; i++){
        frame_samps[i].real = 0;
        frame_samps[i].imag = 0;
    }

    /* Pull out the frame and demod */
    memcpy(&frame_samps[0],&sample_buffer[tdma->sample_sync_offset],slot_samps*sizeof(COMP));

    /* Demodulate the frame */
    fsk_demod(fsk,bit_buf,frame_samps);

    i32 delta,off;
    off = fvhff_search_uw(bit_buf,nbits,TDMA_UW_V,16,&delta);
    i32 f_start = off- (frame_bits-16)/2;
    int f_valid = 0; /* Flag indicating wether or not we've found a UW;

    /* Check frame tolerance and sync state*/
    if(slot->state == rx_sync){
        f_valid = delta <= tdma->settings.frame_sync_tol;
    }else if(slot->state == rx_no_sync){
        f_valid = delta <= tdma->settings.first_sync_tol;
    }

    /* Calculate offset (in samps) from start of frame */
    /* Note: FSK outputs one symbol from the last batch, so we have to account for that */
    i32 target_frame_offset = ((slot_size-frame_size)/2)*Ts;
    i32 frame_offset = ((f_start-bits_per_sym)*(Ts/bits_per_sym)) - target_frame_offset;
    if(f_valid)
        slot->slot_local_frame_offset = frame_offset;
    
    fprintf(stderr,"slot: %d fstart:%d offset: %d delta: %d f1:%.3f EbN0:%f\n",slot_idx,frame_offset,off,delta,fsk->f_est[0],fsk->EbNodB);
    for(i=0; i<nbits; i++){
        if((i>off && i<=off+16) || i==f_start || i==(f_start+frame_bits-1)){
            if(bit_buf[i])  fprintf(stderr,"1̲");
            else            fprintf(stderr,"0̲");
        } else fprintf(stderr,"%d",bit_buf[i]);
    }
    fprintf(stderr,"\n");
    return f_valid;

}
/* We got a new slot's worth of samples. Run the slot modem and try to get slot sync */
/* This will probably also work for the slot_sync state */
void tdma_rx_pilot_sync(tdma_t * tdma){
    struct TDMA_MODE_SETTINGS mode = tdma->settings;
    u32 Rs = mode.sym_rate;
    u32 Fs = mode.samp_rate;
    u32 slot_size = mode.slot_size;
    u32 frame_size = mode.frame_size;
    u32 n_slots = mode.n_slots;
    u32 M = mode.fsk_m;
    u32 Ts = Fs/Rs;
    u32 slot_samps = slot_size*Ts;
    u32 bits_per_sym = M==2?1:2;
    slot_t * slot = tdma_get_slot(tdma,tdma->slot_cur);
    struct FSK * fsk = slot->fsk;
    size_t nbits = (slot_size+1)*bits_per_sym;
    size_t slot_offset = tdma->sample_sync_offset;
    u8 bit_buf[nbits];
    COMP * sample_buffer = tdma->sample_buffer;
    COMP frame_samps[(slot_size+1)*Ts];
    u32 frame_bits = frame_size*bits_per_sym;

    /* Compensate for frame timing offset sliding towards end of buffer */
    /* As in, don't demodulate this slot since we're getting ahead of ourselves */
    if( (slot_offset+slot_samps) > (slot_samps*(n_slots+1)-(slot_samps/4)) ){
        /* Move slot offset back by 1 slot and don't increment slot index. We'll just handle this one on the next batch of samps */
        tdma->sample_sync_offset -= slot_samps;
        fprintf(stderr,"Skipping\n");
        return;
    }

    /* Demod a slot in the sample buffer */
    //int f_valid = tdma_demod_end_slot(tdma,tdma->slot_cur,bit_buf);

    /* Samples that belong to this frame */

        /* Zero out tail end of bit buffer so we can get last symbol out of demod */
    /* TODO: This is a hack. Look into better burst mode support in FSK */
    size_t i;
    for(i = slot_samps; i< (slot_size+1)*Ts; i++){
        frame_samps[i].real = 0;
        frame_samps[i].imag = 0;
    }

    /* Pull out the frame and demod */
    memcpy(&frame_samps[0],&sample_buffer[tdma->sample_sync_offset],slot_samps*sizeof(COMP));

    /* Demodulate the frame */
    fsk_demod(fsk,bit_buf,frame_samps);

    size_t delta,off;
    off = fvhff_search_uw(bit_buf,nbits,TDMA_UW_V,16,&delta);
    i32 f_start = off- (frame_bits-16)/2;
    int f_valid = 0; /* Flag indicating wether or not we've found a UW;

    /* Check frame tolerance and sync state*/
    if(slot->state == rx_sync){
        f_valid = delta <= tdma->settings.frame_sync_tol;
    }else if(slot->state == rx_no_sync){
        f_valid = delta <= tdma->settings.first_sync_tol;
    }

    /* Calculate offset (in samps) from start of frame */
    /* Note: FSK outputs one symbol from the last batch, so we have to account for that */
    i32 target_frame_offset = ((slot_size-frame_size)/2)*Ts;
    i32 frame_offset = ((f_start-bits_per_sym)*(Ts/bits_per_sym)) - target_frame_offset;
    if(f_valid)
        slot->slot_local_frame_offset = frame_offset;
    

    i32 single_slot_offset = slot->slot_local_frame_offset;
    /* Flag indicating wether or not we should call the callback */
    int do_frame_found_call = 0;   

    /* TODO: deal with re-demod for frames just outside of buffer */
    /* Do single slot state machine */
    if( slot->state == rx_sync){
        fprintf(stderr,"Slot %d: sunk\n",tdma->slot_cur);
        if(!f_valid){   /* on bad UW, increment bad uw count and possibly unsync */
            slot->bad_uw_count++;
            if(slot->bad_uw_count > tdma->settings.frame_sync_baduw_tol){
                slot->state = rx_no_sync;
            }else{
                do_frame_found_call = 1;
            }
        }else{ /* Good UW found */
            slot->bad_uw_count = 0;
            do_frame_found_call = 1;
        }
    }else if(slot->state == rx_no_sync){
        fprintf(stderr,"Slot %d: no sync\n",tdma->slot_cur);
        if(f_valid){
            slot->state = rx_sync;
            do_frame_found_call;
        }
    }

    /* Unicode underline for pretty printing */
    char underline[] = {0xCC,0xB2,0x00};

    fprintf(stderr,"slot: %d fstart:%d offset: %d delta: %d f1:%.3f EbN0:%f\n",tdma->slot_cur,frame_offset,off,delta,fsk->f_est[0],fsk->EbNodB);
    for(i=0; i<nbits; i++){
        fprintf(stderr,"%d",bit_buf[i]);
        if((i>off && i<=off+16) || i==f_start || i==(f_start+frame_bits-1)){
            fprintf(stderr,underline);
        }
    }
    fprintf(stderr,"\n");

    /* Update slot offset to compensate for frame centering */
    u32 offset_total = 0;
    u32 offset_slots = 0;
    for( i=0; i<n_slots; i++){
        /* Only check offset from valid frames */
        if(tdma_get_slot(tdma,i)->state == rx_sync){
            i32 local_offset = tdma_get_slot(tdma,i)->slot_local_frame_offset;
            /* Filter out extreme spurious timing offsets */
            if(abs(local_offset)<(slot_samps/4)){
                fprintf(stderr,"Local offset: %d\n",local_offset);
                offset_total+=local_offset;
                offset_slots++;
            }
        }
    }
    offset_total = offset_slots>0 ? offset_total/offset_slots:0;
    fprintf(stderr,"Total Offset:%d\n",offset_total);
    fprintf(stderr,"Slot offset: %d of %d\n",tdma->sample_sync_offset,slot_samps*n_slots);
    //tdma->sample_sync_offset -= (offset_total/4);

    tdma->slot_cur++;
    if(tdma->slot_cur >= n_slots)
        tdma->slot_cur = 0;

    /* Compensate for frame timing offset sliding towards the start of buffer */
    /* Do so by running this again, demodulating two slots, and moving the slot offset one slot forward */
    if( slot_offset < (slot_samps/4) ){
        tdma->sample_sync_offset = tdma->sample_sync_offset + slot_samps;
        tdma_rx_pilot_sync(tdma);
        fprintf(stderr,"Recursing\n");
    }
        
}

void tdma_rx_no_sync(tdma_t * tdma, COMP * samps, u64 timestamp){
    struct TDMA_MODE_SETTINGS mode = tdma->settings;
    u32 Rs = mode.sym_rate;
    u32 Fs = mode.samp_rate;
    u32 slot_size = mode.slot_size;
    //u32 frame_size = mode.frame_size;
    u32 n_slots = mode.n_slots;
    u32 M = mode.fsk_m;
    u32 Ts = Fs/Rs;
    u32 bits_per_sym = M==2?1:2;

    //Number of bits per pilot modem chunk (half a slot)
    u32 n_pilot_bits = (slot_size/2)*bits_per_sym;
    //We look at a full slot for the UW
    u8 pilot_bits[n_pilot_bits];


    /*
    Pseudocode:
        copy into samp buffer 
            (? may want to let the downstream stuff do this; We may not even really need a local samp buffer)
            (could/probably should do this in tdma_rx)
        demod a half slot
        look for UW in slot-wide bit buffer
        if UW found
          change to pilot sync state
          set slot offset to match where a frame should be
          go to slot_sync_rx to try and decode found frame
        else  
          set up for next half slot, repeat 
            next half slot overlaps current slot      
   */
}

void tdma_rx(tdma_t * tdma, COMP * samps,u64 timestamp){

    COMP * sample_buffer = tdma->sample_buffer;
    struct TDMA_MODE_SETTINGS mode = tdma->settings;
    u32 Rs = mode.sym_rate;
    u32 Fs = mode.samp_rate;
    u32 slot_size = mode.slot_size;
    //u32 frame_size = mode.frame_size;
    u32 n_slots = mode.n_slots;
    u32 M = mode.fsk_m;
    u32 Ts = Fs/Rs;
    u32 bits_per_sym = M==2?1:2;
    u32 slot_samps = slot_size*Ts;

    /* Copy samples into the local buffer for some reason */
    /* Move the current samps in the buffer back by a slot or so */
    size_t move_samps = slot_samps*n_slots*sizeof(COMP);
    uintptr_t move_from = ((uintptr_t)sample_buffer) + n_slots*slot_samps*sizeof(COMP);
    uintptr_t move_to = (uintptr_t)sample_buffer; /* Don't really need this, but it's cleaner than doing it all in memmove */
    memmove((void*)move_to,(void*)move_from,move_samps);

    move_samps = slot_samps*sizeof(COMP);
    move_from = (uintptr_t)samps;
    move_to = ((uintptr_t)sample_buffer) + (n_slots)*slot_samps*sizeof(COMP);
    memcpy((void*)move_to,(void*)move_from,move_samps);

    /* Set the timestamp. Not sure if this makes sense */
    tdma->timestamp = timestamp - (slot_samps*(n_slots-1));

    /* Staate machine for TDMA modem */
    switch(tdma->state){
        case no_sync:
            tdma_rx_no_sync(tdma,samps,timestamp);
            break;
        case pilot_sync:
            tdma_rx_pilot_sync(tdma);
            break;
        case slot_sync:
            break;
        case master_sync:
            break;
        default:
            break;
    }
    tdma->state = pilot_sync;
}


void tdma_set_rx_cb(tdma_t * tdma,tdma_cb_rx_frame rx_callback){
    tdma->rx_callback = rx_callback;
}


#pragma GCC diagnostic pop
