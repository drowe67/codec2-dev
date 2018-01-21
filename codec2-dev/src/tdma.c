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
#include "tdma.h"
#include <stdint.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

/* Easy handle to enable/disable a whole slew of debug printouts */
//#define VERY_DEBUG 1

static const uint8_t TDMA_UW_V[] =    {0,1,1,0,0,1,1,1,
                                       1,0,1,0,1,1,0,1};

static const uint8_t TDMA_UW_D[] =    {1,1,1,1,0,0,0,1,
                                       1,1,1,1,1,1,0,0};         
                                       
static const uint8_t * TDMA_UW_LIST_A[] = {&TDMA_UW_V,&TDMA_UW_D};

tdma_t * tdma_create(struct TDMA_MODE_SETTINGS mode){
    tdma_t * tdma;
    
    u32 Rs = mode.sym_rate;
    u32 Fs = mode.samp_rate;
    u32 slot_size = mode.slot_size;
    u32 n_slots = mode.n_slots;
    u32 M = mode.fsk_m;
    u32 P = Fs/Rs;
    u32 Ts = Fs/Rs;
    COMP * samp_buffer = NULL;
    
    size_t i;

    assert( (Fs%Rs)==0 );
    assert( M==2 || M==4);

    /* allocate the modem */
    tdma = (tdma_t *) malloc(sizeof(tdma_t));
    if(tdma == NULL) goto cleanup_bad_alloc;

    /* Symbols over which pilot modem operates */
    u32 pilot_nsyms = slot_size/2;

    /* Set up pilot modem */
    fsk_t * pilot = fsk_create_hbr(Fs,Rs,P,M,Rs,Rs);
    if(pilot == NULL) goto cleanup_bad_alloc;
    fsk_enable_burst_mode(pilot,pilot_nsyms);
    tdma->fsk_pilot = pilot;
    
    tdma->settings = mode;
    tdma->state = no_sync;
    tdma->sample_sync_offset = 960;
    tdma->slot_cur = 0;
    tdma->rx_callback = NULL;
    tdma->tx_callback = NULL;
    tdma->tx_burst_callback = NULL;
    tdma->ignore_rx_on_tx = true;
    tdma->sync_misses = 0;

    /* Set up the UWs we use for this mode. */
    if(mode.frame_type == TDMA_FRAME_A){
        tdma->uw_types = 2;
        tdma->uw_list = TDMA_UW_LIST_A;
        tdma->master_bit_pos = 35;
    }
    /* Allocate buffer for incoming samples */
    /* TODO: We may only need a single slot's worth of samps -- look into this */
    samp_buffer = (COMP *) malloc(sizeof(COMP)*slot_size*Ts*(n_slots+1));
    if(samp_buffer == NULL) goto cleanup_bad_alloc;

    tdma->sample_buffer = samp_buffer;
    for(i=0; i<slot_size*Ts*n_slots; i++){
        tdma->sample_buffer[i].real = 0;
        tdma->sample_buffer[i].imag = 0;
    }

    slot_t * slot;
    slot_t * last_slot;
    fsk_t * slot_fsk;
    last_slot = NULL;
    for(i=0; i<n_slots; i++){
        slot = (slot_t *) malloc(sizeof(slot_t));
        if(slot == NULL) goto cleanup_bad_alloc;
        slot->fsk = NULL;
        tdma->slots = slot;
        slot->next_slot = last_slot;
        slot->slot_local_frame_offset = 0;
        slot->state = rx_no_sync;
        slot->single_tx = true;
        slot->bad_uw_count = 0;
        slot->master_count = 0;
        slot_fsk = fsk_create_hbr(Fs,Rs,P,M,Rs,Rs);
        
        if(slot_fsk == NULL) goto cleanup_bad_alloc;
        fsk_enable_burst_mode(slot_fsk, slot_size+1);
        
        slot->fsk = slot_fsk;
        last_slot = slot;
    }

    return tdma;

    /* Clean up after a failed malloc */
    cleanup_bad_alloc:
    if(tdma == NULL) return NULL;

    slot_t * cleanup_slot = tdma->slots;
    slot_t * cleanup_slot_next;
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

void tdma_destroy(tdma_t * tdma){
    slot_t * slot = tdma->slots;
    slot_t * next_slot;
    while(slot != NULL){
        next_slot = slot->next_slot;
        fsk_destroy(slot->fsk);
        free(slot);
        slot = next_slot;
    }
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

/* Search for a complete UW in a buffer of bits */
size_t tdma_search_uw(tdma_t * tdma, u8 bits[], size_t nbits, size_t * delta_out, size_t * uw_type_out){
/*size_t fvhff_search_uw(const uint8_t bits[],size_t nbits,
    const uint8_t uw[],    size_t uw_len,
    size_t * delta_out,    size_t bits_per_sym){*/

    size_t uw_len = tdma->settings.uw_len;
    size_t bits_per_sym = (tdma->settings.fsk_m==2)?1:2;
    u8 ** uw_list = tdma->uw_list;
    size_t ibits,iuw;
    size_t delta_min = uw_len;
    size_t delta;
    size_t j;
    size_t uw_type_min;
    size_t uw_delta_min = uw_len;
    size_t offset_min,uw_offset_min;
    /* Check each UW */
    for(j = 0; j < tdma->uw_types; j++){
        /* Walk through buffer bits */
        u8 * uw = uw_list[j];
        size_t offset_min = 0;
        for(ibits = 0; ibits < nbits-uw_len; ibits+=bits_per_sym){
            delta = 0;
            for(iuw = 0; iuw < uw_len; iuw++){
                if(bits[ibits+iuw] != uw[iuw]) delta++;
            }
            if( delta < delta_min ){
                delta_min = delta;
                offset_min = ibits;
            }
        }
        if( delta_min < uw_delta_min ){
            uw_delta_min = delta_min;
            uw_offset_min = offset_min;
            uw_type_min = j;
        } 
    }
    if(delta_out != NULL) *delta_out = uw_delta_min;
    if(uw_type_out != NULL) *uw_type_out = uw_type_min;
    return uw_offset_min;
}


/* Convience function to look up a slot from it's index number */
slot_t * tdma_get_slot(tdma_t * tdma, u32 slot_idx){
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


/* Get a TX frame from the callback, add UW, and schedule for transmission */
void tdma_do_tx_frame(tdma_t * tdma, int slot_idx){
    struct TDMA_MODE_SETTINGS mode = tdma->settings;
    size_t frame_size = mode.frame_size;
    size_t slot_size = mode.slot_size;
    size_t bits_per_sym = (mode.fsk_m==2)?1:2;
    slot_t * slot = tdma_get_slot(tdma,slot_idx);
    size_t frame_size_bits = bits_per_sym*frame_size;
    size_t nbits = (slot_size+1)*bits_per_sym;
    size_t n_slots = mode.n_slots;
    u32 Rs = mode.sym_rate;
    u32 Fs = mode.samp_rate;
    u32 Ts = Fs/Rs;

    COMP mod_samps[(slot_size+1)*Ts];
    u8 frame_bits[frame_size_bits];
    u8 mod_bits[nbits];
    u8 uw_type = 0;
    if(slot == NULL) return;

    /* Clear bit buffer */
    memset(&mod_bits[0],0,nbits*sizeof(u8));

    /* Get a frame, or leave blank if call not setup */
    if(tdma->tx_callback != NULL){
        int ret = tdma->tx_callback(frame_bits,slot_idx,slot,tdma,&uw_type,tdma->tx_cb_data);
        if(!ret){
            slot->state = rx_no_sync;
            return;
        }
        if(uw_type > 1)
            uw_type = 0;
    }

    /* Copy frame bits to front of mod bit buffer */
    memcpy(&mod_bits[0],&frame_bits[0],frame_size_bits*sizeof(u8));

    const u8 * uw = &(TDMA_UW_LIST_A[uw_type][0]);
    /* Copy UW into frame */
    size_t uw_offset = (frame_size_bits-mode.uw_len)/2;
    memcpy(&mod_bits[uw_offset],uw,mode.uw_len*sizeof(u8));

    /* Modulate frame */
    fsk_mod_c(slot->fsk,mod_samps,mod_bits);

    /* Calculate TX time and send frame down to radio */
    /* timestamp of head of slot currently being demod'ed */
    i64 tx_timestamp = tdma->timestamp + (i64)tdma->sample_sync_offset;

    /* Figure out how far in future the next instance of 'this' slot will be */
    int delta_slots = 0;
    u32 slot_cur = tdma->slot_cur;
    if(slot_idx > slot_cur){
        delta_slots = slot_idx-slot_cur;
    }else if(slot_idx < slot_cur){
        delta_slots = (slot_cur+slot_idx+2)%n_slots;
    }
    
    /* Point timestamp to next available slot */
    tx_timestamp += delta_slots*slot_size*Ts;

    /* Add multi-slot/frame offset to tx timestamp */
    tx_timestamp += tdma->tx_multislot_delay*slot_size*Ts;

    /* Send frame on to radio if callback is setup */
    if(tdma->tx_burst_callback != NULL){
        tdma->tx_burst_callback(tdma,mod_samps,Ts*frame_size,tx_timestamp,tdma->tx_burst_cb_data);
    }
}

/* Pull TDMA frame out of bit stream and call RX CB if present */
void tdma_deframe_cbcall(u8 demod_bits[], u32 slot_i, tdma_t * tdma, slot_t * slot){
    size_t frame_size = tdma->settings.frame_size;
    size_t slot_size = tdma->settings.slot_size;
    size_t uw_len = tdma->settings.uw_len;
    size_t bits_per_sym = (tdma->settings.fsk_m==2)?1:2;
    size_t n_demod_bits = (slot_size+1)*bits_per_sym;
    size_t frame_size_bits = bits_per_sym*frame_size;
    size_t delta,off;
    i32 f_start;
    u32 master_max = tdma->settings.mastersat_max;

    u8 frame_bits[frame_size_bits];
    /* Re-find UW in demod'ed slice */
    /* Should probably just be left to tdma_rx_pilot_sync */
    //off = fvhff_search_uw(demod_bits,n_demod_bits,TDMA_UW_V,uw_len,&delta,bits_per_sym);
    off = tdma_search_uw(tdma, demod_bits, n_demod_bits, &delta, NULL);
    f_start = off - (frame_size_bits-uw_len)/2;

    /* If frame is not fully in demod bit buffer, there's not much we can do */
    if( (f_start < 0) || ((f_start+frame_size_bits) > n_demod_bits)){
        return;
    }

    /* Extract bits */
    memcpy(&frame_bits[0],&demod_bits[f_start],frame_size_bits*sizeof(u8));

    /* Check to see if this is a master timing frame. */
    bool master_mode = false;
    if(frame_bits[tdma->master_bit_pos])
        master_mode = true;
    
    /* Handle counter fiddling for master timing frames */
    if(master_mode){
        slot->master_count++;
        if(slot->master_count > master_max)
            slot->master_count = master_max;
    }else{
        slot->master_count--;
        if(slot->master_count < 0)
            slot->master_count = 0;
    }

    /* Right now we're not actually deframing the bits */
    /* TODO: actually extract UW type */
    if(tdma->rx_callback != NULL){
        tdma->rx_callback(frame_bits,slot_i,slot,tdma,0,tdma->rx_cb_data);
    }
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
    size_t uw_len = mode.uw_len;
    slot_t * slot = tdma_get_slot(tdma,tdma->slot_cur);
    fsk_t * fsk = slot->fsk;
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
        #ifdef VERY_DEBUG
        fprintf(stderr,"Skipping\n");
        #endif
        return;
    }

    /* Do TX if this is a TX slot */
    if(slot->state == tx_client){
        tdma_do_tx_frame(tdma,tdma->slot_cur);
    }
    /* If we're set up to ignore RX during a TX frame, and we're in a TX frame, ignore RX */
    if(!(tdma->ignore_rx_on_tx && slot->state == tx_client))
    {
        /* Zero out tail end of bit buffer so we can get last symbol out of demod */
        /* TODO: This is a hack. Look into better burst mode support in FSK */
        size_t i;
        for(i = slot_samps; i< (slot_size+1)*Ts; i++){
            frame_samps[i].real = 0;
            frame_samps[i].imag = 0;
        }

        /* Flag to indicate whether or not we should re-do the demodulation */
        bool repeat_demod = false;
        int rdemod_offset = 0;
        size_t delta,off,uw_type;
        i32 f_start;
        i32 frame_offset;
        bool f_valid = false;
        /* Demod section in do-while loop so we can repeat once if frame is just outside of bit buffer */
        do{
            /* Pull out the frame and demod */
            memcpy(&frame_samps[0],&sample_buffer[tdma->sample_sync_offset+rdemod_offset],slot_samps*sizeof(COMP));

            /* Demodulate the frame */
            fsk_demod(fsk,bit_buf,frame_samps);

            off = tdma_search_uw(tdma, bit_buf, nbits, &delta, &uw_type);
            f_start = off- (frame_bits-uw_len)/2;

            /* Check frame tolerance and sync state*/
            if(slot->state == rx_sync){
                f_valid = delta <= tdma->settings.frame_sync_tol;
            }else if(slot->state == rx_no_sync){
                f_valid = delta <= tdma->settings.first_sync_tol;
            }

            /* Calculate offset (in samps) from start of frame */
            /* Note: FSK outputs one symbol from the last batch, so we have to account for that */
            i32 target_frame_offset = ((slot_size-frame_size)/2)*Ts;
            frame_offset = ((f_start-bits_per_sym)*(Ts/bits_per_sym)) - target_frame_offset;

            /* Flag a large frame offset as a bad UW sync */
            if( abs(frame_offset) > (slot_samps/4) )
                f_valid = false;
            
            if(f_valid && !repeat_demod)
                slot->slot_local_frame_offset = frame_offset;


            /* Check to see if the bits are outside of the demod buffer. If so, adjust and re-demod*/
            /* Disabled for now; will re-enable when worked out better in head */
            /* No, this is enabled. I have it more or less worked out */
            if(f_valid && !repeat_demod){
                if((f_start < bits_per_sym) || ((f_start+(bits_per_sym*frame_size)) > (bits_per_sym*(slot_size+1)))){
                    repeat_demod = true;
                }
            }else repeat_demod = false;


            #ifdef VERY_DEBUG
            if(repeat_demod){
                fprintf(stderr,"f_start: %d, Re-demod-ing\n",f_start);
            }
            if(f_valid){
                fprintf(stderr,"Good UW, type %d\n",uw_type);
            }else{
                fprintf(stderr,"Bad UW\n");
            }
            #endif

            rdemod_offset = frame_offset;
            
        }while(repeat_demod);

        i32 single_slot_offset = slot->slot_local_frame_offset;
        /* Flag indicating whether or not we should call the callback */
        bool do_frame_found_call = false;   

        /* Do single slot state machine */
        /* TODO: think about/play around with other fsk_clear_estimators() positions */
        if( slot->state == rx_sync){
            do_frame_found_call = true;
            if(!f_valid){   /* on bad UW, increment bad uw count and possibly unsync */
                slot->bad_uw_count++;
                slot->master_count--;
                if(slot->master_count < 0)
                    slot->master_count = 0;

                #ifdef VERY_DEBUG
                fprintf(stderr,"----BAD UW COUNT %d TOL %d----\n",slot->bad_uw_count,tdma->settings.frame_sync_baduw_tol);
                #endif
                if(slot->bad_uw_count >= tdma->settings.frame_sync_baduw_tol){
                    slot->state = rx_no_sync;
                    slot->master_count = 0;
                    do_frame_found_call = false;

                    #ifdef VERY_DEBUG
                    fprintf(stderr,"----DESYNCING----\n");
                    #endif
                }
            }else{ /* Good UW found */
                slot->bad_uw_count = 0;
                do_frame_found_call = true;
            }

            #ifdef VERY_DEBUG
            fprintf(stderr,"Slot %d: sunk\n",tdma->slot_cur);
            #endif
        }else if(slot->state == rx_no_sync){
            #ifdef VERY_DEBUG
            fprintf(stderr,"Slot %d: no sync\n",tdma->slot_cur);
            #endif
            if(f_valid ){
                slot->state = rx_sync;
                do_frame_found_call = true;
            }else{
                fsk_clear_estimators(fsk);
            }
        }

        if(do_frame_found_call){
            tdma_deframe_cbcall(bit_buf,tdma->slot_cur,tdma,slot);
        }

        #ifdef VERY_DEBUG
        /* Unicode underline for pretty printing */
        char underline[] = {0xCC,0xB2,0x00};

        fprintf(stderr,"slot: %d fstart:%d offset: %d delta: %d f1:%.3f EbN0:%f\n",
            tdma->slot_cur,f_start,off,delta,fsk->f_est[0],fsk->EbNodB);
        for(i=0; i<nbits; i++){
            fprintf(stderr,"%d",bit_buf[i]);
            if((i>off && i<=off+uw_len) || i==f_start || i==(f_start+frame_bits-1)){
                fprintf(stderr,underline);
            }
        }
        fprintf(stderr,"\n");
        #endif

        /* If we are in master sync mode, don't adjust demod timing. we ARE demod timing */
        /* TODO: check if any slots are in TX mode and if any have master sync. If we are
           TXing and don't have master sync, lock sync offset. Otherwise, TX frames jitter
           wildly */

        if(tdma->state != master_sync){
            /* Update slot offset to compensate for frame centering */
            /* Also check to see if any slots are master. If so, take timing from them */
            i32 offset_total = 0;
            i32 offset_slots = 0;
            i32 offset_master = 0;  /* Offset of master slot */
            i32 master_max = 0;     /* Highest 'master count' */
            for( i=0; i<n_slots; i++){
                /* Only check offset from valid frames */
                slot_t * i_slot = tdma_get_slot(tdma,i);
                if(i_slot->state == rx_sync){
                    i32 local_offset = i_slot->slot_local_frame_offset;
                    /* Filter out extreme spurious timing offsets */
                    if(abs(local_offset)<(slot_samps/4)){
                        #ifdef VERY_DEBUG
                        fprintf(stderr,"Local offset: %d\n",local_offset);
                        #endif
                        offset_total+=local_offset;
                        offset_slots++;
                        if(i_slot->master_count > master_max){
                            master_max = i_slot->master_count;
                            offset_master = local_offset;
                        }
                    }
                }
            }
            /* Check for zero here; otherwise we get div-by-zero errors */
            offset_total = offset_slots>0 ? offset_total/offset_slots:0;
            /* Use master slot for timing if available, otherwise take average of all frames */
            if(master_max >= mode.mastersat_min){
                tdma->sample_sync_offset +=  (offset_master/4);
                #ifdef VERY_DEBUG
                fprintf(stderr,"Syncing to master offset %d\n",tdma->sample_sync_offset);
                #endif
            }else{
                tdma->sample_sync_offset +=  (offset_total/4);
                #ifdef VERY_DEBUG
                fprintf(stderr,"Total Offset:%d\n",offset_total);
                fprintf(stderr,"Slot offset: %d of %d\n",tdma->sample_sync_offset,slot_samps*n_slots);
                #endif
            }
            #ifdef VERY_DEBUG
            fprintf(stderr,"\n");
            #endif
        }
    }

    tdma->slot_cur++;
    if(tdma->slot_cur >= n_slots)
        tdma->slot_cur = 0;

    
    /* Compensate for frame timing offset sliding towards the start of buffer */
    /* Do so by running this again, demodulating two slots, and moving the slot offset one slot forward */
    if( slot_offset < (slot_samps/4) ){
        #ifdef VERY_DEBUG
        fprintf(stderr,"Recursing\n");
        #endif
        tdma->sample_sync_offset = tdma->sample_sync_offset + slot_samps;
        tdma_rx_pilot_sync(tdma);
    }
}

/* Attempt at 'plot modem' search for situations where no synchronization is had */
/* This currently preforms worse than just running tdma_rx_pilot_sync on every frame */
/* It may still be needed to acquire first frames -- more testing is needed */
void tdma_rx_no_sync(tdma_t * tdma, COMP * samps, u64 timestamp){
    fprintf(stderr,"searching for pilot\n");
    struct TDMA_MODE_SETTINGS mode = tdma->settings;
    COMP * sample_buffer = tdma->sample_buffer;
    u32 Rs = mode.sym_rate;
    u32 Fs = mode.samp_rate;
    u32 slot_size = mode.slot_size;
    u32 frame_size = mode.frame_size;
    u32 n_slots = mode.n_slots;
    u32 M = mode.fsk_m;
    u32 Ts = Fs/Rs;
    u32 bits_per_sym = M==2?1:2;
    u32 samps_per_slot = slot_size*Ts;
    size_t i, delta, offset, f_start;
    fsk_t * fsk = tdma->fsk_pilot;
    size_t uw_len = mode.uw_len;
    size_t frame_bits = frame_size*bits_per_sym;

    //Number of bits per pilot modem chunk (half a slot)
    u32 n_pilot_bits = (slot_size/2)*bits_per_sym;
    //u32 n_pilot_bits = (slot_size)*bits_per_sym;
    //We look at a full slot for the UW
    u8 pilot_bits[n_pilot_bits];

    /* Start search at the last quarter of the previously rx'ed slot's worth of samples */
    size_t search_offset_i = (3*samps_per_slot)/4;
    size_t best_match_offset;
    u32 best_delta = uw_len;
    /* Search every half slot at quarter slot offsets */
    for(i = 0; i < 4; i++){
        fsk_clear_estimators(fsk);
        fsk_demod(fsk,pilot_bits,&sample_buffer[search_offset_i]);
        fsk_demod(fsk,pilot_bits,&sample_buffer[search_offset_i]);
        
        offset = tdma_search_uw(tdma, pilot_bits, n_pilot_bits, &delta, NULL);
        f_start = offset - (frame_bits-uw_len)/2;

        fprintf(stderr,"delta: %d offset %d so:%d\n",delta,offset,search_offset_i);
        if(delta <= mode.pilot_sync_tol){
        }
        search_offset_i += samps_per_slot/4;
        if(delta<best_delta){
            best_delta = delta;
            best_match_offset = f_start + search_offset_i;
        }
    }
    if(best_delta <= mode.pilot_sync_tol){
        fprintf(stderr,"Pilot got UW delta %d search offset %d\n",best_delta,best_match_offset);
        tdma->sample_sync_offset = best_match_offset;
        tdma_rx_pilot_sync(tdma);
    }

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
    memmove(&sample_buffer[0],&sample_buffer[slot_samps],move_samps);

    move_samps = slot_samps*sizeof(COMP);
    memcpy(&sample_buffer[n_slots*slot_samps],&samps[0],move_samps);

    /* Set the timestamp. Not sure if this makes sense */
    tdma->timestamp = timestamp - (slot_samps*(n_slots-1));

    /* Staate machine for TDMA modem */
    switch(tdma->state){
        case no_sync:
            //tdma_rx_no_sync(tdma,samps,timestamp);
            //break;
        //case pilot_sync:
        case slot_sync:
        case master_sync:
            tdma_rx_pilot_sync(tdma);
            break;
        default:
            tdma->state = no_sync;
            break;
    }

    /* Check to see if we should change overall TDMA state */
    bool have_slot_sync = false;    /* Are any slots sunk? */
    slot_t * slot = tdma->slots;
    while(slot != NULL){
        have_slot_sync = have_slot_sync || (slot->state == rx_sync);
        slot = slot->next_slot;
    }
    /* Reset slot miss counter */
    if(have_slot_sync){
        tdma->sync_misses = 0;
    }

    /* Go from no sync to slot sync if we have slots sunk */
    if(have_slot_sync && tdma->state == no_sync){
        tdma->state = slot_sync;
    }

    /* Deal with case where no slots are sunk but we are coming off of sync */
    if( (!have_slot_sync) && tdma->state == slot_sync){
        tdma->sync_misses++;
        if( tdma->sync_misses > (mode.loss_of_sync_frames*n_slots)){
            tdma->state = no_sync;
        }
    }

    /* If we have no sync and no idea, nudge slot offset a bit so maybe we'll line up with any active frames */
    if( (!have_slot_sync) && (tdma->state == no_sync)){
        tdma->sample_sync_offset += (slot_samps/8);
    }
}


void tdma_set_rx_cb(tdma_t * tdma,tdma_cb_rx_frame rx_callback,void * cb_data){
    tdma->rx_callback = rx_callback;
    tdma->rx_cb_data = cb_data;
}


void tdma_set_tx_cb(tdma_t * tdma,tdma_cb_tx_frame tx_callback,void * cb_data){
    tdma->tx_callback = tx_callback;
    tdma->tx_cb_data = cb_data;
}


void tdma_set_tx_burst_cb(tdma_t * tdma,tdma_cb_tx_burst tx_burst_callback, void * cb_data){
    tdma->tx_burst_callback = tx_burst_callback;
    tdma->tx_burst_cb_data = cb_data;
}

/* Set up TDMA to schedule the transmission of a single frame. The frame itself will be 
    passed in through the tx_frame callback
*/
void tdma_single_frame_tx(tdma_t * tdma, int slot_idx){
    slot_t * slot = tdma_get_slot(tdma,slot_idx);
    if(slot == NULL) return;
    slot->state = tx_client;
    slot->single_tx = true;
}

/* Start transmission of a bunch of frames on a particular slot
*/
void tdma_start_tx(tdma_t * tdma, int slot_idx){
    slot_t * slot = tdma_get_slot(tdma,slot_idx);
    if(slot == NULL) return;
    slot->state = tx_client;
    slot->single_tx = false;
}


/* Stop ongoing transmission of a bunch of frames for some slot
*/
void tdma_stop_tx(tdma_t * tdma, int slot_idx){
    slot_t * slot = tdma_get_slot(tdma,slot_idx);
    if(slot == NULL) return;
    slot->state = rx_no_sync;
    slot->single_tx = false;
}

size_t tdma_nin(tdma_t * tdma){
    struct TDMA_MODE_SETTINGS mode = tdma->settings;
    u32 Rs = mode.sym_rate;
    u32 Fs = mode.samp_rate;
    u32 slot_size = mode.slot_size;
    u32 Ts = Fs/Rs;
    u32 slot_samps = slot_size*Ts;

    return slot_samps;
}

size_t tdma_nout(tdma_t * tdma){
    struct TDMA_MODE_SETTINGS mode = tdma->settings;
    size_t frame_size = mode.frame_size;
    u32 Rs = mode.sym_rate;
    u32 Fs = mode.samp_rate;
    u32 Ts = Fs/Rs;
    return frame_size*Ts;
}

#pragma GCC diagnostic pop
