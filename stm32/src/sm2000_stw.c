/*---------------------------------------------------------------------------*\

  FILE........: sm2000_stw.c
  AUTHOR......: David Rowe
  DATE CREATED: June 2016

  Test program to support SM2000 "set to work" - tetsing and bring up.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2016 David Rowe

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
#include <stdint.h>
#include "new_i2c.h"
#include "si53xx.h"
#include "debugblinky.h"
#include "sm1000_leds_switches.h"
#include "stm32f4_dac.h"
#include "stm32f4_adc.h"
//#include "stm32f4_usb_vcp.h"
#include "fsk.h"
#include "freedv_vhf_framing.h"
#include "golay23.h"
#include "string.h"
#define SINE_SAMPLES   24

/* 32 sample sine wave which at Fs=16kHz will be 500Hz.  Note samples
   are 16 bit 2's complement, the DAC driver convertsto 12 bit
   unsigned. */

short aSine[] = {
     5000,  -2500,  -2500, 5000,  -2500,  -2500, 
     5000,  -2500,  -2500, 5000,  -2500,  -2500, 
     5000,  -2500,  -2500, 5000,  -2500,  -2500, 
     5000,  -2500,  -2500, 5000,  -2500,  -2500
};
short zeros[SINE_SAMPLES];

uint8_t bit_buf[] = { 1,0,1,1,0,0,0,1,
                  1,0,1,1,0,1,0,1,
                  0,1,0,0,1,0,1,1,
                  1,0,0,0,0,0,0,1,
                  1,1,0,0,0,0,0,1,
                  1,0,1,0,0,0,0,1,
                  1,0,0,1,0,0,0,1,
                  1,0,0,0,1,0,0,1,
                  1,0,0,0,0,1,0,1,
                  1,0,0,0,0,0,1,1,
                  1,0,0,0,0,0,1,1,
                  1,0,0,0,0,0,1,1, };
                  
int main(void) {
    int ptt, i;
    uint64_t freq_in_Hz_times_100;
    struct FSK * fsk;
    struct freedv_vhf_deframer * deframer;
    float * mod_buf;
    
    sm1000_leds_switches_init();
    led_pwr(1);
    
    fsk = fsk_create_hbr(96000,1200,10,4,32000-1800,1200);
    fsk_set_est_limits(fsk,29000,35000);
    deframer = fvhff_create_deframer(FREEDV_VHF_FRAME_A,1);
    if(fsk==NULL){
		led_err(1);
		while(1);
	}
    mod_buf = malloc(sizeof(float)*fsk->Nmem);
    init_debug_blinky();
    txrx_12V(0);

    I2C_Setup();
    si5351_init(0, SI5351_CRYSTAL_LOAD_6PF, 0);
    freq_in_Hz_times_100 = 1070000000ULL - 3200000ULL;
    si5351_set_freq(freq_in_Hz_times_100, 0, SI5351_CLK0);

    dac_open(DAC_FS_96KHZ, 2000);
    adc_open(ADC_FS_96KHZ, 2000);

    //usb_vcp_init();
    int mbptr = 0;
	int golay_ctr = 0;
	uint8_t c2_buffer[8];
	int k;
	int spstate;
	ptt = 1;
	int nin = fsk_nin(fsk);
	spstate = 0;
    while(1) {
		if(switch_ptt() && (!spstate)){
			ptt = ptt ? 0 : 1;
			mbptr = 0;
		}
		spstate = switch_ptt();
		
        led_ptt(ptt);
        txrx_12V(ptt);
        
        if(ptt){
			int n = dac1_free();
			if (n) {
				int16_t  buf[n];
				for(i=0; i<n && mbptr<fsk->N; i++,mbptr++) {
					buf[i] = (short)(mod_buf[mbptr]*700);
				}
				dac1_write((short*)buf, i);
			}
			if(mbptr>=fsk->N){
				/* Encode frame index into golay codeword to protect from test BER*/
				k = golay23_encode((golay_ctr)&0x0FFF);
				c2_buffer[5] = (k    )&0xFF;
				c2_buffer[1] = (k>>8 )&0xFF;
				c2_buffer[0] = (k>>16)&0x7F;
				/* Frame the bits */
				fvhff_frame_bits(FREEDV_VHF_FRAME_A, bit_buf, c2_buffer,NULL,NULL);
				/* Mod the FSK */
				fsk_mod(fsk,mod_buf,bit_buf);
				mbptr = 0;
				golay_ctr++;
			}
		}else{
			
			int n = adc1_samps();
			if (n) {
				int16_t buf[n];
				adc1_read(buf,n);
				for(i=0; i<n && mbptr<nin; i++,mbptr++){
					mod_buf[mbptr] = ((float)buf[i]);
				}
				if(mbptr>=nin){
					led_rt(1);
					fsk_demod(fsk,bit_buf,mod_buf);
					led_rt(0);
					if(fvhff_deframe_bits(deframer,c2_buffer,NULL,NULL,bit_buf)){
						led_err(1);		
					}else{
						//led_err(0);
						led_err(0);
					}
					mbptr = 0;
					if(fsk_nin(fsk)!=nin){
						led_ptt(1);
					}else{
						led_ptt(0);
					}
					nin = fsk_nin(fsk);
					
				}
			}
		}
    }
}
