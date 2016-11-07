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
#include <stdio.h>
#include "new_i2c.h"
#include "si53xx.h"
#include "debugblinky.h"
#include "sm1000_leds_switches.h"
#include "stm32f4_dac.h"
#include "stm32f4_adc.h"
//#include "stm32f4_usb_vcp.h"
#include "fsk.h"
#include "freedv_vhf_framing.h"
#include "freedv_api.h"
#include "golay23.h"
#include "string.h"

#include "malloc.h"

extern void initialise_monitor_handles(void);
volatile size_t tos;

size_t stack_size(){
    volatile int x = 0;
    return tos-((size_t)&x);
}

//float modbuf[4040];
//short speechbuf[320];

int main(void) {
    int ptt;
    uint64_t freq_in_Hz_times_100;
    struct freedv * fdv;
    struct FSK * fsk;
    volatile int z = 0;
    tos = (size_t)&z;
    //initialise_monitor_handles();
    //printf("Testing\n");
    sm1000_leds_switches_init();
    led_pwr(1);
    
    init_debug_blinky();
    txrx_12V(0);

    
    I2C_Setup();
    si5351_init(0, SI5351_CRYSTAL_LOAD_6PF, 0);
    freq_in_Hz_times_100 = 1070000000ULL - 3200000ULL;
    // int ret = si5351_set_freq(freq_in_Hz_times_100, 0, SI5351_CLK0);
    si5351_set_freq(freq_in_Hz_times_100, 0, SI5351_CLK0);
    
    /* Open up the FreeDV thing */
	fdv = freedv_open(FREEDV_MODE_2400A);
	/* Set 96k samp rate to allow for 32k fsk signal */
    freedv_set_alt_modem_samp_rate(fdv,96000);
    /* Get the FSK struct and set up freq estimator params */
    fsk = freedv_get_fsk(fdv);
    fsk_set_est_limits(fsk,28000,34000);
    fsk->f1_tx = 32000-1800;
    
    float * modbuf = malloc(sizeof(float)*(freedv_get_n_max_modem_samples(fdv)+10));
    short * speechbuf = malloc(sizeof(short)*freedv_get_n_speech_samples(fdv));
    //short buf[ADC_BUF_SZ*4];
    dac_open(DAC_FS_16KHZ*2, DAC_BUF_SZ*4);
    adc_open(ADC_FS_96KHZ, ADC_BUF_SZ*11);
    //usb_vcp_init();
    int mbptr = 0;
	int spstate;
	ptt = 0;
	int nin = freedv_nin(fdv);
	spstate = 0;
    for(int i=0;i<freedv_get_n_max_modem_samples(fdv);i++){
        modbuf[i]=0;
    }
    while(1) {
		if(switch_ptt() && (!spstate)){
			ptt = ptt ? 0 : 1;
			mbptr = 0;
		}
		spstate = switch_ptt();
		
        txrx_12V(ptt);
        
        if(ptt){
		}else{
			int n = adc1_samps();
            if (n) {
                if((n+mbptr)>nin){
                    n = nin-mbptr;
                }
                //In it's own scope so buf is deallocated after it's used 
                {
                    short buf[n];
                    adc1_read(buf,n);
                    for(int i=0;i<n&&mbptr<nin;i++,mbptr++){
                        modbuf[mbptr] = buf[i];
                        //printf("%d\n",stack_size(tos));
                    }
                }
               
				if(mbptr>=nin){
					mbptr = 0;
					led_rt(1);
					freedv_floatrx(fdv,speechbuf,modbuf);
                    led_rt(0);
    
					if(freedv_get_sync(fdv)){
						led_err(1);
					}else{
						led_err(0);
					}
					if(freedv_nin(fdv)!=nin){
						led_ptt(1);
					}else{
						led_ptt(0);
					}
					nin = freedv_nin(fdv);
					
				}
			}
            if(dac2_free()>freedv_get_n_speech_samples(fdv)){
                dac2_write(speechbuf,freedv_get_n_speech_samples(fdv));
            }
		}
    }
}
