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
#include "new_i2c.h"
#include "si53xx.h"
#include "debugblinky.h"
#include "sm1000_leds_switches.h"
#include "stm32f4_dac.h"
#include "stm32f4_adc.h"
#include "stm32f4_usb_vcp.h"

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

int main(void) {
    int ptt, i;
    uint64_t freq_in_Hz_times_100;

    for(i=0; i<SINE_SAMPLES; i++)
        aSine[i] *= 1.5;

    sm1000_leds_switches_init();
    led_pwr(1);
    led_ptt(0);

    init_debug_blinky();
    txrx_12V(0);

    I2C_Setup();
    si5351_init(0, SI5351_CRYSTAL_LOAD_6PF, 0);
    freq_in_Hz_times_100 = 1070000000ULL - 3200000ULL;
    // int ret = si5351_set_freq(freq_in_Hz_times_100, 0, SI5351_CLK0);
    si5351_set_freq(freq_in_Hz_times_100, 0, SI5351_CLK0);

    //dac_open(DAC_FS_96KHZ, 4*DAC_BUF_SZ);
    adc_open(ADC_FS_96KHZ, 10*ADC_BUF_SZ);
    
    usb_vcp_init();

    while(1) {
        ptt = switch_ptt();
        led_ptt(ptt);
        txrx_12V(ptt);

        /*
        if (ptt)
            dac1_write((short*)aSine, SINE_SAMPLES);
        else
            dac1_write((short*)zeros, SINE_SAMPLES);
        */

        /* read another buffer from USB when DAC empties */
        /* note assumes USB host write two bytes at a time */
        /* and assumes enough bytes are available, need to test that
           host is throttled appropriately */

        int n = adc1_samps();
        if (n) {
            uint16_t  buf[n];
            uint8_t   b;
            adc1_read((short*)buf,n);
            for(i=0; i<n; i++) {
				b = (uint8_t)(buf[i]&0xFF);
				VCP_put_char(b);
				b = (uint8_t)((buf[i]>>8)&0xFF);
				VCP_put_char(b);
                //VCP_get_char(&b);
                //uint16_t s = b << 8;
                //VCP_get_char(&b);
                //s += b;
                //buf[i] = s;
            }
            //dac1_write((short*)buf, n);
        }
                
    }
}
