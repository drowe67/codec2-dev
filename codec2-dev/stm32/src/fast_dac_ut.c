/*---------------------------------------------------------------------------*\

  FILE........: dac_ut.c
  AUTHOR......: David Rowe
  DATE CREATED: May 31 2013

  Plays a 500 Hz sine wave sampled at 16 kHz out of PA5 on a Discovery board,
  or the speaker output of the SM1000.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2013 David Rowe

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
#include <stdlib.h>
#include "stm32f4_dacduc.h"
#include "iir_duc.h"
#include "stm32f4xx.h"
#include <stm32f4xx_tim.h>
#include <stm32f4xx_rcc.h>
#include "gdb_stdio.h"
#include "comp.h"
#include <string.h> 
//#include "gmsk_test_dat_m4.h"
#define SINE_SAMPLES  32


/* 32 sample sine wave which at Fs=16kHz will be 500Hz.  Note samples
   are 16 bit 2's complement, the DAC driver convertsto 12 bit
   unsigned. */

short aWave[] = {4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,
	4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,};

short aSine[] = {1600, 3200, 1601, 0, 1600, 3200, 1601, 0, 1600, 3200, 1601, 0, 1600, 3200, 1601, 0, 1600, 3200, 1601, 0, 1600, 3200, 1601, 0, 1600, 3200, 1600, 0, 1600, 3200, 1601, 0
};

//Sine at Fs/4
float f4sine[] = {1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,
1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,
1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,
1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,1,0,-1,0,};

//Intermediate 80k real before tx
int tx_imm[DUC_N];

//Complex input to chain
#define COMP_IN_SZ (DUC_48N)
COMP comp_in[COMP_IN_SZ];

unsigned short outbuf[DAC_DUC_BUF_SZ];

void setup_timer()
{
    RCC_APB1PeriphClockCmd(RCC_APB1Periph_TIM2, ENABLE);

    TIM_TimeBaseInitTypeDef timerInitStructure;
    timerInitStructure.TIM_Prescaler = 84;
    timerInitStructure.TIM_CounterMode = TIM_CounterMode_Up;
    timerInitStructure.TIM_Period = 0x8FFFFFFF;
    timerInitStructure.TIM_ClockDivision = 0;
    timerInitStructure.TIM_RepetitionCounter = 0;
    TIM_TimeBaseInit(TIM2, &timerInitStructure);
    TIM_Cmd(TIM2, ENABLE);
}

int main(void) {
    int tstart,tup,tend,cyc,i;

    memset((void*)outbuf,0,sizeof(short)*DAC_DUC_BUF_SZ);
    setup_timer();
    fast_dac_open(2*DAC_DUC_BUF_SZ,2*DAC_BUF_SZ);
    tstart=tend=tup=cyc=0;
    //Initalize complex input with signal at zero
    for(i=0;i<COMP_IN_SZ;i++){
        comp_in[i].real=1;
        comp_in[i].imag=0;
    }
    while (1) {
	cyc++;
        //if(cyc>GMSK_TEST_LEN)
        //    cyc=0;
	if(cyc%10000==0){
                printf("48c80r takes %d uSecs\n",tup-tstart);
		printf("iir upconvert takes %d uSecs\n",tend-tup);
	}
        tstart = TIM_GetCounter(TIM2);

        upconv_48c_80r(comp_in,tx_imm,1);

	tup = TIM_GetCounter(TIM2);

	iir_upconv_fixp(tx_imm,outbuf);

	tend = TIM_GetCounter(TIM2);

        //Sit and spin until we can get more samples into the dac
	while(dac1_write((short*)outbuf,DAC_DUC_BUF_SZ)<0);
    }

}
