
/*---------------------------------------------------------------------------*\

  FILE........: stm32f4_timer.c
  AUTHOR......: David Rowe
  DATE CREATED: May 2 2013

  STM32F4 implementation of the machine dependant timer functions.

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

#include <string.h>
#include "machdep.h"
#include "gdb_stdio.h"

volatile unsigned int *DWT_CYCCNT   = (volatile unsigned int *)0xE0001004; 
volatile unsigned int *DWT_CONTROL  = (volatile unsigned int *)0xE0001000;
volatile unsigned int *SCB_DEMCR    = (volatile unsigned int *)0xE000EDFC;
 
#define CORE_CLOCK 168E6
#define BUF_SZ     4096

static char buf[BUF_SZ];

void machdep_timer_init(void)
{
    static int enabled = 0;
 
    if (!enabled) {
        *SCB_DEMCR = *SCB_DEMCR | 0x01000000;
        *DWT_CYCCNT = 0; // reset the counter
        *DWT_CONTROL = *DWT_CONTROL | 1 ; // enable the counter
 
        enabled = 1;
    }
    *buf = 0;
}

void machdep_timer_reset(void)
{
    *DWT_CYCCNT = 0; // reset the counter
}

unsigned int machdep_timer_sample(void) {
    return *DWT_CYCCNT;
}

/* log to a buffer, we only call printf after timing finished as it is slow */

unsigned int machdep_timer_sample_and_log(unsigned int start, char s[])
{
    char tmp[80];

    unsigned int dwt = *DWT_CYCCNT - start;
    sprintf(tmp, "%s %5.2f msecs\n",s,1000.0*(float)dwt/CORE_CLOCK);
    if ((strlen(buf) + strlen(tmp)) < BUF_SZ)
        strcat(buf, tmp);
    return *DWT_CYCCNT;
}

void machdep_timer_print_logged_samples(void)
{
    gdb_stdio_printf("%s", buf);
    *buf = 0;
}
