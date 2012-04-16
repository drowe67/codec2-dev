/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: fdmdv.c
  AUTHOR......: David Rowe                                                          
  DATE CREATED: April 14 2012
                                                                             
  Functions that implement a Frequency Divison Multiplexed Modem for
  Digital Voice (FDMDV) over HF channels.
                                                                             
\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2012 David Rowe

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

/*---------------------------------------------------------------------------*\
                                                                             
                               DEFINES

\*---------------------------------------------------------------------------*/

#define FS                    8000  /* sample rate in Hz                                                    */
#define T                   (1/Fs)  /* sample period in seconds                                             */
#define RS                      50  /* symbol rate in Hz                                                    */
#define NC                      14  /* number of carriers                                                   */
#define NB                       2  /* Bits/symbol for QPSK modulation                                      */
#define RB              (Nc*Rs*Nb)  /* bit rate                                                             */
#define M                    Fs/Rs  /* oversampling factor                                                  */
#define NSYM                     4  /* number of symbols to filter over                                     */
#define FSEP                    75  /* Separation between carriers (Hz)                                     */
#define FCENTRE               1200  /* Centre frequency, Nc/2 carriers below this, Nc/2 carriers above (Hz) */
#define NT                       5  /* number of symbols we estimate timing over                            */
#define P                        4  /* oversample factor used for initial rx symbol filtering               */
#define NFILTER            (NSYM*M) /* size of tx/rx filters at sampel rate M                               */
#define NFILTERTIMING (M+Nfilter+M) /* filter memory used for resampling after timing estimation            */

#define NTEST_BITS        (Nc*Nb*4) /* length of test bit sequence */

/*---------------------------------------------------------------------------*\
                                                                             
                               STRUCT for States

\*---------------------------------------------------------------------------*/

struct FDMDV {
    int current_test_bit;
};

/*---------------------------------------------------------------------------*\
                                                                             
                               INCLUDES

\*---------------------------------------------------------------------------*/

#include "fdmdv.h"
#include "rn.h"
#include "testbits.h"

/*---------------------------------------------------------------------------*\
                                                                             
                               FUNCTIONS

\*---------------------------------------------------------------------------*/

static COMP cneg(COMP a)
{
    COMP res;

    res.real = -a.real;
    res.imag = -a.imag;

    return res;
}

static COMP cmult(COMP a, COMP b)
{
    COMP res;

    res.real = a.real*b.real - a.imag*b.imag;
    res.imag = a.real*b.imag + a.imag*b.real;

    return res;
}

static COMP cadd(COMP a, COMP b)
{
    COMP res;

    res.real = a.real + b.real;
    res.imag = a.imag + b.imag;

    return res;
}

static COMP cdot(COMP a[], COMP b[], int n)
{
    COMP res;
    int  i;
    
    for(i=0; i<n; i++) 
	res = cadd(res, cmult(a,b));

    return res;
}

static void cbuf_shift_update(COMP buf[], COMP update[], int buflen, int updatelen)
{
    int  i,j;
    
    for(i=0; i<buflen-updatelen; i++) 
	buf[i] = buf[updatelen+i];
    for(j=0; j<updatelen; j++) 
	buf[i] = update[j];
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: fdmdv_create	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 16/4/2012 

  Create and initialise an instance of the modem.  Returns a pointer
  to the modem states or NULL on failure.  One set of states is
  sufficient for a full duuplex modem.

\*---------------------------------------------------------------------------*/

struct FDMDV *fdmdv_create(void)
{
    struct FDMDV *fdmdv;

    fdmdv = (struct FDMDV*)malloc(sizeof(struct FDMDV));
    if (fdmdv == NULL)
	return NULL;
    
    return fdmdv;
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: fdmdv_destroy	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 16/4/2012

  Destroy an instance of the modem.

\*---------------------------------------------------------------------------*/

void codec2_destroy(struct FDMDV *fdmdv)
{
    assert(fdmdv != NULL);
    free(fdmdv);
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: bits_to_dqpsk_symbols()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 16/4/2012

  Generate Nc+1 QPSK symbols from vector of (1,Nc*Nb) input tx_bits.
  The Nc+1 symbol is the +1 -1 +1 .... BPSK sync carrier.

\*---------------------------------------------------------------------------*/

void bits_to_dqpsk_symbols(COMP tx_symbols[], COMP prev_tx_symbols[], int tx_bits[], int *pilot_bit)
{
    int c, msb, lsb;
    COMP j = {0.0,1.0};
    COMP minusj = {0.0,-1.0};

    /* map tx_bits to to Nc DQPSK symbols */

    for(c=0; c<NC; c++) {
	msb = tx_bits[2*c]; 
	lsb = tx_bits[2*c+1];
	if ((msb == 0) && (lsb == 0))
	    tx_symbols[c] = prev_tx_symbols[c];
	if ((msb == 0) && (lsb == 1))
	    tx_symbols[c] = cmult(j, prev_tx_symbols[c]);
	if ((msb == 1) && (lsb == 0))
	    tx_symbols[c] = cneg(prev_tx_symbols[c]);
	if ((msb == 1) && (lsb == 1))
	    tx_symbols[c] = cmult(cneg(j),prev_tx_symbols[c]);
    }

    /* +1 -1 +1 -1 BPSK sync carrier, once filtered becomes (roughly)
       two spectral lines at +/- Rs/2 */
 
    if (*pilot_bit)
	tx_symbols[Nc] = cneg(prev_tx_symbols[Nc]);
    else
	tx_symbols[Nc] = prev_tx_symbols[Nc];

    if (*pilot_bit) 
	*pilot_bit = 0;
    else
	*pilot_bit = 1;
}

