/*---------------------------------------------------------------------------*\

  FILE........: codec2_fifo.c
  AUTHOR......: David Rowe
  DATE CREATED: Oct 15 2012

  A FIFO design useful in gluing the FDMDV modem and codec together in
  integrated applications.  The unittest/tfifo indicates these
  routines are thread safe without the need for syncronisation
  object, e.g. a different thread can read and write to a fifo at the
  same time.

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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef __STDC_NO_ATOMICS__
#include <stdint.h>
#include <stdatomic.h>
#endif /* !__STDC_NO_ATOMICS__ */

#include "codec2_fifo.h"

struct FIFO {
    short *buf;
#ifdef __STDC_NO_ATOMICS__
    short *pin;
    short *pout;
#else
    atomic_intptr_t pin;
    atomic_intptr_t pout;
#endif /* __STDC_NO_ATOMICS__ */
    int    nshort;
};

// standard create function
struct FIFO *codec2_fifo_create(int nshort) {
    short *buf = (short*)malloc(sizeof(short)*nshort);
    assert(buf != NULL);
    return codec2_fifo_create_buf(nshort, buf);
}

// alternate create function where buffer is externally supplied
struct FIFO *codec2_fifo_create_buf(int nshort, short* buf) {
    struct FIFO *fifo;
    assert(buf != NULL);
    fifo = (struct FIFO *)malloc(sizeof(struct FIFO));
    assert(fifo != NULL);

    fifo->buf = buf;
#ifdef __STDC_NO_ATOMICS__
    fifo->pin = fifo->buf;
    fifo->pout = fifo->buf;
#else
    atomic_store(&fifo->pin, (intptr_t)fifo->buf);
    atomic_store(&fifo->pout, (intptr_t)fifo->buf);
#endif /* __STDC_NO_ATOMICS__ */
    fifo->nshort = nshort;

    return fifo;
}

void codec2_fifo_destroy(struct FIFO *fifo) {
    assert(fifo != NULL);
    free(fifo->buf);
    free(fifo);
}

static int codec2_fifo_used_internal(short* pin, short* pout, int nshort)
{
    int used = 0;
    
    if (pin >= pout)
        used = pin - pout;
    else
        used = nshort + (unsigned int)(pin - pout);

    return used;
}

static int codec2_fifo_free_internal(short* pin, short* pout, int nshort)
{
    // available storage is one less than nshort as prd == pwr
    // is reserved for empty rather than full

    return nshort - codec2_fifo_used_internal(pin, pout, nshort) - 1;
}

int codec2_fifo_write(struct FIFO *fifo, short data[], int n) {
    assert(fifo != NULL);
    assert(data != NULL);

#ifdef __STDC_NO_ATOMICS__
    short         *pin = fifo->pin;
    short         *pin = fifo->pout;
#else
    short         *pin = (short*)atomic_load(&fifo->pin);
    short         *pout = (short*)atomic_load(&fifo->pout);
#endif /* __STDC_NO_ATOMICS__ */
        
    int free = codec2_fifo_free_internal(pin, pout, fifo->nshort);
    if (n > free) {
        return -1;
    }
    else {
        short         *pdata = data;
        if ((pin + n) >= (fifo->buf + fifo->nshort))
        {
            int firstSamples = fifo->buf + fifo->nshort - pin;
            memcpy(pin, pdata, firstSamples * sizeof(short));
            n -= firstSamples;
            pin = fifo->buf;
            pdata += firstSamples;
        }
        memcpy(pin, pdata, n * sizeof(short));
        pin += n;
#ifdef __STDC_NO_ATOMICS__
        fifo->pin = pin;
#else
        atomic_store(&fifo->pin, (intptr_t)pin);
#endif /* __STDC_NO_ATOMICS__ */
    }

    return 0;
}

int codec2_fifo_read(struct FIFO *fifo, short data[], int n)
{
    assert(fifo != NULL);
    assert(data != NULL);
    
#ifdef __STDC_NO_ATOMICS__
    short         *pin = fifo->pin;
    short         *pin = fifo->pout;
#else
    short         *pin = (short*)atomic_load(&fifo->pin);
    short         *pout = (short*)atomic_load(&fifo->pout);
#endif /* __STDC_NO_ATOMICS__ */

    int used = codec2_fifo_used_internal(pin, pout, fifo->nshort);
    if (n > used) {
        return -1;
    }
    else {
        short         *pdata = data;
        if ((pout + n) >= (fifo->buf + fifo->nshort))
        {
            int firstSamples = fifo->buf + fifo->nshort - pout;
            memcpy(pdata, pout, firstSamples * sizeof(short));
            n -= firstSamples;
            pout = fifo->buf;
            pdata += firstSamples;
        }
        memcpy(pdata, pout, n * sizeof(short));
        pout += n;
#ifdef __STDC_NO_ATOMICS__
        fifo->pout = pout;
#else
	atomic_store(&fifo->pout, (intptr_t)pout);
#endif /* __STDC_NO_ATOMICS__ */
    }

    return 0;
}

int codec2_fifo_used(const struct FIFO * const fifo)
{
#ifdef __STDC_NO_ATOMICS__
    short         *pin = fifo->pin;
    short         *pout = fifo->pout;
#else
    short         *pin = (short*)atomic_load(&fifo->pin);
    short         *pout = (short*)atomic_load(&fifo->pout);
#endif /* __STDC_NO_ATOMICS__ */
    
    return codec2_fifo_used_internal(pin, pout, fifo->nshort);
}

int codec2_fifo_free(const struct FIFO * const fifo)
{
#ifdef __STDC_NO_ATOMICS__
    short         *pin = fifo->pin;
    short         *pout = fifo->pout;
#else
    short         *pin = (short*)atomic_load(&fifo->pin);
    short         *pout = (short*)atomic_load(&fifo->pout);
#endif /* __STDC_NO_ATOMICS__ */
    
    return codec2_fifo_free_internal(pin, pout, fifo->nshort);
}
