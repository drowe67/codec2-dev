/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: fifo.h
  AUTHOR......: David Rowe
  DATE CREATED: Oct 15 2012
                                                                             
  A FIFO design useful in gluing the FDMDV modem and codec together in
  integrated applications.

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

#ifndef __FIFO__
#define __FIFO__

struct FIFO;

struct FIFO *fifo_create(int nshort);
void fifo_destroy(struct FIFO *fifo);
int fifo_write(struct FIFO *fifo, short data[], int n);
int fifo_read(struct FIFO *fifo, short data[], int n);
int fifo_n(struct FIFO *fifo);

#endif
