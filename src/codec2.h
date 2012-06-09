/*---------------------------------------------------------------------------*\

  FILE........: codec2.h
  AUTHOR......: David Rowe
  DATE CREATED: 21 August 2010

  Codec 2 fully quantised encoder and decoder functions.  If you want use 
  Codec 2, these are the functions you need to call.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2010 David Rowe

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

#ifdef __cplusplus
  extern "C" {
#endif

#ifndef __CODEC2__
#define  __CODEC2__

#ifdef __WIN32__
#ifdef __BUILDING_DLL__
#define WIN32SUPPORT __declspec(dllexport) __stdcall
#else
#define WIN32SUPPORT __declspec(dllimport) __stdcall
#endif
#else
#define WIN32SUPPORT
#endif

#define CODEC2_MODE_2400 0
#define CODEC2_MODE_1400 1
#define CODEC2_MODE_1200 2

struct CODEC2;

struct CODEC2 * WIN32SUPPORT codec2_create(int mode);
void WIN32SUPPORT codec2_destroy(struct CODEC2 *codec2_state);
void WIN32SUPPORT codec2_encode(struct CODEC2 *codec2_state, unsigned char * bits, short speech_in[]);
void WIN32SUPPORT codec2_decode(struct CODEC2 *codec2_state, short speech_out[], const unsigned char *bits);
int  WIN32SUPPORT codec2_samples_per_frame(struct CODEC2 *codec2_state);
int  WIN32SUPPORT codec2_bits_per_frame(struct CODEC2 *codec2_state);

#endif

#ifdef __cplusplus
}
#endif

