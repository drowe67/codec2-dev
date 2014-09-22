/*---------------------------------------------------------------------------*\

  FILE........: sm1000_leds_switches.h
  AUTHOR......: David Rowe
  DATE CREATED: 18 July 2014

  Functions for controlling LEDs and reading switches on SM1000.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2014 David Rowe

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

#ifndef __LEDS_SWITCHES__
#define __LEDS_SWITCHES__

void sm1000_leds_switches_init(void);

void led_pwr(int state);
void led_ptt(int state);
void led_rt(int state);
void led_err(int state);
void not_cptt(int state);

int switch_ptt(void);
int switch_select(void);
int switch_back(void);

void ColorfulRingOfDeath(int code);

#endif
