/*---------------------------------------------------------------------------*\

  FILE........: modem_probe.h
  AUTHOR......: Brady O'Brien
  DATE CREATED: 9 January 2016

  Library to easily extract debug traces from modems during development

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

#ifndef __MODEMPROBE_H
#define __MODEMPROBE_H
#include <stdint.h>
#include <stdlib.h>
#include "comp.h"


#ifdef MODEMPROBE_ENABLE

void modem_probe_init_int(char *modname, char *runname);
void modem_probe_close_int();

void modem_probe_samp_i_int(char * tracename,int samp[],size_t cnt);
void modem_probe_samp_f_int(char * tracename,float samp[],size_t cnt);
void modem_probe_samp_c_int(char * tracename,COMP samp[],size_t cnt);


static inline void modem_probe_init(char *modname,char *runname){
	modem_probe_init_int(modname,runname);
}

static inline void modem_probe_close(){
	modem_probe_close_int();
}

static inline void modem_probe_samp_i(char *tracename,int samp[],size_t cnt){
	modem_probe_samp_i_int(tracename,samp,cnt);
}


static inline void modem_probe_samp_f(char *tracename,float samp[],size_t cnt){
	modem_probe_samp_f_int(tracename,samp,cnt);
}	

static inline void modem_probe_samp_c(char *tracename,COMP samp[],size_t cnt){
	modem_probe_samp_c_int(tracename,samp,cnt);
}

#else

static inline void modem_probe_init(char *modname,char *runname){
	return;
}

static inline void modem_probe_close(){
	return;
}

static inline void modem_probe_samp_i(char *name,int samp[],size_t sampcnt){
	return;
}

static inline void modem_probe_samp_f(char *name,float samp[],size_t cnt){
	return;
}

static inline void modem_probe_samp_c(char *name,COMP samp[],size_t cnt){
	return;
}

#endif

#endif
