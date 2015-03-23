 /*---------------------------------------------------------------------------*\

  FILE........: iir_duc.h
  AUTHOR......: Brady O'Brien
  DATE CREATED: 6 Mar 2015

  Interapolator/Filter for IF upconversion

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2015 Brady O'Brien

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

#ifndef __IIR_DUC_H
#define __IIR_DUC_H
#include "comp.h"

void iir_upconv(float modin[],unsigned short dac_out[]);
void iir_upconv_fixp(int modin[], unsigned short dac_out[]);
void upconv_48c_80r(COMP comp_8[],int real_80[],int count);
void upconv_8c_80r(COMP comp_8[],float real_80[],int count);

#endif
