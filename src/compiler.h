/*
  Copyright (C) 2020 Jon Beniston, M7RCE
  All rights reserved.

  Compiler specific definitions to support compilers without C99.

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

#ifndef COMPILER_H
#define COMPILER_H

#ifdef NO_C99
/* Use alloca instead of variable length arrays */
#ifdef _MSC_VER
#include <malloc.h>
#define alloca _alloca
#else
#include <stdlib.h>
#endif
#endif

#include <complex.h>

#ifdef _MSC_VER
typedef _Fcomplex complexf_t;
#define complexf_0                      _FCbuild(0.0f, 0.0f)
#define complexf(real,imag)             _FCbuild(real, imag)
#define complexf_init(real,imag)        {real, imag}
#define complexf_mul(a,b)               _FCmulcc(a,b)
#define complexf_mulr(a,b)              _FCmulcr(a,b)
#define complexf_add(a,b)               _FCbuild(crealf(a)+crealf(b), cimagf(a)+cimagf(b))
#define complexf_addr(a,b)              _FCbuild(crealf(a)+(b), cimagf(a))
#define complexf_sub(a,b)               _FCbuild(crealf(a)-crealf(b), cimagf(a)-cimagf(b))
#define complexf_divr(a,b)              _FCbuild(crealf(a)/(b), cimagf(a)/(b))
#else
typedef complex float complexf_t;
#define complexf_0                      0.0f
#define complexf(real,imag)             ((real)+(imag)*I)
#define complexf_init(real,imag)        ((real)+(imag)*I)
#define complexf_mul(a,b)               ((a)*(b))
#define complexf_mulr(a,b)              ((a)*(b))
#define complexf_add(a,b)               ((a)+(b))
#define complexf_addr(a,b)              ((a)+(b))
#define complexf_sub(a,b)               ((a)-(b))
#define complexf_divr(a,b)              ((a)/(b))
#endif

#if __GNUC__
#define ATTRIBUTE_UNUSED __attribute__((unused))
#else
#define ATTRIBUTE_UNUSED
#endif

#endif /* COMPILER_H */
