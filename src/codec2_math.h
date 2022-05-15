#ifndef CODEC2_MATH_H
#define CODEC2_MATH_H

//==========================================================================
// Name:            codec2_math.h
//
// Purpose:         A wrapper around architecture specific math libraries 
//                  used on embedded devices to improve Codec2 performance.
// Created:         May 15, 2022
// Authors:         Mooneer Salem
//
// License:
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License version 2.1,
//  as published by the Free Software Foundation.  This program is
//  distributed in the hope that it will be useful, but WITHOUT ANY
//  WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program; if not, see <http://www.gnu.org/licenses/>.
//
//==========================================================================

#if defined(__ARM_ARCH)
#include "arm_math.h"

#define codec2_dot_product_f32(left, right, len, result) (arm_dot_prod_f32(left, right, len, result))
#define codec2_complex_dot_product_f32(left, right, len, resultReal, resultImag) (arm_cmplx_dot_prod_f32(left, right, len, resultReal, resultImag))

#else

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

// Your embedded project must define these two methods if it is not running on ARM.
extern void codec2_dot_product_f32(float* left, float* right, size_t len, float* result);
extern void codec2_complex_dot_product_f32(float* left, float* right, size_t len, float* resultReal, float* resultImag);

#ifdef __cplusplus
}
#endif // __cplusplus
    
#endif // defined(__ARM_ARCH)

#endif // CODEC2_MATH_H