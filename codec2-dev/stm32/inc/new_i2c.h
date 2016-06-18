/* 
 * File:   new_i2c.h
 * Author: leon (zs6lmg@gmail.com or leon@lrlabs.com)
 *
 * Created on March 17, 2016, 6:09 PM
 * 
 * GNU license apply.
 * 
 */

/*
  Copyright (C) 2016 Leon

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

#ifndef NEWI2C_H
#define	NEWI2C_H

#include <stdint.h>

//
// I2C1 device Enable
//
#define I2C_D1 0
// I2C1 only PORTB
// SCL - B6 or B8
// SDA - B7 or B9
#if I2C_D1
#define I2C_DEVICE I2C1
#define I2C_DX_P_SCK GPIOB
#define I2C_DX_P_SDA GPIOB
#define I2C_DX_SCK 6
#define I2C_DX_SDA 9
// could not get macro expansion in gcc to play with us
#define I2C_DX_CLK_SCK RCC_AHB1Periph_GPIOB
#define I2C_DX_CLK_SDA RCC_AHB1Periph_GPIOB
#endif

//
// I2C2 device Enable
//
#define I2C_D2 0
// I2C2 on PORTB, PORTF and PORTH
// SCL - B10 F1 H4
// SDA - B11 F0 H5
#if I2C_D2
#define I2C_DEVICE I2C2
#define I2C_DX_P_SCK GPIOB
#define I2C_DX_P_SDA GPIOB
#define I2C_DX_SCK 10
#define I2C_DX_SDA 11
// could not get macro expansion in gcc to play with us
#define I2C_DX_CLK_SCK RCC_AHB1Periph_GPIOB
#define I2C_DX_CLK_SDA RCC_AHB1Periph_GPIOB
#endif

//
// I2C3 device Enable
//
#define I2C_D3 1
// I2C2 on PORTA and PORTH
// SCL - A8 H7
// SDA - C9 H8
#if I2C_D3
#define I2C_DEVICE I2C3
#define I2C_DX_P_SCK GPIOA
#define I2C_DX_P_SDA GPIOC
#define I2C_DX_SCK 8
#define I2C_DX_SDA 9
// could not get macro expansion in gcc to play with us
#define I2C_DX_CLK_SCK RCC_AHB1Periph_GPIOA
#define I2C_DX_CLK_SDA RCC_AHB1Periph_GPIOC
#endif




#define I2C_SPEED                100000
#define I2C_STIMEOUT             ((uint32_t)0x1000)
#define I2C_LTIMEOUT             ((uint32_t)(300 * 0x1000))

// software I2C model for testing
// 
#define I2Cmodel 1

//void I2C_GPIO_Init(void);
void I2C_Setup(void);

uint32_t I2C_NewWriteRegister(uint8_t Addr,uint8_t Register,uint8_t Value);
uint32_t I2C_NewWriteRegisterN(uint8_t Addr,uint8_t Register,uint8_t *Value,uint8_t N);
uint32_t I2C_NewReadRegister(uint8_t Addr,uint8_t Register);
uint32_t I2C_NewReadRegisterN(uint8_t Addr,uint8_t Register,uint8_t *buffer, uint8_t N);

#endif	/* NEWI2C_H */

