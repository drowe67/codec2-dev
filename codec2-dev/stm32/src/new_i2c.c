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

#include "stm32f4xx.h"
#include "stm32f4xx_conf.h"
#include "new_i2c.h"
//#include "cds_i2s.h"

//
//
//
#define CONCAT_(x,y) x##y
#define CONCAT(x,y) CONCAT_(x,y)
//
//
//

//
//
//
//
#define I2C_SDA_AFPIN_ID CONCAT(GPIO_PinSource,I2C_DX_SDA)
#define I2C_SCK_AFPIN_ID CONCAT(GPIO_PinSource,I2C_DX_SCK)
#define I2C_SDA_PIN_ID CONCAT(GPIO_Pin_,I2C_DX_SDA)
#define I2C_SCK_PIN_ID CONCAT(GPIO_Pin_,I2C_DX_SCK)


/* definition to expand macro then apply to pragma message */
// Verify that my shortcuts are actually working
//#define VALUE_TO_STRING(x) #x
//#define VALUE(x) VALUE_TO_STRING(x)
//
//#define VAR_NAME_VALUE(var) #var "="  VALUE(var)
//#pragma message(VAR_NAME_VALUE(I2C_SCK_AFPIN_ID))
//#pragma message(VAR_NAME_VALUE(I2C_SDA_AFPIN_ID))

/*
 * @brief Setup I2C interface pins
 * PB6 is SCL
 * PB9 is SDA
 * 
 * @param       None
 * 
 * @return      None
 */
void I2C_Setup(void) {
    I2C_InitTypeDef I2C_Init_S;
    GPIO_InitTypeDef GPIO_Init_P;
    // setup the clock for the GPIO device
    RCC_AHB1PeriphClockCmd(I2C_DX_CLK_SCK|I2C_DX_CLK_SDA, ENABLE);
    // enable clock for the i2c device
#if I2C_D1
    RCC_APB1PeriphClockCmd(RCC_APB1Periph_I2C1, ENABLE);
    // Performing a Reset
    //RCC_APB2PeriphClockCmd(RCC_APB2Periph_SYSCFG, ENABLE);
    //RCC_APB1PeriphResetCmd(RCC_APB1Periph_I2C1, ENABLE);
    //RCC_APB1PeriphResetCmd(RCC_APB1Periph_I2C1, DISABLE);
#elif I2C_D3
    RCC_APB1PeriphClockCmd(RCC_APB1Periph_I2C3, ENABLE);
#else
#error "Problem on I2C"
#endif
    // setup the pins
    GPIO_Init_P.GPIO_Pin=I2C_SDA_PIN_ID;
    GPIO_Init_P.GPIO_Mode=GPIO_Mode_AF;
    GPIO_Init_P.GPIO_Speed=GPIO_Speed_50MHz;
    GPIO_Init_P.GPIO_OType=GPIO_OType_OD;
    GPIO_Init_P.GPIO_PuPd=GPIO_PuPd_NOPULL;
    GPIO_Init(I2C_DX_P_SDA, &GPIO_Init_P);
    //
    //
    GPIO_Init_P.GPIO_Pin=I2C_SCK_PIN_ID;
    GPIO_Init(I2C_DX_P_SCK, &GPIO_Init_P);
    // assign alternate functions
#if I2C_D1
    GPIO_PinAFConfig(I2C_DX_P_SDA, I2C_SDA_AFPIN_ID, GPIO_AF_I2C1);
    GPIO_PinAFConfig(I2C_DX_P_SCK, I2C_SCK_AFPIN_ID, GPIO_AF_I2C1);
#elif I2C_D3
    GPIO_PinAFConfig(I2C_DX_P_SDA, I2C_SDA_AFPIN_ID, GPIO_AF_I2C3);
    GPIO_PinAFConfig(I2C_DX_P_SCK, I2C_SCK_AFPIN_ID, GPIO_AF_I2C3);
#else
#error "Device not defined"
#endif
    //
    I2C_Cmd(I2C_DEVICE, DISABLE);
    I2C_DeInit(I2C_DEVICE);
    //
    //
    I2C_Init_S.I2C_Mode=I2C_Mode_I2C;
    I2C_Init_S.I2C_DutyCycle=I2C_DutyCycle_2;
    I2C_Init_S.I2C_Ack=I2C_Ack_Enable;
    I2C_Init_S.I2C_OwnAddress1=0x00;
    I2C_Init_S.I2C_AcknowledgedAddress=I2C_AcknowledgedAddress_7bit;
    I2C_Init_S.I2C_ClockSpeed=I2C_SPEED;
    I2C_Cmd(I2C_DEVICE, ENABLE);
    I2C_Init(I2C_DEVICE, &I2C_Init_S);
}



uint32_t I2C_timeout(uint32_t retval) {
    I2C_GenerateSTOP(I2C_DEVICE, ENABLE);
    I2C_SoftwareResetCmd(I2C_DEVICE, ENABLE);
    I2C_SoftwareResetCmd(I2C_DEVICE, DISABLE);
    I2C_DeInit(I2C_DEVICE);
    I2C_Setup();
    return(retval);
}



/**
 * @brief  Writes a Byte to a given register through the control interface (I2C)
 * @param  Addr: I2C address to write to.
 * @param  Register: Register inside device to write to.
 * @param  Value: the Byte value to be written into destination register.
 * @return 0 all ok.
 *  0x101 i2c busy.
 *  0x102 master mode not selected.
 *  0x103 master transmitter mode.
 *  0x104 send data failed.
 *  0x105 no reply from device.
 */
uint32_t I2C_NewWriteRegister(uint8_t Addr, uint8_t Register, uint8_t Value) {
    uint32_t result = 0;

    /* check if bus is busy */
    __IO uint32_t Timeout = I2C_LTIMEOUT;
    while (I2C_GetFlagStatus(I2C_DEVICE, I2C_FLAG_BUSY )) {
        if ((Timeout--) == 0) return 0x101;
    }
    I2C_GenerateSTART(I2C_DEVICE, ENABLE);
    Timeout = I2C_STIMEOUT;
    while (!I2C_CheckEvent(I2C_DEVICE,I2C_EVENT_MASTER_MODE_SELECT)) {
        if ((Timeout--) == 0) return 0x102;
    }
    I2C_Send7bitAddress(I2C_DEVICE,Addr,I2C_Direction_Transmitter);
    Timeout = I2C_STIMEOUT;
    while (!I2C_CheckEvent(I2C_DEVICE,I2C_EVENT_MASTER_TRANSMITTER_MODE_SELECTED)) {
        if ((Timeout--) == 0) return 0x103;
    }
    I2C_SendData(I2C_DEVICE, Register);
    Timeout = I2C_STIMEOUT;
    while (!I2C_CheckEvent(I2C_DEVICE,I2C_EVENT_MASTER_BYTE_TRANSMITTING)) {
        if ((Timeout--) == 0) return 0x104;
    }
    I2C_SendData(I2C_DEVICE, Value);
    // Wait for reply
    Timeout = I2C_LTIMEOUT;
    while (!I2C_GetFlagStatus(I2C_DEVICE, I2C_FLAG_BTF)) {
        if ((Timeout--) == 0) return 0x105;
    }
    I2C_GenerateSTOP(I2C_DEVICE, ENABLE);
    return result;
}

/**
 * @brief  Writes more than 1 bytes to a given register through the control interface (I2C)
 * @param  Addr: I2C address to write to.
 * @param  Register: Register inside device to write to.
 * @param  Value: the Byte(s) to be written into destination register.
 * @param  N: Number of byte(s) to write
 * @return 0 all ok.
 *  0x101 i2c busy.
 *  0x102 master mode not selected.
 *  0x103 master transmitter mode.
 *  0x104 send data failed.
 *  0x105 no reply from device.
 */
uint32_t I2C_NewWriteRegisterN(uint8_t Addr, uint8_t Register, uint8_t *Value, uint8_t N) {
    uint32_t result = 0;

    /* check if bus is busy */
    __IO uint8_t i;
    __IO uint32_t Timeout = I2C_LTIMEOUT;
    while (I2C_GetFlagStatus(I2C_DEVICE, I2C_FLAG_BUSY )) {
        if ((Timeout--) == 0) return 0x101;
    }
    I2C_GenerateSTART(I2C_DEVICE, ENABLE);
    Timeout = I2C_STIMEOUT;
    while (!I2C_CheckEvent(I2C_DEVICE,I2C_EVENT_MASTER_MODE_SELECT)) {
        if ((Timeout--) == 0) return 0x102;
    }
    I2C_Send7bitAddress(I2C_DEVICE,Addr,I2C_Direction_Transmitter);
    Timeout = I2C_STIMEOUT;
    while (!I2C_CheckEvent(I2C_DEVICE,I2C_EVENT_MASTER_TRANSMITTER_MODE_SELECTED)) {
        if ((Timeout--) == 0) return 0x103;
    }
    I2C_SendData(I2C_DEVICE, Register);
    Timeout = I2C_STIMEOUT;
    while (!I2C_CheckEvent(I2C_DEVICE,I2C_EVENT_MASTER_BYTE_TRANSMITTING)) {
        if ((Timeout--) == 0) return 0x104;
    }
    for (i = 0; i < N; i++) {
        I2C_SendData(I2C_DEVICE, *(Value++));
        // Wait for reply
        Timeout = I2C_LTIMEOUT;
        while (!I2C_GetFlagStatus(I2C_DEVICE, I2C_FLAG_BTF)) {
            if ((Timeout--) == 0) return 0x105;
        }
    }
    I2C_GenerateSTOP(I2C_DEVICE, ENABLE);
    return result;
}



/**
 * @brief  Reads a Byte to a given register via I2C
 * @param  Addr: I2C address.
 * @param  Register: Register inside device to read.
 * @return value read.
 *  0x101 i2c busy.
 *  0x102 master mode not selected.
 *  0x103 master transmitter mode.
 *  0x104 send data failed.
 *  0x105 no reply from device.
 *  0x106 direction receive.
 *  0x107 data register not empty.
 *  0x108 waiting data.
 */
uint32_t I2C_NewReadRegister(uint8_t Addr,uint8_t Register) {
    uint32_t result = 0;
    __IO uint32_t Timeout = I2C_LTIMEOUT;
    while (I2C_GetFlagStatus(I2C_DEVICE,I2C_FLAG_BUSY)) {
        if ((Timeout--) == 0) return I2C_timeout(0x101);
    }
    I2C_GenerateSTART(I2C_DEVICE, ENABLE);
#if (I2Cmodel == 1)
    Timeout = I2C_STIMEOUT;
    while (!I2C_GetFlagStatus(I2C_DEVICE,I2C_FLAG_SB)) {
        if ((Timeout--) == 0) return I2C_timeout(0x102);
    }
    I2C_AcknowledgeConfig(I2C_DEVICE, DISABLE);
    I2C_Send7bitAddress(I2C_DEVICE, Addr, I2C_Direction_Transmitter);
    Timeout = I2C_STIMEOUT;
    while (!I2C_GetFlagStatus(I2C_DEVICE, I2C_FLAG_ADDR)) {
        if ((Timeout--) == 0) return I2C_timeout(0x103);
    }
    (void) I2C_DEVICE->SR2;
    Timeout = I2C_STIMEOUT;
    while (!I2C_GetFlagStatus(I2C_DEVICE, I2C_FLAG_TXE)) {
        if ((Timeout--) == 0) return I2C_timeout(0x104);
    }
    I2C_SendData(I2C_DEVICE, Register);
    Timeout = I2C_STIMEOUT;
    while ( (!I2C_GetFlagStatus(I2C_DEVICE, I2C_FLAG_TXE)) ||
            (!I2C_GetFlagStatus(I2C_DEVICE, I2C_FLAG_BTF))
          )  {
        if ((Timeout--) == 0) return I2C_timeout(0x105);
    }
    I2C_GenerateSTART(I2C_DEVICE, ENABLE);
    Timeout = I2C_STIMEOUT;
    while (!I2C_GetFlagStatus(I2C_DEVICE, I2C_FLAG_SB)) {
        if ((Timeout--) == 0) return I2C_timeout(0x106);
    }
    I2C_Send7bitAddress(I2C_DEVICE, Addr, I2C_Direction_Receiver);
    Timeout = I2C_STIMEOUT;
    while (I2C_GetFlagStatus(I2C_DEVICE, I2C_FLAG_ADDR) == RESET) {
        if ((Timeout--) == 0) return I2C_timeout(0x107);
    }
    (void) I2C_DEVICE->SR2;
    while (!I2C_GetFlagStatus(I2C_DEVICE, I2C_FLAG_RXNE)) {
        if ((Timeout--) == 0) return I2C_timeout(0x108);
    }
    I2C_GenerateSTOP(I2C_DEVICE, ENABLE);
    Timeout = I2C_STIMEOUT;
    result = I2C_ReceiveData(I2C_DEVICE );
    Timeout = I2C_STIMEOUT;
    while (I2C_DEVICE ->CR1 & I2C_CR1_STOP ) {
        if ((Timeout--) == 0) return I2C_timeout(0x109);
    }
    I2C_AcknowledgeConfig(I2C_DEVICE, ENABLE);
    I2C_ClearFlag(I2C_DEVICE, I2C_FLAG_AF );
#else
    Timeout = I2C_STIMEOUT;
    while (!I2C_CheckEvent(I2C_DEVICE,I2C_EVENT_MASTER_MODE_SELECT)) {
        if ((Timeout--) == 0) return 0x102;
    }
    I2C_Send7bitAddress(I2C_DEVICE,Addr,I2C_Direction_Transmitter);
    Timeout = I2C_STIMEOUT;
    while (!I2C_CheckEvent(I2C_DEVICE,I2C_EVENT_MASTER_TRANSMITTER_MODE_SELECTED)) {
        if ((Timeout--) == 0) return 0x103;
    }
    I2C_SendData(I2C_DEVICE, Register);
    Timeout = I2C_STIMEOUT;
    while (I2C_GetFlagStatus(I2C_DEVICE,I2C_FLAG_BTF) == RESET) {
            if ((Timeout--) == 0) return 0x104;
    }
    I2C_GenerateSTART(I2C_DEVICE, ENABLE);
    Timeout = I2C_STIMEOUT;
    while (!I2C_CheckEvent(I2C_DEVICE,I2C_EVENT_MASTER_MODE_SELECT)) {
        if ((Timeout--) == 0) return 0x105;
    }
    I2C_Send7bitAddress(I2C_DEVICE,Addr,I2C_Direction_Receiver);
    Timeout = I2C_STIMEOUT;
    while (I2C_GetFlagStatus(I2C_DEVICE, I2C_FLAG_ADDR) == RESET) {
        if ((Timeout--) == 0) return 0x106;
    }
    I2C_AcknowledgeConfig(I2C_DEVICE,DISABLE);
    (void) I2C_DEVICE->SR2;
    I2C_GenerateSTOP(I2C_DEVICE, ENABLE);
    Timeout = I2C_STIMEOUT;
    while (I2C_GetFlagStatus(I2C_DEVICE,I2C_FLAG_RXNE) == RESET) {
        if ((Timeout--) == 0) return 0x107;
    }
    result = I2C_ReceiveData(I2C_DEVICE );
    Timeout = I2C_STIMEOUT;
    while (I2C_DEVICE ->CR1 & I2C_CR1_STOP ) {
        if ((Timeout--) == 0) return 0x108;
    }
    I2C_AcknowledgeConfig(I2C_DEVICE, ENABLE);
    I2C_ClearFlag(I2C_DEVICE, I2C_FLAG_AF );
#endif
    return result;
}

/**
 * @brief  Reads more than one byte from a given register via I2C
 * @param  Addr: I2C address.
 * @param  Register: Register inside device to read.
 * @return value read.
 *  0x101 i2c busy.
 *  0x102 master mode not selected.
 *  0x103 master transmitter mode.
 *  0x104 send data failed.
 *  0x105 no reply from device.
 *  0x106 direction receive.
 *  0x107 data register not empty.
 *  0x108 waiting data.
 */
uint32_t I2C_NewReadRegisterN(uint8_t Addr,uint8_t Register,uint8_t *buffer, uint8_t N) {
    uint32_t result = 0;
    uint8_t     cnt;
    __IO uint32_t Timeout = I2C_LTIMEOUT;
    if ( (N==0)||(N>64) ) return 0x109;
    while (I2C_GetFlagStatus(I2C_DEVICE,I2C_FLAG_BUSY)) {
        if ((Timeout--) == 0) return 0x101;
    }
    I2C_GenerateSTART(I2C_DEVICE, ENABLE);
    Timeout = I2C_STIMEOUT;
    while (!I2C_CheckEvent(I2C_DEVICE,I2C_EVENT_MASTER_MODE_SELECT)) {
        if ((Timeout--) == 0) return 0x102;
    }
    if(N==1) {
        // ACK disable
        I2C_AcknowledgeConfig(I2C_DEVICE, DISABLE);
    } else {
        // ACK enable
        I2C_AcknowledgeConfig(I2C_DEVICE, ENABLE);
    }
    I2C_Send7bitAddress(I2C_DEVICE,Addr,I2C_Direction_Transmitter);
    Timeout = I2C_STIMEOUT;
    while (!I2C_CheckEvent(I2C_DEVICE,I2C_EVENT_MASTER_TRANSMITTER_MODE_SELECTED)) {
        if ((Timeout--) == 0) return 0x103;
    }
    (void) I2C_DEVICE->SR2;
    Timeout = I2C_STIMEOUT;
    while (!I2C_GetFlagStatus(I2C_DEVICE, I2C_FLAG_TXE)) {
        if ((Timeout--) == 0) return I2C_timeout(0x104);
    }
    I2C_SendData(I2C_DEVICE, Register);
    Timeout = I2C_STIMEOUT;
    while ( (!I2C_GetFlagStatus(I2C_DEVICE, I2C_FLAG_TXE)) ||
            (!I2C_GetFlagStatus(I2C_DEVICE, I2C_FLAG_BTF))
          )  {
        if ((Timeout--) == 0) return I2C_timeout(0x105);
    }
    I2C_GenerateSTART(I2C_DEVICE, ENABLE);
    Timeout = I2C_STIMEOUT;
    while (!I2C_GetFlagStatus(I2C_DEVICE, I2C_FLAG_SB)) {
        if ((Timeout--) == 0) return I2C_timeout(0x106);
    }
    I2C_Send7bitAddress(I2C_DEVICE, Addr, I2C_Direction_Receiver);
    Timeout = I2C_STIMEOUT;
    while (I2C_GetFlagStatus(I2C_DEVICE, I2C_FLAG_ADDR) == RESET) {
        if ((Timeout--) == 0) return I2C_timeout(0x107);
    }
    (void) I2C_DEVICE->SR2;
    for (cnt=0; cnt<N; cnt++) {
        if((cnt+1) >= N) {
            // ACK disable
            I2C_AcknowledgeConfig(I2C_DEVICE, DISABLE);
            // Stop-Sequenz
            I2C_GenerateSTOP(I2C_DEVICE, ENABLE);
        }
        Timeout = I2C_STIMEOUT;
        while (!I2C_GetFlagStatus(I2C_DEVICE, I2C_FLAG_RXNE)) {
            if ((Timeout--) == 0) return I2C_timeout(0x108);
        }
        result=I2C_ReceiveData(I2C_DEVICE);
        Timeout = I2C_STIMEOUT;
        while (I2C_DEVICE ->CR1 & I2C_CR1_STOP ) {
            if ((Timeout--) == 0) return I2C_timeout(0x109);
        }
        *buffer=(uint8_t)(result & 0xff);
        buffer++;
        Timeout = I2C_STIMEOUT;
    }
    I2C_AcknowledgeConfig(I2C_DEVICE, ENABLE);
    //I2C_ClearFlag(I2C_DEVICE, I2C_FLAG_AF );
    if (result<0x100) result=0;
    return result;
}

