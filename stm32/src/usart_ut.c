/*
  From:
    http://stm32projectconsulting.blogspot.com/2013/04/stm32f4-discovery-usart-example.html
*/

#include <stm32f4xx.h>
#include <stm32f4xx_usart.h>

void init_usart(void){

 GPIO_InitTypeDef GPIO_InitStructure;
 USART_InitTypeDef USART_InitStructure;

 /* enable peripheral clock for USART3 */
 RCC_APB1PeriphClockCmd(RCC_APB1Periph_USART3, ENABLE);

 /* GPIOB clock enable */
 RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_GPIOB, ENABLE);

 /* GPIOA Configuration:  USART3 TX on PB10 */
 GPIO_InitStructure.GPIO_Pin = GPIO_Pin_10;
 GPIO_InitStructure.GPIO_Mode = GPIO_Mode_AF;
 GPIO_InitStructure.GPIO_Speed = GPIO_Speed_50MHz;
 GPIO_InitStructure.GPIO_OType = GPIO_OType_PP;
 GPIO_InitStructure.GPIO_PuPd = GPIO_PuPd_UP ;
 GPIO_Init(GPIOB, &GPIO_InitStructure);

 /* Connect USART3 pins to AF2 */
 // TX = PB10
 GPIO_PinAFConfig(GPIOB, GPIO_PinSource10, GPIO_AF_USART3);

 USART_InitStructure.USART_BaudRate = 9600;
 USART_InitStructure.USART_WordLength = USART_WordLength_8b;
 USART_InitStructure.USART_StopBits = USART_StopBits_1;
 USART_InitStructure.USART_Parity = USART_Parity_No;
 USART_InitStructure.USART_HardwareFlowControl = USART_HardwareFlowControl_None;
 USART_InitStructure.USART_Mode = USART_Mode_Tx;
 USART_Init(USART3, &USART_InitStructure);

 USART_Cmd(USART3, ENABLE); // enable USART3

}

void Delay(__IO uint32_t nCount)
{
  while(nCount--)
  {
  }
}

int main(void){

 init_usart();

 while(1){
  USART_SendData(USART3, 'h'); // defined in stm32f4xx_usart.h
  Delay(0x3FFFFF);
 }
}
