/*---------------------------------------------------------------------------*\

  FILE........: usb_vcp_ut.c
  AUTHOR......: xenovacivus
  DATE CREATED: 31 August 2014

  USB Virtual COM Port (VCP) unit test that I found here:

    https://github.com/xenovacivus/STM32DiscoveryVCP

  Remarkably, it compiled and ran first time, and even the LEDs blink
  as advertised, they just happen to match the LEDs on the SM1000!

  When I fired up Minicom, it echoed characters:

    $ sudo minicom -D /dev/ttyACM0

\*---------------------------------------------------------------------------*/

#include <stm32f4xx.h>
#include <stm32f4xx_gpio.h>
#include "stm32f4_usb_vcp.h"
#include "sm1000_leds_switches.h"

volatile uint32_t ticker, downTicker;

int main(void) {

    sm1000_leds_switches_init();
    usb_vcp_init();
    SysTick_Config(SystemCoreClock/1000);

    while (1) {

        /* Blink the orange LED at 1Hz */

        if (500 == ticker) {
            GPIOD->BSRRH = GPIO_Pin_13;
        }
        else if (1000 == ticker) {
            ticker = 0;
            GPIOD->BSRRL = GPIO_Pin_13;
        }


        /* If there's data on the virtual serial port:
         *  - Echo it back
         *  - Turn the green LED on for 10ms
         */
        uint8_t theByte;
        if (VCP_get_char(&theByte)) {
            VCP_put_char(theByte);


            GPIOD->BSRRL = GPIO_Pin_12;
            downTicker = 10;
        }
        if (0 == downTicker) {
            GPIOD->BSRRH = GPIO_Pin_12;
        }
    }

    return 0;
}

/*
 * Interrupt Handler
 */

void SysTick_Handler(void)
{
	ticker++;
	if (downTicker > 0)
	{
		downTicker--;
	}
}

