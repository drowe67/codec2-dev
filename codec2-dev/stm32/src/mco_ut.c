/*
     mco_ut.c
     
     Slightly modified version of : 
       http://stm32f4-discovery.net/2014/10/library-40-output-clocks-stm32f4/
 
     Outputs the HSI oscillator on MCO1, pin PA8, which is SCL3 on the SM2000, can
     be probed around R124.  Used to track down a USB boot loader bug, suspect the
     HSI on this particular STM32F407 is out so the boot loader is failing on USB
     discovery/enumeration and forcing a reset.
*/

/* Include core modules */
#include "stm32f4xx.h"
/* Include my libraries here */
#include "defines.h"
#include "tm_stm32f4_mco_output.h"
 
int main(void) {
    /* Initialize system */
    SystemInit();
 
    /* Initialize MCO1 output, pin PA8 */
    TM_MCOOUTPUT_InitMCO1();
    
    /* Initialize MCO2 output, pin PC9 */
    TM_MCOOUTPUT_InitMCO2();
    
    /* Set MCO1 output = HSI with prescaler 2 = 16MHz / 2 = 8MHz*/
    TM_MCOOUTPUT_SetOutput1(TM_MCOOUTPUT1_Source_HSI, TM_MCOOUTPUT_Prescaler_2);
    
    /* Set MCO2 output = SYSCLK / 4 */
    TM_MCOOUTPUT_SetOutput2(TM_MCOOUTPUT2_Source_SYSCLK, TM_MCOOUTPUT_Prescaler_4);
 
    while (1) {
        
    }
}

