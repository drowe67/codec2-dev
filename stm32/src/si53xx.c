/*------------------------------------------------------------------------------
 * 
 * Ported to stm32F4xx non c++ by Leon Lessing leon@lrlabs.com or zs6lmg@gmail.com
 * 
 * 
 * Copyright (C) 2015-2016 Jason Milldrum <milldrum@gmail.com>
 *                         Dana H. Myers <k6jq@comcast.net>
 *
 * Many defines derived from clk-si5351.h in the Linux kernel.
 * Sebastian Hesselbarth <sebastian.hesselbarth@gmail.com>
 * Rabeeh Khoury <rabeeh@solid-run.com>
 *
 * do_div() macro derived from /include/asm-generic/div64.h in
 * the Linux kernel.
 * Copyright (C) 2003 Bernardo Innocenti <bernie@develer.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * 
 *------------------------------------------------------------------------------
*/

#include "stm32f4xx.h"
#include "stm32f4xx_conf.h"
//#include "system.h"

#include "new_i2c.h"
//#include "inc/xprintf.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "si53xx.h"


void si5351_write(uint8_t REGaddr, uint8_t data) {
    // ignoring errors
    // Waiting for the bite
    Si5351_Config.I2C_ErrorCode=I2C_NewWriteRegister(Si5351_Config.I2C_add, REGaddr, data);

    /*
    printf("  si5351_write: REGaddr: 0x%02x data...: 0x%02x", REGaddr, data);
    if (Si5351_Config.I2C_ErrorCode > 0xff)
        printf("  ErrorCode: 0x%02x\n", Si5351_Config.I2C_ErrorCode);
    else
        printf("\n");
    */
}

void si5351_write_bulk(uint8_t REGaddr, uint8_t bytes, uint8_t *data) {
    Si5351_Config.I2C_ErrorCode=I2C_NewWriteRegisterN(Si5351_Config.I2C_add, REGaddr, data, bytes);
}


uint8_t si5351_read(uint8_t REGaddr) {
    uint8_t reg_val;
    Si5351_Config.I2C_ErrorCode=I2C_NewReadRegister(Si5351_Config.I2C_add,REGaddr);
    //printf("  si5351_read.: REGaddr: 0x%02x", REGaddr);
    if (Si5351_Config.I2C_ErrorCode>0xff) {
        reg_val=0;
        //printf(" ErrorCode: 0x%02x\n", Si5351_Config.I2C_ErrorCode);
    } else {
        reg_val=(uint8_t)(Si5351_Config.I2C_ErrorCode & 0xff);
        Si5351_Config.I2C_ErrorCode=0;
        //printf(" reg_val: 0x%02x\n", reg_val);
    }
    return reg_val;
}

/*------------------------------------------------------------------------------
 * @brief si5351_init(uint8_t xtal_load_c, uint32_t ref_osc_freq)
 * Setup communications to the Si5351 and set the crystal
 * load capacitance.
 *
 * @param I2C_Address - enter I2C address here, use 0 to use .h file defined value
 * 
 * @param xtal_load_c - Crystal load capacitance. Use the SI5351_CRYSTAL_LOAD_*PF
 * defines in the header file
 * 
 * @param ref_osc_freq - Crystal/reference oscillator frequency in 1 Hz increments.
 * Defaults to 25000000 if a 0 is used here.
 * 
 *------------------------------------------------------------------------------
 */
void si5351_init(uint8_t I2C_Address, uint8_t xtal_load_c, uint32_t ref_osc_freq) {
    //printf("si5351_init\n");

    Si5351_Config.clk0_freq=0;
    Si5351_Config.lock_plla = SI5351_CLKNONE;
    Si5351_Config.lock_pllb = SI5351_CLKNONE;
    Si5351_Config.clk0_int_mode = 0;
    Si5351_Config.clk1_int_mode = 0;
    Si5351_Config.clk2_int_mode = 0;
    Si5351_Config.plla_freq = 0;
    Si5351_Config.pllb_freq = 0;
    Si5351_Config.clk0_freq = 0;
    Si5351_Config.clk1_freq = 0;
    Si5351_Config.clk2_freq = 0;
    Si5351_Config.xtal_freq = SI5351_XTAL_FREQ;
    Si5351_Config.I2C_add = SI5351_BUS_BASE_ADDR << 1;
    if (I2C_Address != 0) Si5351_Config.I2C_add=I2C_Address << 1;
    // Set crystal load capacitance
    uint8_t reg_val = 0x12; // 0b010010 reserved value bits
    reg_val |= xtal_load_c;
    si5351_write(SI5351_CRYSTAL_LOAD, reg_val);

    // DR: test of I2C
    reg_val = si5351_read(SI5351_CRYSTAL_LOAD);
    //printf("reg_val: 0x%02x\n", reg_val);

    // Change the ref osc freq if different from default
    // Divide down if greater than 30 MHz
    if (ref_osc_freq != 0) {
        uint8_t reg_val;
        reg_val = si5351_read(SI5351_PLL_INPUT_SOURCE);
        //
        // Clear the bits first
        reg_val &= ~(SI5351_CLKIN_DIV_MASK);
        if(ref_osc_freq <= 30000000) {
            Si5351_Config.xtal_freq = ref_osc_freq;
            reg_val |= SI5351_CLKIN_DIV_1;
        } else if(ref_osc_freq > 30000000 && ref_osc_freq <= 60000000) {
            Si5351_Config.xtal_freq = ref_osc_freq / 2;
            reg_val |= SI5351_CLKIN_DIV_2;
        } else if(ref_osc_freq > 60000000 && ref_osc_freq <= 100000000) {
            Si5351_Config.xtal_freq = ref_osc_freq / 4;
            reg_val |= SI5351_CLKIN_DIV_4;
        }
        si5351_write(SI5351_PLL_INPUT_SOURCE, reg_val);
    }
    // Initialize the CLK outputs according to flowchart in datasheet
    // First, turn them off
    si5351_write(16, 0x80);
    si5351_write(17, 0x80);
    si5351_write(18, 0x80);

    // Turn the clocks back on...
    si5351_write(16, 0x0c);
    si5351_write(17, 0x0c);
    si5351_write(18, 0x0c);

    // Then reset the PLLs
    si5351_pll_reset(SI5351_PLLA);
    si5351_pll_reset(SI5351_PLLB);
}




/*------------------------------------------------------------------------------
 * @brief si5351_set_freq(uint64_t freq, uint64_t pll_freq, enum si5351_clock output)
 * Sets the clock frequency of the specified CLK output
 *
 * @param freq - Output frequency in 0.01 Hz increments, so for 10MHz 
 *               use 1000000000ULL
 * 
 * @param pll_freq - Frequency of the PLL driving the Multisynth
 *   Use a 0 to have the function choose a PLL frequency
 * 
 * @param clk - Clock output
 *   (use the si5351_clock enum)
 *------------------------------------------------------------------------------
 */
uint8_t si5351_set_freq(uint64_t freq, uint64_t pll_freq, enum si5351_clock clk) {
    struct Si5351RegSet ms_reg;
    enum si5351_pll target_pll;
    uint8_t write_pll = 0;
    uint8_t r_div = SI5351_OUTPUT_CLK_DIV_1;
    uint8_t int_mode = 0;
    uint8_t div_by_4 = 0;

    //printf("si5351_set_freq:\n");

    // PLL bounds checking
    if(pll_freq != 0) {
        if ((pll_freq < SI5351_PLL_VCO_MIN * SI5351_FREQ_MULT)
                || (pll_freq > SI5351_PLL_VCO_MAX * SI5351_FREQ_MULT)) {
            return 1;
        }
    }

    //printf("freq: 0x%0x 0x%0x\n", (uint32_t)(freq >> 32), (uint32_t)(freq & 0xffffffff));

    // Lower bounds check
    if(freq < SI5351_CLKOUT_MIN_FREQ * SI5351_FREQ_MULT) {
        freq = SI5351_CLKOUT_MIN_FREQ * SI5351_FREQ_MULT;
    }

    // Upper bounds check
    if(freq > SI5351_MULTISYNTH_MAX_FREQ * SI5351_FREQ_MULT) {
        freq = SI5351_MULTISYNTH_MAX_FREQ * SI5351_FREQ_MULT;
    }

    //printf("freq: 0x%0x 0x%0x\n", (uint32_t)(freq >> 32), (uint32_t)(freq & 0xffffffff));

    // Select the proper R div value
    r_div = si5351_select_r_div(&freq);

    // Calculate the synth parameters
    // If pll_freq is 0 and freq < 150 MHz, let the algorithm pick a PLL frequency
    if ((pll_freq) && (freq < SI5351_MULTISYNTH_DIVBY4_FREQ * SI5351_FREQ_MULT)) {
        si5351_multisynth_calc(freq, pll_freq, &ms_reg);
        write_pll = 0;
        div_by_4 = 0;
        int_mode = 0;

        switch(clk)      {
            case SI5351_CLK0:
                    Si5351_Config.clk0_freq = freq;
                    break;
            case SI5351_CLK1:
                    Si5351_Config.clk1_freq = freq;
                    break;
            case SI5351_CLK2:
                    Si5351_Config.clk2_freq = freq;
                    break;
            default:
                    break;
        }
    } else {
        // The PLL must be calculated and set by firmware when 150 MHz <= freq <= 160 MHz
        if(freq >= SI5351_MULTISYNTH_DIVBY4_FREQ * SI5351_FREQ_MULT) {
            pll_freq = si5351_multisynth_calc(freq, 0, &ms_reg);
            write_pll = 1;
            div_by_4 = 1;
            int_mode = 1;
        }

        // Determine which PLL to use
        // CLK0 gets PLLA, CLK1 gets PLLB
        // CLK2 gets PLLB if necessary
        // Only good for Si5351A3 variant at the moment
        switch(clk) {
            case SI5351_CLK0:
                    //printf("case SI5351_CLK0\n");
                    pll_freq = si5351_multisynth_calc(freq, 0, &ms_reg);
                    target_pll = SI5351_PLLA;
                    write_pll = 1;
                    si5351_set_ms_source(SI5351_CLK0, SI5351_PLLA);

                    //printf("pll_freq: 0x%0x 0x%0x\n", (uint32_t)(pll_freq >> 32), (uint32_t)(pll_freq & 0xffffffff));
                    //printf("freq: 0x%0x 0x%0x\n", (uint32_t)(freq >> 32), (uint32_t)(freq & 0xffffffff));

                    Si5351_Config.plla_freq = pll_freq;
                    Si5351_Config.clk0_freq = freq;
                    break;
            case SI5351_CLK1:
                    // Check to see if PLLB is locked due to other output being < 1.024 MHz or >= 112.5 MHz
                    if(Si5351_Config.lock_pllb == SI5351_CLK2) {
                        // We can't have a 2nd output < 1.024 MHz or >= 112.5 MHz on the same PLL unless exact same freq, so exit
                        if((freq >= SI5351_MULTISYNTH_SHARE_MAX * SI5351_FREQ_MULT
                                || freq < SI5351_CLKOUT_MIN_FREQ * SI5351_FREQ_MULT * 128)
                                && freq != Si5351_Config.clk2_freq) {
                            Si5351_Config.clk1_freq = 0;
                            return 1;
                        } else {
                            // Else, set multisynth to same PLL freq as CLK2
                            pll_freq = Si5351_Config.pllb_freq;
                            si5351_multisynth_calc(freq, pll_freq, &ms_reg);
                            write_pll = 0;
                            si5351_set_ms_source(SI5351_CLK1, SI5351_PLLB);
                        }
                    } else {
                        Si5351_Config.pllb_freq = pll_freq;
                        pll_freq = si5351_multisynth_calc(freq, 0, &ms_reg);
                        write_pll = 1;
                        si5351_set_ms_source(SI5351_CLK1, SI5351_PLLB);
                    }

                    if(freq >= SI5351_MULTISYNTH_SHARE_MAX * SI5351_FREQ_MULT
                            || freq < SI5351_CLKOUT_MIN_FREQ * SI5351_FREQ_MULT * 128) {
                        Si5351_Config.lock_pllb = SI5351_CLK1;

                        // Recalc and rewrite the multisynth parameters on CLK2
                        if(Si5351_Config.clk2_freq != 0) {
                            struct Si5351RegSet ms_temp_reg;
                            r_div = si5351_select_r_div(&Si5351_Config.clk2_freq);
                            si5351_multisynth_calc(Si5351_Config.clk2_freq, \
                                    Si5351_Config.pllb_freq, &ms_temp_reg);
                            si5351_set_ms(SI5351_CLK2, ms_temp_reg, 0, r_div, 0);
                        }
                    } else {
                        Si5351_Config.lock_pllb = SI5351_CLKNONE;
                    }

                    target_pll = SI5351_PLLB;
                    Si5351_Config.clk1_freq = freq;
                    break;
            case SI5351_CLK2:
                    // Check to see if PLLB is locked due to other output being < 1.024 MHz or >= 112.5 MHz
                    if(Si5351_Config.lock_pllb == SI5351_CLK1) {
                        // We can't have a 2nd output < 1.024 MHz  or >= 112.5 MHz on the same PLL unless exact same freq, so exit
                        if((freq >= SI5351_MULTISYNTH_SHARE_MAX * SI5351_FREQ_MULT
                                || freq < SI5351_CLKOUT_MIN_FREQ * SI5351_FREQ_MULT * 128)
                                && freq != Si5351_Config.clk2_freq) {
                            Si5351_Config.clk2_freq = 0;
                            return 1;
                        } else {
                            // Else, set multisynth to same PLL freq as CLK1
                            pll_freq = Si5351_Config.pllb_freq;
                            si5351_multisynth_calc(freq, pll_freq, &ms_reg);
                            write_pll = 0;
                            si5351_set_ms_source(SI5351_CLK2, SI5351_PLLB);
                        }
                    } else {
                        // need to account for CLK2 set before CLK1
                        Si5351_Config.pllb_freq = pll_freq;
                        pll_freq = si5351_multisynth_calc(freq, 0, &ms_reg);
                        write_pll = 1;
                        si5351_set_ms_source(SI5351_CLK2, SI5351_PLLB);
                    }

                    if(freq >= SI5351_MULTISYNTH_SHARE_MAX * SI5351_FREQ_MULT 
                            || freq < SI5351_CLKOUT_MIN_FREQ * SI5351_FREQ_MULT * 128) {
                        Si5351_Config.lock_pllb = SI5351_CLK2;

                        if(Si5351_Config.clk1_freq != 0) {
                            // Recalc and rewrite the multisynth parameters on CLK1
                            struct Si5351RegSet ms_temp_reg;
                            r_div = si5351_select_r_div(&Si5351_Config.clk1_freq);
                            si5351_multisynth_calc(Si5351_Config.clk1_freq, Si5351_Config.pllb_freq, &ms_temp_reg);
                            si5351_set_ms(SI5351_CLK1, ms_temp_reg, 0, r_div, 0);
                        }
                    } else {
                        Si5351_Config.lock_pllb = SI5351_CLKNONE;
                    }

                    target_pll = SI5351_PLLB;
                    Si5351_Config.clk2_freq = freq;
                    break;
            default:
                    return 1;
        }
    }

    // Set multisynth registers (MS must be set before PLL)
    si5351_set_ms(clk, ms_reg, int_mode, r_div, div_by_4);

    // Set PLL if necessary
    if (write_pll == 1) {
        si5351_set_pll(pll_freq, target_pll);
    }

    return 0;
}

/*------------------------------------------------------------------------------
 * @brief si5351_set_pll(uint64_t pll_freq, enum si5351_pll target_pll)
 *
 * Set the specified PLL to a specific oscillation frequency
 *
 * @param pll_freq - Desired PLL frequency
 * 
 * @param target_pll - Which PLL to set
 *     (use the si5351_pll enum)
 *------------------------------------------------------------------------------
 */
void si5351_set_pll(uint64_t pll_freq, enum si5351_pll target_pll) {
  struct Si5351RegSet pll_reg;

  si5351_pll_calc(pll_freq, &pll_reg, Si5351_Config.ref_correction);

  // Derive the register values to write

  // Prepare an array for parameters to be written to
  uint8_t *params;
  uint8_t buffer[20];
  uint8_t i = 0;
  uint8_t temp;
  params = (uint8_t *)&buffer;

  // Registers 26-27
  temp = ((pll_reg.p3 >> 8) & 0xFF);
  params[i++] = temp;

  temp = (uint8_t)(pll_reg.p3  & 0xFF);
  params[i++] = temp;

  // Register 28
  temp = (uint8_t)((pll_reg.p1 >> 16) & 0x03);
  params[i++] = temp;

  // Registers 29-30
  temp = (uint8_t)((pll_reg.p1 >> 8) & 0xFF);
  params[i++] = temp;

  temp = (uint8_t)(pll_reg.p1  & 0xFF);
  params[i++] = temp;

  // Register 31
  temp = (uint8_t)((pll_reg.p3 >> 12) & 0xF0);
  temp += (uint8_t)((pll_reg.p2 >> 16) & 0x0F);
  params[i++] = temp;

  // Registers 32-33
  temp = (uint8_t)((pll_reg.p2 >> 8) & 0xFF);
  params[i++] = temp;

  temp = (uint8_t)(pll_reg.p2  & 0xFF);
  params[i++] = temp;

  // Write the parameters
  if (target_pll == SI5351_PLLA) {
    si5351_write_bulk(SI5351_PLLA_PARAMETERS, i, params);
  } else if (target_pll == SI5351_PLLB) {
    si5351_write_bulk(SI5351_PLLB_PARAMETERS, i, params);
  }
}




/*------------------------------------------------------------------------------
 * @brief  si5351_pll_reset(enum si5351_pll target_pll)
 * target_pll - Which PLL to reset
 *     (use the si5351_pll enum)
 *
 * Apply a reset to the indicated PLL.
 * @param SI5351_PLL_RESET_A or SI5351_PLL_RESET_B
 * 
 * @return none
 *------------------------------------------------------------------------------
 */
void si5351_pll_reset(enum si5351_pll target_pll) {
    if(target_pll == SI5351_PLLA) {
        si5351_write(SI5351_PLL_RESET, SI5351_PLL_RESET_A);
    } else if(target_pll == SI5351_PLLB) {
        si5351_write(SI5351_PLL_RESET, SI5351_PLL_RESET_B);
    }
}

/*------------------------------------------------------------------------------
 * @brief set_int(enum si5351_clock clk, uint8_t int_mode)
 * Set the indicated multisynth into integer mode.
 *
 * @param clk - Clock output
 *   (use the si5351_clock enum)
 * 
 * @param enable - Set to 1 to enable, 0 to disable
 *
 *------------------------------------------------------------------------------
 */
void si5351_set_int(enum si5351_clock clk, uint8_t enable){
    uint8_t reg_val;
    reg_val = si5351_read(SI5351_CLK0_CTRL + (uint8_t)clk);

    if(enable == 1) {
            reg_val |= (SI5351_CLK_INTEGER_MODE);
    } else {
            reg_val &= ~(SI5351_CLK_INTEGER_MODE);
    }

    si5351_write(SI5351_CLK0_CTRL + (uint8_t)clk, reg_val);

    // Integer mode indication
    switch(clk) {
	case SI5351_CLK0:
		Si5351_Config.clk0_int_mode = enable;
		break;
	case SI5351_CLK1:
		Si5351_Config.clk1_int_mode = enable;
		break;
	case SI5351_CLK2:
		Si5351_Config.clk2_int_mode = enable;
		break;
	default:
		break;
    }
}


/*------------------------------------------------------------------------------
 * @brief si5351_set_clock_pwr(enum si5351_clock clk, uint8_t pwr)
 * Enable or disable power to a clock output (a power
 * saving feature).
 *
 * @param clk - Clock output
 *   (use the si5351_clock enum)
 * 
 * @param pwr - Set to 1 to enable, 0 to disable
 *
 *------------------------------------------------------------------------------
 */
void si5351_set_clock_pwr(enum si5351_clock clk, uint8_t pwr) {
    uint8_t reg_val;
    reg_val = si5351_read(SI5351_CLK0_CTRL + (uint8_t)clk);
    if(pwr == 1) {
        reg_val &= 0b01111111;
    } else {
        reg_val |= 0b10000000;
    }
    si5351_write(SI5351_CLK0_CTRL + (uint8_t)clk, reg_val);
}

/*------------------------------------------------------------------------------
 * @brief si5351_set_clock_invert(enum si5351_clock clk, uint8_t inv)
 * Enable to invert the clock output waveform.
 *
 * @param clk - Clock output
 *   (use the si5351_clock enum)
 * 
 * @param inv - Set to 1 to enable, 0 to disable
 *
 * -----------------------------------------------------------------------------
 */
void si5351_set_clock_invert(enum si5351_clock clk, uint8_t inv) {
    uint8_t reg_val;
    reg_val = si5351_read(SI5351_CLK0_CTRL + (uint8_t)clk);
    if(inv == 1) {
        reg_val |= (SI5351_CLK_INVERT);
    } else {
        reg_val &= ~(SI5351_CLK_INVERT);
    }
    si5351_write(SI5351_CLK0_CTRL + (uint8_t)clk, reg_val);
}


/*------------------------------------------------------------------------------
 * @brief si5351_set_ms(enum si5351_clock clk, struct Si5351RegSet ms_reg,
 *  uint8_t int_mode, uint8_t r_div, uint8_t div_by_4)
 *
 * Set the specified multisynth parameters. Not normally needed, but public for advanced users.
 *
 * @param clk - Clock output
 *   (use the si5351_clock enum)
 * 
 * @param int_mode - Set integer mode
 *  Set to 1 to enable, 0 to disable
 * 
 * @param r_div - Desired r_div ratio
 * 
 * @param div_by_4 - Set Divide By 4 mode
 *   Set to 1 to enable, 0 to disable
 *------------------------------------------------------------------------------
 */
void si5351_set_ms(enum si5351_clock clk, struct Si5351RegSet ms_reg, uint8_t int_mode, uint8_t r_div, uint8_t div_by_4) {
    uint8_t *params;
    uint8_t buffer[20];
    uint8_t i = 0;
    uint8_t temp;
    uint8_t reg_val;

    params = (uint8_t *)&buffer;
    // Registers 42-43 for CLK0
    temp = (uint8_t)((ms_reg.p3 >> 8) & 0xFF);
    params[i++] = temp;

    temp = (uint8_t)(ms_reg.p3  & 0xFF);
    params[i++] = temp;

    // Register 44 for CLK0
    reg_val = si5351_read((SI5351_CLK0_PARAMETERS + 2) + (clk * 8));
    reg_val &= ~(0x03);
    temp = reg_val | ((uint8_t)((ms_reg.p1 >> 16) & 0x03));
    params[i++] = temp;

    // Registers 45-46 for CLK0
    temp = (uint8_t)((ms_reg.p1 >> 8) & 0xFF);
    params[i++] = temp;

    temp = (uint8_t)(ms_reg.p1  & 0xFF);
    params[i++] = temp;

    // Register 47 for CLK0
    temp = (uint8_t)((ms_reg.p3 >> 12) & 0xF0);
    temp += (uint8_t)((ms_reg.p2 >> 16) & 0x0F);
    params[i++] = temp;

    // Registers 48-49 for CLK0
    temp = (uint8_t)((ms_reg.p2 >> 8) & 0xFF);
    params[i++] = temp;

    temp = (uint8_t)(ms_reg.p2  & 0xFF);
    params[i++] = temp;

    // Write the parameters
    switch (clk) {
        case SI5351_CLK0:
                si5351_write_bulk(SI5351_CLK0_PARAMETERS, i, params);
                break;
        case SI5351_CLK1:
                si5351_write_bulk(SI5351_CLK1_PARAMETERS, i, params);
                break;
        case SI5351_CLK2:
                si5351_write_bulk(SI5351_CLK2_PARAMETERS, i, params);
                break;
        case SI5351_CLK3:
                si5351_write_bulk(SI5351_CLK3_PARAMETERS, i, params);
                break;
        case SI5351_CLK4:
                si5351_write_bulk(SI5351_CLK4_PARAMETERS, i, params);
                break;
        case SI5351_CLK5:
                si5351_write_bulk(SI5351_CLK5_PARAMETERS, i, params);
                break;
        case SI5351_CLK6:
                si5351_write_bulk(SI5351_CLK6_PARAMETERS, i, params);
                break;
        case SI5351_CLK7:
                si5351_write_bulk(SI5351_CLK7_PARAMETERS, i, params);
                break;
        case SI5351_CLKNONE:
                return;
    }

    si5351_set_int(clk, int_mode);
    si5351_ms_div(clk, r_div, div_by_4);
}

void si5351_ms_div(enum si5351_clock clk, uint8_t r_div, uint8_t div_by_4) {
    uint8_t reg_val, reg_addr;

    switch(clk) {
        case SI5351_CLK0:
                reg_addr = SI5351_CLK0_PARAMETERS + 2;
                break;
        case SI5351_CLK1:
                reg_addr = SI5351_CLK1_PARAMETERS + 2;
                break;
        case SI5351_CLK2:
                reg_addr = SI5351_CLK2_PARAMETERS + 2;
                break;
        case SI5351_CLK3:
                reg_addr = SI5351_CLK3_PARAMETERS + 2;
                break;
        case SI5351_CLK4:
                reg_addr = SI5351_CLK4_PARAMETERS + 2;
                break;
        case SI5351_CLK5:
                reg_addr = SI5351_CLK5_PARAMETERS + 2;
                break;
        case SI5351_CLK6:
                return;
        case SI5351_CLK7:
                return;
        case SI5351_CLKNONE:
                return;
	default:
		return;
    }

    reg_val = si5351_read(reg_addr);

    // Clear the relevant bits
    reg_val &= ~(0x7c);

    if(div_by_4 == 0) {
            reg_val &= ~(SI5351_OUTPUT_CLK_DIVBY4);
    } else {
            reg_val |= (SI5351_OUTPUT_CLK_DIVBY4);
    }

    reg_val |= (r_div << SI5351_OUTPUT_CLK_DIV_SHIFT);

    si5351_write(reg_addr, reg_val);
}

/*------------------------------------------------------------------------------
 * @brief set_ms_source(enum si5351_clock clk, enum si5351_pll pll)
 * Set the desired PLL source for a multisynth.
 *
 * @param clk - Clock output
 *   (use the si5351_clock enum)
 * 
 * @param pll - Which PLL to use as the source
 *     (use the si5351_pll enum)
 *
 *------------------------------------------------------------------------------
 */
void si5351_set_ms_source(enum si5351_clock clk, enum si5351_pll pll) {
    uint8_t reg_val;

    reg_val = si5351_read(SI5351_CLK0_CTRL + (uint8_t)clk);

    if(pll == SI5351_PLLA) {
            reg_val &= ~(SI5351_CLK_PLL_SELECT);
    } else if (pll == SI5351_PLLB) {
            reg_val |= SI5351_CLK_PLL_SELECT;
    }
    si5351_write(SI5351_CLK0_CTRL + (uint8_t)clk, reg_val);
}

uint64_t si5351_multisynth_calc(uint64_t freq, uint64_t pll_freq, struct Si5351RegSet *reg) {
    uint64_t lltmp;
    uint32_t a, b, c, p1, p2, p3;
    uint8_t divby4;
    uint8_t ret_val = 0;

    // Multisynth bounds checking
    if (freq > SI5351_MULTISYNTH_MAX_FREQ * SI5351_FREQ_MULT) {
        freq = SI5351_MULTISYNTH_MAX_FREQ * SI5351_FREQ_MULT;
    }
    if (freq < SI5351_MULTISYNTH_MIN_FREQ * SI5351_FREQ_MULT) {
        freq = SI5351_MULTISYNTH_MIN_FREQ * SI5351_FREQ_MULT;
    }

    divby4 = 0;
    if (freq >= SI5351_MULTISYNTH_DIVBY4_FREQ * SI5351_FREQ_MULT) {
        divby4 = 1;
    }

    if(pll_freq == 0) {
        // Find largest integer divider for max
        // VCO frequency and given target frequency
        if(divby4 == 0) {
            lltmp = SI5351_PLL_VCO_MAX * SI5351_FREQ_MULT;
            do_div(lltmp, freq);
            a = (uint32_t)lltmp;
        } else {
            a = 4;
        }

        b = 0;
        c = 1;
        pll_freq = a * freq;
    } else {
        // Preset PLL, so return the actual freq for these params instead of PLL freq
        ret_val = 1;
        // Determine integer part of feedback equation
        a = pll_freq / freq;
        if (a < SI5351_MULTISYNTH_A_MIN) {
            freq = pll_freq / SI5351_MULTISYNTH_A_MIN;
        }
        if (a > SI5351_MULTISYNTH_A_MAX) {
            freq = pll_freq / SI5351_MULTISYNTH_A_MAX;
        }
        b = (pll_freq % freq * RFRAC_DENOM) / freq;
        c = b ? RFRAC_DENOM : 1;
    }

    // Calculate parameters
    if (divby4 == 1) {
        p3 = 1;
        p2 = 0;
        p1 = 0;
    } else {
        p1 = 128 * a + ((128 * b) / c) - 512;
        p2 = 128 * b - c * ((128 * b) / c);
        p3 = c;
    }

    reg->p1 = p1;
    reg->p2 = p2;
    reg->p3 = p3;

    if(ret_val == 0) {
        return pll_freq;
    } else {
        return freq;
    }
}

uint64_t si5351_pll_calc(uint64_t freq, struct Si5351RegSet *reg, int32_t correction) {
    uint64_t ref_freq = Si5351_Config.xtal_freq * SI5351_FREQ_MULT;
    uint32_t a, b, c, p1, p2, p3;
    uint64_t lltmp, denom;


    // Factor calibration value into nominal crystal frequency
    // Measured in parts-per-billion

    ref_freq = ref_freq + (int32_t)((((((int64_t)correction) << 31) / 1000000000LL) * ref_freq) >> 31);

    // PLL bounds checking
    if (freq < SI5351_PLL_VCO_MIN * SI5351_FREQ_MULT) {
        freq = SI5351_PLL_VCO_MIN * SI5351_FREQ_MULT;
    }
    if (freq > SI5351_PLL_VCO_MAX * SI5351_FREQ_MULT) {
        freq = SI5351_PLL_VCO_MAX * SI5351_FREQ_MULT;
    }

    // Determine integer part of feedback equation
    a = freq / ref_freq;

    if (a < SI5351_PLL_A_MIN) {
        freq = ref_freq * SI5351_PLL_A_MIN;
    }
    if (a > SI5351_PLL_A_MAX) {
        freq = ref_freq * SI5351_PLL_A_MAX;
    }

    // Find best approximation for b/c = fVCO mod fIN
    denom = 1000ULL * 1000ULL;
    lltmp = freq % ref_freq;
    lltmp *= denom;
    do_div(lltmp, ref_freq);

    b = (((uint64_t)(freq % ref_freq)) * RFRAC_DENOM) / ref_freq;
    c = b ? RFRAC_DENOM : 1;

	// Calculate parameters
    p1 = 128 * a + ((128 * b) / c) - 512;
    p2 = 128 * b - c * ((128 * b) / c);
    p3 = c;

    // Recalculate frequency as fIN * (a + b/c)
    lltmp  = ref_freq;
    lltmp *= b;
    do_div(lltmp, c);
    freq = lltmp;
    freq += ref_freq * a;

    reg->p1 = p1;
    reg->p2 = p2;
    reg->p3 = p3;

    return freq;
}

uint8_t si5351_select_r_div(uint64_t *freq) {
    uint8_t r_div = SI5351_OUTPUT_CLK_DIV_1;

    // Choose the correct R divider
    if((*freq >= SI5351_CLKOUT_MIN_FREQ * SI5351_FREQ_MULT) 
            && (*freq < SI5351_CLKOUT_MIN_FREQ * SI5351_FREQ_MULT * 2)) {
        r_div = SI5351_OUTPUT_CLK_DIV_128;
        *freq *= 128ULL;
    } else if ((*freq >= SI5351_CLKOUT_MIN_FREQ * SI5351_FREQ_MULT * 2)
            && (*freq < SI5351_CLKOUT_MIN_FREQ * SI5351_FREQ_MULT * 4)) {
        r_div = SI5351_OUTPUT_CLK_DIV_64;
        *freq *= 64ULL;
    } else if ((*freq >= SI5351_CLKOUT_MIN_FREQ * SI5351_FREQ_MULT * 4)
            && (*freq < SI5351_CLKOUT_MIN_FREQ * SI5351_FREQ_MULT * 8)) {
        r_div = SI5351_OUTPUT_CLK_DIV_32;
        *freq *= 32ULL;
    } else if ((*freq >= SI5351_CLKOUT_MIN_FREQ * SI5351_FREQ_MULT * 8)
            && (*freq < SI5351_CLKOUT_MIN_FREQ * SI5351_FREQ_MULT * 16)) {
        r_div = SI5351_OUTPUT_CLK_DIV_16;
        *freq *= 16ULL;
    } else if ((*freq >= SI5351_CLKOUT_MIN_FREQ * SI5351_FREQ_MULT * 16)
            && (*freq < SI5351_CLKOUT_MIN_FREQ * SI5351_FREQ_MULT * 32)) {
        r_div = SI5351_OUTPUT_CLK_DIV_8;
        *freq *= 8ULL;
    } else if ((*freq >= SI5351_CLKOUT_MIN_FREQ * SI5351_FREQ_MULT * 32)
            && (*freq < SI5351_CLKOUT_MIN_FREQ * SI5351_FREQ_MULT * 64)) {
        r_div = SI5351_OUTPUT_CLK_DIV_4;
        *freq *= 4ULL;
    } else if ((*freq >= SI5351_CLKOUT_MIN_FREQ * SI5351_FREQ_MULT * 64)
            && (*freq < SI5351_CLKOUT_MIN_FREQ * SI5351_FREQ_MULT * 128)) {
        r_div = SI5351_OUTPUT_CLK_DIV_2;
        *freq *= 2ULL;
    }

    return r_div;
}

