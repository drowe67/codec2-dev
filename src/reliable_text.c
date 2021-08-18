//==========================================================================
// Name:            reliable_text.c
//
// Purpose:         Handles reliable text (e.g. text with FEC).
// Created:         August 15, 2020
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "freedv_api.h"
#include "reliable_text.h"
#include "ldpc_codes.h"

#define LDPC_TOTAL_SIZE_BITS (112)

#define RELIABLE_TEXT_MAX_LENGTH (8)
#define RELIABLE_TEXT_CRC_LENGTH (1)
#define RELIABLE_TEXT_MAX_RAW_LENGTH (RELIABLE_TEXT_MAX_LENGTH + RELIABLE_TEXT_CRC_LENGTH)

/* Two bytes of text/CRC equal four bytes of LDPC(112,56). */
#define RELIABLE_TEXT_BYTES_PER_ENCODED_SEGMENT (8) 
#define RELIABLE_TEXT_ENCODED_SEGMENT_LENGTH (16)

/* Size of (RELIABLE_TEXT_MAX_LENGTH + RELIABLE_TEXT_CRC_LENGTH) once encoded. */
#define RELIABLE_TEXT_NUMBER_SEGMENTS (RELIABLE_TEXT_MAX_RAW_LENGTH / RELIABLE_TEXT_BYTES_PER_ENCODED_SEGMENT)
#define RELIABLE_TEXT_TOTAL_ENCODED_LENGTH (RELIABLE_TEXT_NUMBER_SEGMENTS * RELIABLE_TEXT_ENCODED_SEGMENT_LENGTH)

/* Internal definition of reliable_text_t. */
typedef struct 
{
    on_text_rx_t text_rx_callback;
    
    char tx_text[RELIABLE_TEXT_TOTAL_ENCODED_LENGTH];
    int tx_text_index;
    int tx_text_length;
    
    char inbound_pending_chars[RELIABLE_TEXT_TOTAL_ENCODED_LENGTH + RELIABLE_TEXT_ENCODED_SEGMENT_LENGTH];
    char* current_inbound_queue_pos;
    int has_successfully_decoded;
    
    struct LDPC ldpc;
} reliable_text_impl_t;

// 6 bit character set for text field use:
// 0: ASCII null
// 1-9: ASCII 38-47
// 10-19: ASCII '0'-'9'
// 20-46: ASCII 'A'-'Z'
// 47: ASCII ' '
static void convert_callsign_to_ota_string_(const char* input, char* output, int maxLength)
{
    int outidx = 0;
    for (size_t index = 0; index < maxLength; index++)
    {
        if (input[index] == 0) break;
        
        if (input[index] >= 38 && input[index] <= 47)
        {
            output[outidx++] = input[index] - 37;
        }
        else if (input[index] >= '0' && input[index] <= '9')
        {
            output[outidx++] = input[index] - '0' + 10;
        }
        else if (input[index] >= 'A' && input[index] <= 'Z')
        {
            output[outidx++] = input[index] - 'A' + 20;
        }
        else if (input[index] >= 'a' && input[index] <= 'z')
        {
            output[outidx++] = toupper(input[index]) - 'A' + 20;
        }
    }
    output[outidx] = 0;
}

static void convert_ota_string_to_callsign_(const char* input, char* output, int maxLength)
{
    int outidx = 0;
    for (size_t index = 0; index < maxLength; index++)
    {
        if (input[index] == 0) break;
        
        if (input[index] >= 1 && input[index] <= 9)
        {
            output[outidx++] = input[index] + 37;
        }
        else if (input[index] >= 10 && input[index] <= 19)
        {
            output[outidx++] = input[index] - 10 + '0';
        }
        else if (input[index] >= 20 && input[index] <= 46)
        {
            output[outidx++] = input[index] - 20 + 'A';
        }
    }
    output[outidx] = 0;
}

static char calculateCRC8_(char* input, int length)
{
    unsigned char generator = 0x1D;
    unsigned char crc = 0; /* start with 0 so first byte can be 'xored' in */

    while (length > 0)
    {
        unsigned char ch = *input++;
        length--;

        // Ignore 6-bit carriage return and sync characters.
        if (ch == 63 || ch == 48) continue;
        
        crc ^= ch; /* XOR-in the next input byte */
        
        for (int i = 0; i < 8; i++)
        {
            if ((crc & 0x80) != 0)
            {
                crc = (unsigned char)((crc << 1) ^ generator);
            }
            else
            {
                crc <<= 1;
            }
        }
    }

    return crc;
}

static void convertDigitToASCII_(char* dest, unsigned char digit)
{
    if (digit >= 0 && digit <= 9)
    {
        *dest = digit + 10; // using 6 bit character set defined above.
    }
    else if (digit >= 0xA && digit <= 0xF)
    {
        *dest = (digit - 0xA) + 20; // using 6 bit character set defined above.
    }
    else
    {
        // Should not reach here.
        *dest = 10;
    }
}

static unsigned char convertHexStringToDigit_(char* src)
{
    unsigned char ret = 0;
    for (int i = 0; i < 2; i++)
    {
        ret <<= 4;
        unsigned char temp = 0;
        if (*src >= '0' && *src <= '9')
        {
            temp = *src - '0';
        }
        else if (*src >= 'A' && *src <= 'F')
        {
            temp = *src - 'A' + 0xA;
        }
        ret |= temp;
        src++;
    }
    
    return ret;
}

static int reliable_text_ldpc_decode(reliable_text_impl_t* obj, char* dest, char* src)
{
    float incomingData[LDPC_TOTAL_SIZE_BITS];
    float llr[LDPC_TOTAL_SIZE_BITS];
    char output[LDPC_TOTAL_SIZE_BITS];
    int parityCheckCount = 0;
    
    for (int bitIndex = 0; bitIndex < LDPC_TOTAL_SIZE_BITS; bitIndex++)
    {
        // Map to value expected by sd_to_llr().
        int bit = (src[bitIndex / 7] & (1 << (bitIndex % 7))) != 0;
        incomingData[bitIndex] = 1.0 - 2.0 * bit;
    }
    
    sd_to_llr(llr, incomingData, LDPC_TOTAL_SIZE_BITS);
    
    int iter = run_ldpc_decoder(&obj->ldpc, output, llr, &parityCheckCount);
    
    // Data is valid if BER < 0.1.
    float ber_est = (float)(obj->ldpc.NumberParityBits - parityCheckCount)/obj->ldpc.NumberParityBits;
    int result = (ber_est < 0.2);
        
    if (result)
    {        
        memset(dest, 0, RELIABLE_TEXT_BYTES_PER_ENCODED_SEGMENT);
        
        for (int bitIndex = 0; bitIndex < 8; bitIndex++)
        {
            if (output[bitIndex])
                dest[0] |= 1 << bitIndex;
        }
        for (int bitIndex = 8; bitIndex < (LDPC_TOTAL_SIZE_BITS / 2); bitIndex++)
        {
            int bitsSinceCrc = bitIndex - 8;
            if (output[bitIndex])
                dest[1 + (bitsSinceCrc / 6)] |= (1 << (bitsSinceCrc % 6));
        }
    }
    
    return result;
}

static void reliable_text_freedv_callback_rx(void *state, char chr)
{
    fprintf(stderr, "received char: %d\n", (chr & 0x3F));
    
    reliable_text_impl_t* obj = (reliable_text_impl_t*)state;
        
    // No need to further process if we got a valid string already.
    if (obj->has_successfully_decoded)
    {
        return;
    }
    
    // Queue up received character
    int increment_queue_pos = 1;
    if (obj->current_inbound_queue_pos == &obj->inbound_pending_chars[RELIABLE_TEXT_TOTAL_ENCODED_LENGTH + RELIABLE_TEXT_ENCODED_SEGMENT_LENGTH - 1])
    {
        // Shift all received characters left by 1.
        memmove(obj->inbound_pending_chars, &obj->inbound_pending_chars[1], RELIABLE_TEXT_TOTAL_ENCODED_LENGTH + RELIABLE_TEXT_ENCODED_SEGMENT_LENGTH - 1);
        increment_queue_pos = 0;
    }
    *obj->current_inbound_queue_pos = chr;
    if (increment_queue_pos) 
    {
        obj->current_inbound_queue_pos++;
    }
    
    for (int bufferIndex = 0; bufferIndex < RELIABLE_TEXT_ENCODED_SEGMENT_LENGTH; bufferIndex++)
    {
        char decodedStr[RELIABLE_TEXT_MAX_LENGTH + 1];
        char rawStr[RELIABLE_TEXT_BYTES_PER_ENCODED_SEGMENT + 1];
        memset(rawStr, 0, RELIABLE_TEXT_BYTES_PER_ENCODED_SEGMENT + 1);
        
        // Only proceed if BER is low enough.
        if (reliable_text_ldpc_decode(obj, rawStr, &obj->inbound_pending_chars[bufferIndex]) == 0) continue;
        
        convert_ota_string_to_callsign_(&rawStr[1], &decodedStr[1], RELIABLE_TEXT_BYTES_PER_ENCODED_SEGMENT);
        decodedStr[0] = rawStr[0]; // CRC
        
        // Decoded string must be at least more than one character to proceed.
        if (decodedStr[0] == 0 || decodedStr[1] == 0) continue;
        
        // Get expected and actual CRC.
        unsigned char receivedCRC = decodedStr[0];
        unsigned char calcCRC = calculateCRC8_(&rawStr[RELIABLE_TEXT_CRC_LENGTH], strlen(&rawStr[RELIABLE_TEXT_CRC_LENGTH]));
        
        if (receivedCRC == calcCRC)
        {
            // We got a valid string. Call assigned callback.
            obj->has_successfully_decoded = 1;
            obj->text_rx_callback(&decodedStr[RELIABLE_TEXT_CRC_LENGTH], strlen(&decodedStr[RELIABLE_TEXT_CRC_LENGTH]));
            
            // Resync so that byte 0 of the queue is the data immediately after the packet we found.
            memmove(
                &obj->inbound_pending_chars[0],
                &obj->inbound_pending_chars[bufferIndex + RELIABLE_TEXT_TOTAL_ENCODED_LENGTH], 
                RELIABLE_TEXT_ENCODED_SEGMENT_LENGTH  - bufferIndex);
            obj->current_inbound_queue_pos = &obj->inbound_pending_chars[RELIABLE_TEXT_ENCODED_SEGMENT_LENGTH - bufferIndex];
            break;
        }
    }
}

static char reliable_text_freedv_callback_tx(void *state)
{
    reliable_text_impl_t* obj = (reliable_text_impl_t*)state;
    
    char ret = obj->tx_text[obj->tx_text_index];
    obj->tx_text_index = (obj->tx_text_index + 1) % (obj->tx_text_length + 1); // to ensure the null at the end is sent
    
    fprintf(stderr, "sent char: %d\n", ret);
    return ret;
}

reliable_text_t reliable_text_create()
{
    reliable_text_impl_t* ret = calloc(sizeof(reliable_text_impl_t), 1);
    ret->current_inbound_queue_pos = ret->inbound_pending_chars;
    
    // Load LDPC code into memory.
    int code_index = ldpc_codes_find("HRA_56_56");
    memcpy(&ret->ldpc, &ldpc_codes[code_index], sizeof(struct LDPC));
    
    return (reliable_text_t)ret;
}

void reliable_text_destroy(reliable_text_t ptr)
{
    free(ptr);
}

void reliable_text_reset(reliable_text_t ptr)
{
    reliable_text_impl_t* impl = (reliable_text_impl_t*)ptr;
    impl->current_inbound_queue_pos = impl->inbound_pending_chars;
    impl->has_successfully_decoded = 0;
}

void reliable_text_set_string(reliable_text_t ptr, const char* str, int strlength)
{
    reliable_text_impl_t* impl = (reliable_text_impl_t*)ptr;
    
    char tmp[RELIABLE_TEXT_MAX_RAW_LENGTH];
    
    convert_callsign_to_ota_string_(str, &tmp[RELIABLE_TEXT_CRC_LENGTH], strlength);
    
    int raw_length = strlen(&tmp[RELIABLE_TEXT_CRC_LENGTH]) + RELIABLE_TEXT_CRC_LENGTH;
    impl->tx_text_length = RELIABLE_TEXT_ENCODED_SEGMENT_LENGTH;
    impl->tx_text_index = 0;
    unsigned char crc = calculateCRC8_(&tmp[RELIABLE_TEXT_CRC_LENGTH], raw_length - RELIABLE_TEXT_CRC_LENGTH);
    tmp[0] = crc;

    // Encode block of text using LDPC(112,56).
    char ibits[LDPC_TOTAL_SIZE_BITS / 2];
    char pbits[LDPC_TOTAL_SIZE_BITS / 2];
    memset(ibits, 0, LDPC_TOTAL_SIZE_BITS / 2);
    memset(pbits, 0, LDPC_TOTAL_SIZE_BITS / 2);
    for (int index = 0; index < 8; index++)
    {
        if (tmp[0] & (1 << index)) ibits[index] = 1;
    }

    // Pack 6 bit characters into single LDPC block.
    for (int ibitsBitIndex = 8; ibitsBitIndex < (LDPC_TOTAL_SIZE_BITS / 2); ibitsBitIndex++)
    {
        int bitsFromCrc = ibitsBitIndex - 8;
        if (tmp[1 + bitsFromCrc / 6] & (1 << (bitsFromCrc % 6)))
            ibits[ibitsBitIndex] = 1;
    }
    
    encode(&impl->ldpc, ibits, pbits);  
    
    // Split LDPC encoded bits into 7 bit characters due to varicode size limit.
    memset(impl->tx_text, 0, RELIABLE_TEXT_ENCODED_SEGMENT_LENGTH);
    for(int ibitsBitIndex = 0; ibitsBitIndex < (LDPC_TOTAL_SIZE_BITS / 2); ibitsBitIndex++)
    {
        if (ibits[ibitsBitIndex])
        {
            impl->tx_text[ibitsBitIndex / 7] |= (1 << (ibitsBitIndex % 7));
        }
    }
    
    // Split parity bits into 7 bit characters as well and append to the end.
    for(int pbitsBitIndex = 0; pbitsBitIndex < (LDPC_TOTAL_SIZE_BITS / 2); pbitsBitIndex++)
    {
        if (pbits[pbitsBitIndex])
        {
            impl->tx_text[(LDPC_TOTAL_SIZE_BITS / 2) / 7 + (pbitsBitIndex / 7)] |= (1 << (pbitsBitIndex % 7));
        }
    }
}

void reliable_text_use_with_freedv(reliable_text_t ptr, struct freedv* fdv, on_text_rx_t text_rx_fn)
{
    reliable_text_impl_t* impl = (reliable_text_impl_t*)ptr;
    
    impl->text_rx_callback = text_rx_fn;
    freedv_set_callback_txt(fdv, reliable_text_freedv_callback_rx, reliable_text_freedv_callback_tx, impl);
}
