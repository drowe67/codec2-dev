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

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "freedv_api.h"
#include "reliable_text.h"
#include "golay23.h"

#define RELIABLE_TEXT_MAX_LENGTH (8)
#define RELIABLE_TEXT_CRC_LENGTH (2)
#define RELIABLE_TEXT_MAX_RAW_LENGTH (RELIABLE_TEXT_MAX_LENGTH + RELIABLE_TEXT_CRC_LENGTH)

/* Two bytes of text/CRC equal four bytes of Golay(23,12). */
#define RELIABLE_TEXT_BYTES_PER_ENCODED_SEGMENT (2) 
#define RELIABLE_TEXT_ENCODED_SEGMENT_LENGTH (4)

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

static void reliable_text_freedv_callback_rx(void *state, char chr)
{
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
        char decodedStr[RELIABLE_TEXT_MAX_RAW_LENGTH + 2];
        char fullRawStr[RELIABLE_TEXT_MAX_RAW_LENGTH + 2];
        char* decodedStrPtr = &decodedStr[0];
        memset(decodedStr, 0, RELIABLE_TEXT_MAX_RAW_LENGTH + 1);
        memset(fullRawStr, 0, RELIABLE_TEXT_MAX_RAW_LENGTH + 1);
        
        for (int charIndex = bufferIndex; (charIndex + RELIABLE_TEXT_ENCODED_SEGMENT_LENGTH - 1) < RELIABLE_TEXT_TOTAL_ENCODED_LENGTH + RELIABLE_TEXT_ENCODED_SEGMENT_LENGTH && strlen(fullRawStr) < RELIABLE_TEXT_MAX_RAW_LENGTH; charIndex += RELIABLE_TEXT_ENCODED_SEGMENT_LENGTH)
        {
            int encodedInput =
                ((obj->inbound_pending_chars[charIndex] & 0x3F) << 18) |
                ((obj->inbound_pending_chars[charIndex + 1] & 0x3F) << 12) |
                ((obj->inbound_pending_chars[charIndex + 2] & 0x3F) << 6) |
                (obj->inbound_pending_chars[charIndex + 3] & 0x3F);
            encodedInput &= 0x7FFFFF;
            
            int rawOutput = golay23_decode(encodedInput) >> 11;
            char rawStr[RELIABLE_TEXT_BYTES_PER_ENCODED_SEGMENT + 1];
            
            rawStr[0] = ((unsigned int)rawOutput) >> 6;
            rawStr[1] = ((unsigned int)rawOutput) & 0x3F;
            rawStr[2] = 0;
            strncat(fullRawStr, rawStr, RELIABLE_TEXT_BYTES_PER_ENCODED_SEGMENT);
            
            convert_ota_string_to_callsign_(rawStr, decodedStrPtr, RELIABLE_TEXT_BYTES_PER_ENCODED_SEGMENT);
            decodedStrPtr += RELIABLE_TEXT_BYTES_PER_ENCODED_SEGMENT;
        }
        
        // Decoded string must be at least more than one character to proceed.
        if (decodedStr[0] == 0 || decodedStr[1] == 0 || decodedStr[2] == 0) continue;
        
        // Get expected and actual CRC.
        unsigned char receivedCRC = convertHexStringToDigit_(&decodedStr[0]);
        unsigned char calcCRC = calculateCRC8_(&fullRawStr[RELIABLE_TEXT_CRC_LENGTH], strlen(&fullRawStr[RELIABLE_TEXT_CRC_LENGTH]));
        
        if (receivedCRC == calcCRC)
        {
            // We got a valid string. Call assigned callback.
            obj->text_rx_callback(&decodedStr[RELIABLE_TEXT_CRC_LENGTH], strlen(&decodedStr[RELIABLE_TEXT_CRC_LENGTH]));
            obj->has_successfully_decoded = 1;
        }
    }
}

static char reliable_text_freedv_callback_tx(void *state)
{
    reliable_text_impl_t* obj = (reliable_text_impl_t*)state;
    
    char ret = obj->tx_text[obj->tx_text_index];
    obj->tx_text_index = (obj->tx_text_index + 1) % (obj->tx_text_length + 1); // to ensure the null at the end is sent
    return ret;
}

reliable_text_t reliable_text_create()
{
    reliable_text_impl_t* ret = calloc(sizeof(reliable_text_impl_t), 1);
    ret->current_inbound_queue_pos = ret->inbound_pending_chars;
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
    impl->tx_text_length = raw_length * RELIABLE_TEXT_ENCODED_SEGMENT_LENGTH;
    impl->tx_text_index = 0;
    unsigned char crc = calculateCRC8_(&tmp[RELIABLE_TEXT_CRC_LENGTH], raw_length - RELIABLE_TEXT_CRC_LENGTH);
    unsigned char crcDigit1 = crc >> 4;
    unsigned char crcDigit2 = crc & 0xF;
    convertDigitToASCII_(&tmp[0], crcDigit1);
    convertDigitToASCII_(&tmp[1], crcDigit2);
    
    for(size_t index = 0, encodedIndex = 0; index < raw_length; index += RELIABLE_TEXT_BYTES_PER_ENCODED_SEGMENT, encodedIndex += RELIABLE_TEXT_ENCODED_SEGMENT_LENGTH)
    {
        // Encode the character as four bytes with parity bits.
        int inputRaw = ((tmp[index] & 0x3F) << 6) | (tmp[index+1] & 0x3F);
        unsigned int outputEncoding = golay23_encode(inputRaw) & 0x7FFFFF;
        impl->tx_text[encodedIndex] = (unsigned int)(outputEncoding >> 18) & 0x3F;
        impl->tx_text[encodedIndex + 1] = (unsigned int)(outputEncoding >> 12) & 0x3F;
        impl->tx_text[encodedIndex + 2] = (unsigned int)(outputEncoding >> 6) & 0x3F;
        impl->tx_text[encodedIndex + 3] = (unsigned int)outputEncoding & 0x3F;
    }
}

void reliable_text_use_with_freedv(reliable_text_t ptr, struct freedv* fdv, on_text_rx_t text_rx_fn)
{
    reliable_text_impl_t* impl = (reliable_text_impl_t*)ptr;
    
    impl->text_rx_callback = text_rx_fn;
    freedv_set_callback_txt(fdv, reliable_text_freedv_callback_rx, reliable_text_freedv_callback_tx, impl);
}
