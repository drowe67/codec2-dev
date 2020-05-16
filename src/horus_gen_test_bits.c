/*---------------------------------------------------------------------------*\

  FILE........: horus_gen_tx_bits.c
  AUTHOR......: Mark Jessop
  DATE CREATED: May 2020

  Horus dummy packet generation, for use with fsk_demod.

  Build:
  gcc horus_gen_test_bits.c horus_l2.c golay23.c -o horus_get_test_bits -Wall -DSCRAMBLER -DINTERLEAVER

  \*---------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <getopt.h>

#include "horus_l2.h"

// TODO: Move these packet format definitions to somehwere common.

/* Horus Mode 0 (Legacy 22-byte) Binary Packet */
struct TBinaryPacket0
{
    uint8_t     PayloadID;
    uint16_t	Counter;
    uint8_t	Hours;
    uint8_t	Minutes;
    uint8_t	Seconds;
    float	Latitude;
    float	Longitude;
    uint16_t  	Altitude;
    uint8_t     Speed;       // Speed in Knots (1-255 knots)
    uint8_t     Sats;
    int8_t      Temp;        // Twos Complement Temp value.
    uint8_t     BattVoltage; // 0 = 0.5v, 255 = 2.0V, linear steps in-between.
    uint16_t    Checksum;    // CRC16-CCITT Checksum.
}  __attribute__ ((packed));

/* Horus Mode 1 (32-byte) Binary Packet */
struct TBinaryPacket1
{
    uint16_t     PayloadID;
    uint16_t	Counter;
    uint8_t	Hours;
    uint8_t	Minutes;
    uint8_t	Seconds;
    float	Latitude;
    float	Longitude;
    uint16_t  	Altitude;
    uint8_t     Speed;       // Speed in Knots (1-255 knots)
    uint8_t     Sats;
    int8_t      Temp;        // Twos Complement Temp value.
    uint8_t     BattVoltage; // 0 = 0.5v, 255 = 2.0V, linear steps in-between.
    uint8_t     dummy1;      // Dummy values for user-configurable section.
    uint8_t     dummy2;
    uint8_t     dummy3;
    uint8_t     dummy4;
    uint8_t     dummy5;
    uint8_t     dummy6;
    uint8_t     dummy7;
    uint8_t     dummy8;
    uint8_t     dummy9;
    uint16_t    Checksum;    // CRC16-CCITT Checksum.
}  __attribute__ ((packed));

/* Horus Mode 2 (16-byte) Binary Packet */
struct TBinaryPacket2
{
    uint8_t     PayloadID;
    uint8_t	Counter;
    uint16_t	BiSeconds;
    uint8_t	  LatitudeMSB;
    uint16_t	Latitude;
    uint8_t	  LongitudeMSB;
    uint16_t	Longitude;
    uint16_t  	Altitude;
    uint8_t     BattVoltage; // 0 = 0.5v, 255 = 2.0V, linear steps in-between.
    uint8_t     flags;      // Dummy values for user-configurable section.
    uint16_t    Checksum;    // CRC16-CCITT Checksum.
}  __attribute__ ((packed));




int main(int argc,char *argv[]) {
    int i, framecnt;
    int horus_mode = 0;

    char usage[] = "usage: %s horus_mode numFrames\nMode 0 = Legacy 22-byte Golay FEC\nMode 1 = 32-byte LDPC FEC\nMode 2 = 16-byte LDPC FEC\n";

    if (argc < 3) {
        fprintf(stderr, usage, argv[0]);
        exit(1);
    }

    horus_mode = atoi(argv[1]);
    fprintf(stderr, "Using Horus Mode %d.\n", horus_mode);

    framecnt = atoi(argv[2]);
    fprintf(stderr, "Generating %d frames.\n", framecnt);

    if(horus_mode == 0){
      int nbytes = sizeof(struct TBinaryPacket0);
      struct TBinaryPacket0 input_payload;
      int num_tx_data_bytes = horus_l2_get_num_tx_data_bytes(nbytes);
      unsigned char tx[num_tx_data_bytes];

      /* all zeros is nastiest sequence for demod before scrambling */

      memset(&input_payload, 0, nbytes);
      input_payload.Checksum = horus_l2_gen_crc16((unsigned char*)&input_payload, nbytes-2);

      horus_l2_encode_tx_packet(tx, (unsigned char*)&input_payload, nbytes);

      int b;
      uint8_t tx_bit;
      while(framecnt > 0){
          for(i=0; i<num_tx_data_bytes; i++) {
              for(b=0; b<8; b++) {
                  tx_bit = (tx[i] >> (7-b)) & 0x1; /* msb first */
                  fwrite(&tx_bit,sizeof(uint8_t),1,stdout);
                  fflush(stdout);
              }
          }
          framecnt -= 1;
      }

    } else if(horus_mode == 1){
      // 32-Byte LDPC Encoded mode.
      int nbytes = sizeof(struct TBinaryPacket1);
      struct TBinaryPacket1 input_payload;

      // TODO: Add Calculation of expected number of TX bytes based on LDPC code.
      int num_tx_data_bytes = nbytes;
      unsigned char tx[num_tx_data_bytes];

      /* all zeros is nastiest sequence for demod before scrambling */
      memset(&input_payload, 0, nbytes);
      input_payload.Checksum = horus_l2_gen_crc16((unsigned char*)&input_payload, nbytes-2);


      // TODO: Replaced with LDPC Encoding
      memcpy(tx, (unsigned char*)&input_payload, nbytes);

      int b;
      uint8_t tx_bit;
      while(framecnt > 0){
          for(i=0; i<num_tx_data_bytes; i++) {
              for(b=0; b<8; b++) {
                  tx_bit = (tx[i] >> (7-b)) & 0x1; /* msb first */
                  fwrite(&tx_bit,sizeof(uint8_t),1,stdout);
                  fflush(stdout);
              }
          }
          framecnt -= 1;
      }
    } else if(horus_mode == 2){
      // 16-Byte LDPC Encoded mode.
      int nbytes = sizeof(struct TBinaryPacket2);
      struct TBinaryPacket2 input_payload;

      // TODO: Add Calculation of expected number of TX bytes based on LDPC code.
      int num_tx_data_bytes = nbytes;
      unsigned char tx[num_tx_data_bytes];

      /* all zeros is nastiest sequence for demod before scrambling */
      memset(&input_payload, 0, nbytes);
      input_payload.Checksum = horus_l2_gen_crc16((unsigned char*)&input_payload, nbytes-2);


      // TODO: Replaced with LDPC Encoding
      memcpy(tx, (unsigned char*)&input_payload, nbytes);

      int b;
      uint8_t tx_bit;
      while(framecnt > 0){
          for(i=0; i<num_tx_data_bytes; i++) {
              for(b=0; b<8; b++) {
                  tx_bit = (tx[i] >> (7-b)) & 0x1; /* msb first */
                  fwrite(&tx_bit,sizeof(uint8_t),1,stdout);
                  fflush(stdout);
              }
          }
          framecnt -= 1;
      }
    } else {
      fprintf(stderr, "Unknown Mode!");
    }

    return 0;
}