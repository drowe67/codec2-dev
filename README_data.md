# README_data.md

# Introduction

FreeDV can be used to send data over radio channels.  Two APis are supported:
+ VHF packet data channel which uses Ethernet style framing.
+ Raw frames of modem data

## Credits

The VHF data channel was developed by Jeroen Vreeken.

## Quickstart

VHF packet data API:

1. Simple test using mode 2400A and VHF packet data

   ```sh
   $ cd ~/codec2/build_linux
   $ ./src/freedv_data_tx 2400A - --frames 15 | ./src/freedv_data_rx 2400A -  

   ```
   You can listen to the modem signal using:
   ```sh
   $ ./src/freedv_data_tx 2400A - --frames 15 | aplay -f S16_LE -r 48000
   
   ```
  
2. Same for 2400B and 800XA

   ```sh
   $ ./src/freedv_data_tx 2400B - --frames 15 | ./src/freedv_data_rx 2400B -  
   $ ./src/freedv_data_tx 800XA - --frames 15 | ./src/freedv_data_rx 800XA -  

   ```

3. Using a different callsign and secondary station id

   ```sh
   $ ./src/freedv_data_tx 2400A - --callsign T3ST --ssid 15 --frames 15 | src/freedv_data_rx 2400A -  
   ```

Raw modem frame API:

1. Let's send two modem frames of 14 bytes using FreeDV 700D:
   ```sh
   $ head -c 28 </dev/urandom > binaryIn.bin
   $ ./src/freedv_data_raw_tx 700D binaryIn.bin - |  ./src/freedv_data_raw_rx 700D - - | hexdump
   bits_per_modem_frame: 112 bytes_per_modem_frame: 14
   frames processed: 4  output bytes: 28
   $ hexdump binaryIn.bin 
   0000000 4325 0363 ce1f fb88 8102 7d76 c487 e092
   0000010 2ded bc06 7689 eb67 5dfe 43df          
   ```

# VHF Packet Data Channel

The FreeDV VHF data channel operates on a packet level. The FreeDV modems however typically operate on a fixed frame base. This means that data packets have to be sent in multiple frames.

The packet format is modeled after Ethernet. As a result, any protocol that is compatible with Ethernet can potentially be used over a FreeDV data link. (There are of course practical limits. Browsing the world wide web with just a few hundred bits per second will not be a pleasant experience.)

## Header optimization

When there are no packets available for transmission a small 'filler' packet with just the sender's address will be sent.
When there is a packet available not all of the header needs to be sent. The sender's address can often be left out if it was already sent in a previous frame. Likewise when the packet has no specific destination but is targeted at a multicast address, this can also be transmitted in a single bit as opposed to a 6 byte broadcast address.


## Addressing

Since the format is based on Ethernet, a 6 byte sender and destination address is used. It is possible to encode an ITU compatible callsign in these bytes. See http://dmlinking.net/eth_ar.html for more info. Or have a look at freedv_data_tx.c and freedv_data_rx.c for an actual implementation.

## Packet types

The 2 byte EtherType field is used to distinguish between various protocols.

## Checks

Not all channels are perfect, and especially since a packet is split up over multiple frames, bits might get lost. Each packet therefore has a CRC which is checked before it is accepted.  Note there is No FEC on 2400A/2400B/800XA.

## Available modes

The data channel is available for modes 2400A, 2400B and 800XA.

## API

The data channel is part of the regular FreeDV API.

### Initialization

After creating a new freedv instance with freedv_open(), a few more calls need to be done before the data channel is usable.

  ```
   void freedv_set_data_header             (struct freedv *freedv, unsigned char *header);
  ```

The address that will be used for 'filler' packets must be set. The freedv_set_data_header() function must be called with a 6 byte header.

  ```
   typedef void (*freedv_callback_datarx)(void *, unsigned char *packet, size_t size);
   typedef void (*freedv_callback_datatx)(void *, unsigned char *packet, size_t *size);
   void freedv_set_callback_data           (struct freedv *freedv, freedv_callback_datarx datarx, freedv_callback_datatx datatx, void *callback_state);
  ```

Using freedv_set_callback_data() two callback functions can be provided. The datarx callback will be used whenever a new data packet has been successfully received. The datatx callback will be used when a new data packet is required for transmission.

### Operation

  ```
   void freedv_datatx  (struct freedv *f, short mod_out[]);
  ```

During normal operation the freedv_datatx() function can be used whenever a data frame has to be sent. If no data is available it will request new data using the datatx callback. The callback function is allowed to set 'size' to zero if no data is available or if it wishes to send an address frame.

For reception the regular freedv_rx() functions can be used as received data will automatically be reported using the datarx callback. Be aware that these functions return the actual number of received speech samples. When a data frame is received the return value will be zero. This may lead to 'gaps' in the audio stream which will have to be filled with silence.

### Examples

The freedv_data_tx and freedv_data_rx test programs implement the minimum needed to send and receive data packets.

## Mixing voice and data

Encoding only voice data is easy with the FreeDV API. Simply use the freedv_tx() function and provide it with speech samples.
Likewise encoding only data is also easy. Make sure to provide a source of data frames using the freedv_set_callback_data() function, and use the freedv_datatx() function to generate frames.

However there are many use cases where one would like to transmit a mix of voice and data. For example one might want to transmit their callsign in a machine readable format, or a short position report. There are a few ways to do this:

### Data bursts at start and/or end of transmission

This method simply transmits voice frames during the transmission, except for a few moments. For example when the user keys the radio the software uses the freedv_datatx() function for a number of frames before switching to regular voice frames.
Likewise when the user releases the key the software may hold it for a number of frames to transmit data before it releases the actual radio.

Be careful though: depending on your setup (radio, PC, soundcard, etc) the generated frames and the keying of your radio might not be perfectly in sync and the first or last frames might be lost in the actual transmission. Make sure to take this into account when using this method.

### Data and voice interleaved

Another method is to generate a mixed stream of frames. Compared to a small burst at the beginning or end a lot more data can be sent. We only need a way to choose between voice or data such that the recovered speech at the other side is not impacted.

#### Detect voice activity

When it is possible to determine activity in the voice signal (and it almost always is) this presence can be used to insert a data frame by calling freedv_datatx() instead of freedv_tx()/freedv_codectx(). This method is used in the freedv_mixed_tx demo program. When the option --codectx is given the codec2 library is used to determine the activity.

  ```
   $ ./src/freedv_mixed_tx 2400A ../raw/hts1a.raw - --codectx | src/freedv_data_rx 2400A -
   $ ./src/freedv_mixed_tx 2400A ../raw/hts1a.raw - | src/freedv_data_rx 2400A -
  ```

The advantage of this method is that the audio is not distorted, there was nothing (or near nothing) to distort. A drawback is that constant voice activity may mean there are insufficient frames for data.

### Receiving mixed voice and data

Receiving and decoding a mixed voice and data stream is (almost) as easy as receiving a regular voice-only transmission.
One simply uses the regular API calls for reception of speech samples. In addition, the callback functions are used for data.
There is one caveat though: when a data frame is received the API functions (like freedv_rx) will return zero as this is the amount of codec/voice data received.
For proper playback silence (or comfort noise) should be inserted for the duration of a frame to restore the timing of the original source speech samples.
An example of how this is done is provided in freedv_mixed_rx

  ```
   $ ./src/freedv_mixed_tx 2400A ../raw/hts1a.raw - | src/freedv_mixed_rx 2400A - ./hts1a_out.raw
  ```

### Insert a data frame periodically

This is a very simple method, simply insert a data frame every n frames, (e.g. once every 10 seconds). Since single FreeDV frames are relatively short (tens of milliseconds) the effect on received audio will be minor. The advantage of this method is that one can create a guaranteed amount of data bandwidth. A drawback is some interruption in the audio that may be noticed.

### Combination of the above.

A combination of the two methods may also be used. Send data when no voice is active and insert a frame when this does not occur for a long time.

# Raw Data using the FreeDV API

The demo programs [freedv_data_raw_tx.c](src/freedv_data_raw_tx.c) and [freedv_data_raw_rx.c](src/freedv_data_raw_rx.c) show how to use the raw data API.

The following FreeDV modes are recommended for *preliminary development* using the raw data API.  These modes were originally designed for streaming voice rather than data and are not suitable for production HF data applications.  They have small payloads, and acquisition algorithms not suitable for packet data over real world HF channels.  New modes are being design for HF data at the time of writing (June 2020).

| FreeDV Mode | RF bandwidth (Hz) | Payload data rate bits/s | bytes/frame | FEC | Min SNR (dB, AWGN) |
| :-: | :-: | :-: | :-: | :-: | :-: |
| 700C | 1100 | 700 | 7 | none | 2 |
| 700D | 1100 | 700 | 14 | rate 0.8 | -2 |
| 2020 | 1500 | 1733 | 39 | rate 0.6 | 2 |

The API can signal if there were uncorrected bit errors in the frame returned using ```freedv_get_uncorrected_errors(freedv)```, these frames are usually discarded in data applications.

The raw data API may lose frames due to channel impairments, loss of sync, or acquisition delays.  The caller must handle these situations.  The caller is also responsible for segmentation/re-assembly of the modem frames into larger blocks of data.
