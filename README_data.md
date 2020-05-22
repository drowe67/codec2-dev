# README_data.md

## Introduction

Documentation for FreeDV data channel.

## Quickstart

1. Simple test using mode 2400A

  ```
   $ ./src/freedv_data_tx 2400A - --frames 15 | ./src/freedv_data_rx 2400A -  

  ```

2. Same for 2400B and 800XA

  ```
   $ ./src/freedv_data_tx 2400B - --frames 15 | ./src/freedv_data_rx 2400B -  

  ```

  ```
   $ ./src/freedv_data_tx 800XA - --frames 15 | ./src/freedv_data_rx 800XA -  

  ```

3. Using a different callsign and secondary station id

  ```
   $ ./src/freedv_data_tx 2400A - --callsign T3ST --ssid 15 --frames 15 | src/freedv_data_rx 2400A -  
  ```

## Packets

The FreeDV data channel operates on a packet level. The FreeDV modems however typically operate on a fixed frame base. This means that data packets have to be send in multiple frames.

The packet format is modelled after EtherNet. As a result any protocol that is compatible with EtherNet can potentially be used over a FreeDV data link. (There are of course practical limits. Browsing the world wide web with just a few hundred bits per second will not be a pleasant experience.)

### Header optimization

When there are no packets available for transmission small 'filler' packet with just the sender's address will be send.
When there is a packet available not all of the header needs to be send. The sender's address can often be left out if it was already send in a previous frame. Likewise when the packet has no specific destination but is targeted at a multicast address, this can also be transmitted in a single bit as oposed to a 6 byte broadcast address.

### Addressing

Since the format is based on EtherNet a 6 byte sender and destination address is used. It is possible to encode a ITU compatible callsign in these bytes. See http://dmlinking.net/eth_ar.html for more info. Or have a look at freedv_data_tx.c and freedv_data_rx.c for an actual implementation.

### Packet types

The 2 byte EtherType field is used to distinguish between various protocols.

### Checks

Not all channels are perfect and especially since a packet is split up over multiple frames bits might get lost. Each packet therefore has a CRC which is checked before it is accepted.

## Available modes

The data channel is available for modes 2400A, 2400B and 800Xa.

## API

The data channel is part of the regular FreeDV API

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

During normal operation the freedv_datatx() function can be used whenever a data frame has to be sent. If no data is available it will request new data using the datatx callback. The callback function is allowed to set 'size' to zero if no data is available or if it whishes to send an address frame.

For reception the regular freedv_rx() functions can be used as received data will automatically be reported using the datarx callback. Be aware that these functions return the actual number of received speech samples. When a data frame is received the return value will be zero. This may lead to 'gaps' in the audio stream which will have to be filled with silence.

### Examples

The freedv_data_tx and freedv_data_rx test programs implement the minimum needed to send and receive data packets.

## Mixing voice and data

Encoding only voice data is easy with the FreeDV API. Simply use the freedv_tx() function and provide it with speech samples.
Likewise encoding only data is also easy. Make sure to provide a source of data frames using the freedv_set_callback_data() function, and use the freedv_datatx() function to generate frames.

However there are many usecases where one would like to transmit a mix of voice and data. For example one might want to transmit their callsign in a machine readable format, or a short position report. There are a few ways to do this:

### Data bursts at start and/or end of transmission

This method simply transmits voice frames during the transmission, except for a few moments. For example when the user keys the radio the software uses the freedv_datatx() function for a number of frames before switching to regular voice frames.
Likewise when the user releases the key the software may hold it for a number of frames to transmit data before it releases the actual radio.

Be carefull though: depending on your setup (radio, pc, soundcard, etc) the generated frames and the keying of your radio might not be perfectly in sync and the first or last frames might be lost in the actual transmission. Make sure to take this into account when using this method.

### Data and voice interleaved

Another method is to generate a mixed stream of frames. Compared to a small burst at begin or end a lot more data can be send. We only need a way to choose between voice or data such that the recovered speech at the other side is not impacted.

#### Detect voice activity

When it is possible to determine activity in the voice signal (and it almost always is) this presence can be used to insert a data frame by calling freedv_datatx() instead of freedv_tx()/freedv_codectx(). This method is used in th freedv_tx demo program when the options --codectx and --datatx are given.

  ```
   $ ./src/freedv_tx 2400A ../raw/hts1a.raw - --codectx --datatx | src/freedv_data_rx 2400A -
  ```

The advantage of this method is that the audio is not distorted, there was nothing (or near nothing) to distort. A drawback is that constant voice activity may mean there are insufficient frames for data.

#### Insert a data frame periodically

This is a very simple method, simply insert a data frame every n frames, (e.g. once every 10 seconds). Since single FreeDV frames are relativly short (tens of milliseconds) the effect on received audio will be minor. The advantage of this method is that one can create a guaranteed amount of data bandwidth. A drawback is some interruption in the audio that may be noticed.

#### Combination of the above.

A combination of the two methods may also be used. Send data when no voice is active and insert a frame when this does not occur for a long time.

