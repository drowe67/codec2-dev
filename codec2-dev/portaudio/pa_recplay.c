/* 
   pa_recplay.c
   David Rowe
   July 8 2012

   Echos audio from sound card input to sound card output. Samples at
   48 kHz, converts to 8 kHz, converts back to 48kHz, and plays using
   the default sound device.  Used as an intermediate step in
   Portaudio integration.

   Modified from paex_record.c Portaudio example. Original author
   author Phil Burk http://www.softsynth.com
*/

/*
 * $Id: paex_record.c 1752 2011-09-08 03:21:55Z philburk $
 *
 * This program uses the PortAudio Portable Audio Library.
 * For more information see: http://www.portaudio.com
 * Copyright (c) 1999-2000 Ross Bencina and Phil Burk
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files
 * (the "Software"), to deal in the Software without restriction,
 * including without limitation the rights to use, copy, modify, merge,
 * publish, distribute, sublicense, and/or sell copies of the Software,
 * and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
 * ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "portaudio.h"
#include "fdmdv.h"

#define SAMPLE_RATE  48000         /* 48 kHz sampling rate rec. as we
				      can trust accuracy of sound
				      card */
#define N8           160           /* processing buffer size at 8 kHz */
#define N48          (N8*FDMDV_OS) /* processing buffer size at 48 kHz */
#define MEM8 (FDMDV_OS_TAPS/FDMDV_OS)
#define NUM_CHANNELS 2             /* I think most sound cards prefer
				      stereo, we will convert to mono
				      as we sample */

/* state information passed to call back */

typedef struct {
    float               in48k[FDMDV_OS_TAPS + N48];
    float               in8k[MEM8 + N8];
} paTestData;


/* 
   This routine will be called by the PortAudio engine when audio is
   required.  It may be called at interrupt level on some machines so
   don't do anything that could mess up the system like calling
   malloc() or free().
*/

static int callback( const void *inputBuffer, void *outputBuffer,
                           unsigned long framesPerBuffer,
                           const PaStreamCallbackTimeInfo* timeInfo,
                           PaStreamCallbackFlags statusFlags,
                           void *userData )
{
    paTestData *data = (paTestData*)userData;
    int         i;
    short      *rptr = (short*)inputBuffer;
    short      *wptr = (short*)outputBuffer;
    float      *in8k = data->in8k;
    float      *in48k = data->in48k;
    float       out8k[N8];
    float       out48k[N48];
    short       out48k_short[N48];

    (void) timeInfo;
    (void) statusFlags;

    assert(inputBuffer != NULL);

    /* just use left channel */

    for(i=0; i<framesPerBuffer; i++,rptr+=2)
	data->in48k[i+FDMDV_OS_TAPS] = *rptr; 

    /* downsample and update filter memory */

    fdmdv_48_to_8(out8k, &in48k[FDMDV_OS_TAPS], N8);
    for(i=0; i<FDMDV_OS_TAPS; i++)
	in48k[i] = in48k[i+framesPerBuffer];

    /* play side, back up to 8k */

    for(i=0; i<N8; i++)
	in8k[MEM8+i] = out8k[i];

    /* upsample and update filter memory */

    fdmdv_8_to_48(out48k, &in8k[MEM8], N8);
    for(i=0; i<MEM8; i++)
	in8k[i] = in8k[i+N8];

    assert(outputBuffer != NULL);

    /* write signal to both channels */

    for(i=0; i<N48; i++)
	out48k_short[i] = (short)out48k[i];
    for(i=0; i<framesPerBuffer; i++,wptr+=2) {
	wptr[0] = out48k_short[i]; 
	wptr[1] = out48k_short[i]; 
    }

    return paContinue;
}

int main(int argc, char *argv[])
{
    PaStreamParameters  inputParameters, outputParameters;
    PaStream*           stream;
    PaError             err = paNoError;
    paTestData          data;
    int                 i;

    for(i=0; i<MEM8; i++)
	data.in8k[i] = 0.0;
    for(i=0; i<FDMDV_OS_TAPS; i++)
	data.in48k[i] = 0.0;

    err = Pa_Initialize();
    if( err != paNoError ) goto done;

    inputParameters.device = Pa_GetDefaultInputDevice(); /* default input device */
    if (inputParameters.device == paNoDevice) {
        fprintf(stderr,"Error: No default input device.\n");
        goto done;
    }
    inputParameters.channelCount = NUM_CHANNELS;         /* stereo input */
    inputParameters.sampleFormat = paInt16;
    inputParameters.suggestedLatency = Pa_GetDeviceInfo( inputParameters.device )->defaultLowInputLatency;
    inputParameters.hostApiSpecificStreamInfo = NULL;

    outputParameters.device = Pa_GetDefaultOutputDevice(); /* default output device */
    if (outputParameters.device == paNoDevice) {
        fprintf(stderr,"Error: No default output device.\n");
        goto done;
    }
    outputParameters.channelCount = NUM_CHANNELS;         /* stereo output */
    outputParameters.sampleFormat = paInt16;
    outputParameters.suggestedLatency = Pa_GetDeviceInfo( outputParameters.device )->defaultLowOutputLatency;
    outputParameters.hostApiSpecificStreamInfo = NULL;

    /* Play some audio --------------------------------------------- */

    err = Pa_OpenStream(
              &stream,
	      &inputParameters,
              &outputParameters,
              SAMPLE_RATE,
              N48,
              paClipOff,      
              callback,
              &data );
    if( err != paNoError ) goto done;

    err = Pa_StartStream( stream );
    if( err != paNoError ) goto done;

    while( ( err = Pa_IsStreamActive( stream ) ) == 1 )
    {
        Pa_Sleep(100);
    }
    if( err < 0 ) goto done;

    err = Pa_CloseStream( stream );
    if( err != paNoError ) goto done;


done:
    Pa_Terminate();
    if( err != paNoError )
    {
        fprintf( stderr, "An error occured while using the portaudio stream\n" );
        fprintf( stderr, "Error number: %d\n", err );
        fprintf( stderr, "Error message: %s\n", Pa_GetErrorText( err ) );
        err = 1;          /* Always return 0 or 1, but no other return codes. */
    }
    return err;
}

