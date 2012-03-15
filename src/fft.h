/*---------------------------------------------------------------------------*\
  
  FILE........: fft.h
  AUTHOR......: Bruce Robertson
  DATE CREATED: 29/11/2010

  Bridge between existing code and kiss_fft.

\*---------------------------------------------------------------------------*/

#ifndef __FFT__
#define __FFT__
void fft(float x[], int n, int isign);

void cleanup_fft(void)
#ifdef __GNUC__
    __attribute__ ((destructor))
#endif
    ;

#endif	/* __FFT__ */

