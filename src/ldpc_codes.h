/*
  FILE...: ldpc_codes.h
  AUTHOR.: David Rowe
  CREATED: July 2020

  Array of LDPC codes used for various Codec2 waveforms.
*/

#ifndef __LDPC_CODES__
#include "mpdecode_core.h"

extern struct LDPC ldpc_codes[];
int ldpc_codes_num(void);
void ldpc_codes_list();
int ldpc_codes_find(char name[]);
void ldpc_codes_setup(struct LDPC *ldpc, char name[]);

#endif
