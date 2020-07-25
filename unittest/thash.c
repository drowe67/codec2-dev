/*---------------------------------------------------------------------------*\

  FILE........: thash.c
  AUTHOR......: David Rowe
  DATE CREATED: July 2020

  Simple test program for freeDV API get hash function

\*---------------------------------------------------------------------------*/

#include <stdio.h>
#include "freedv_api.h"

#define MAX_STR 16

int main(void) { 
    char hash[MAX_STR];
    freedv_get_hash(hash, MAX_STR);
    printf("%s\n", hash);
    return 0;
}


