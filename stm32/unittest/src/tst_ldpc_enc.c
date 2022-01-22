/*
  FILE...: ldpc_enc.c
  AUTHOR.: Bill Cowley, David Rowe
  CREATED: Sep 2016

  STM32 Version: Aug 2018 - Don Reid
*/

/* This is a unit test implementation of the LDPC encode function.
 *
 * Typical run:

    ofdm_gen_test_bits stm_in.raw 6 --rand --ldpc

    ldpc_enc stm_in.raw ref_out.raw --code HRA_112_112

    <Load stm32 and run>

    cmp -l ref_out.raw stm_out.raw

 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>

#include "mpdecode_core.h"

#include "semihosting.h"
#include "stm32f4xx_conf.h"
#include "stm32f4xx.h"
#include "machdep.h"
#include "ldpc_codes.h"

static __attribute__ ((section (".ccm"))) char fin_buffer[8*8192];
char fout_buffer[1024];

int main(int argc, char *argv[])
{
    semihosting_init();

    printf("LDPC encode test and profile\n");

    PROFILE_VAR(ldpc_encode);

    machdep_profile_init();

    char code_name[255] = "H_2064_516_sparse";
    struct LDPC   ldpc;
    ldpc_codes_setup(&ldpc, code_name);
    int num_data_bits = ldpc.NumberRowsHcols;
    int num_parity_bits = ldpc.NumberParityBits;
    unsigned char ibits[num_data_bits];
    unsigned char pbits[num_parity_bits];
    
    FILE* fin = fopen("stm_in.raw", "rb");
    if (fin == NULL) {
        printf("Error opening input file\n");
        exit(1);
    }
    setvbuf(fin, fin_buffer,_IOFBF,sizeof(fin_buffer));

    FILE* fout = fopen("stm_out.raw", "wb");
    if (fout == NULL) {
        printf("Error opening output file\n");
        exit(1);
    }

    while (fread(ibits, sizeof(char), num_data_bits, fin) == num_data_bits) {

        PROFILE_SAMPLE(ldpc_encode);
        encode(&ldpc, ibits, pbits);
        PROFILE_SAMPLE_AND_LOG2(ldpc_encode, "  ldpc_encode");

        fwrite(ibits, sizeof(char) , num_data_bits, fout);
        fwrite(pbits, sizeof(char) , num_parity_bits, fout);
    }

    fclose(fin);
    fclose(fout);
    
    fflush(stdout);
    stdout = freopen("stm_profile", "w", stdout);
    machdep_profile_print_logged_samples();

    fclose(stdout);
    fclose(stderr);

    return 0;
}
