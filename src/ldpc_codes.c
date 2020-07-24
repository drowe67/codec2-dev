/*
  FILE...: ldpc_codes.c
  AUTHOR.: David Rowe
  CREATED: July 2020

  Array of LDPC codes used for various Codec2 waveforms.
*/

#include <stdio.h>
#include <string.h>
#include "ldpc_codes.h"
#include "H2064_516_sparse.h"  
#include "HRA_112_112.h"  
#include "HRAb_396_504.h"
#include "H_256_768_22.h"
#include "H_256_512_4.h"
#include "HRAa_1536_512.h"

struct LDPC ldpc_codes[] = {

    /* default Wenet High Alitiude Balloon rate 0.8 code */
    {
        "H2064_516_sparse",
        MAX_ITER,
        0,
        1,
        1,
        CODELENGTH,
        NUMBERPARITYBITS,
        NUMBERROWSHCOLS,
        MAX_ROW_WEIGHT,
        MAX_COL_WEIGHT,
        H_rows,
        H_cols
    },

    /* short rate 1/2 code for FreeDV 700D */
    {
        "HRA_112_112",
        HRA_112_112_MAX_ITER,
        0,
        1,
        1,
        HRA_112_112_CODELENGTH,
        HRA_112_112_NUMBERPARITYBITS,
        HRA_112_112_NUMBERROWSHCOLS,
        HRA_112_112_MAX_ROW_WEIGHT,
        HRA_112_112_MAX_COL_WEIGHT,
        (uint16_t *)HRA_112_112_H_rows,
        (uint16_t *)HRA_112_112_H_cols
    },

    /* rate 0.8 code used for FreeDV 2020 */
    {
        "HRAb_396_504",
        HRAb_396_504_MAX_ITER,
        0,
        1,
        1,
        HRAb_396_504_CODELENGTH,
        HRAb_396_504_NUMBERPARITYBITS,
        HRAb_396_504_NUMBERROWSHCOLS,
        HRAb_396_504_MAX_ROW_WEIGHT,
        HRAb_396_504_MAX_COL_WEIGHT,
        (uint16_t *)HRAb_396_504_H_rows,
            (uint16_t *)HRAb_396_504_H_cols
    },

    /* rate 1/3 code, works at raw BER of 14% */
    {
        "H_256_768",
        H_256_768_22_MAX_ITER,
        0,
        1,
        1,
        H_256_768_22_CODELENGTH,
        H_256_768_22_NUMBERPARITYBITS,
        H_256_768_22_NUMBERROWSHCOLS,
        H_256_768_22_MAX_ROW_WEIGHT,
        H_256_768_22_MAX_COL_WEIGHT,
        (uint16_t *)H_256_768_22_H_rows,
        (uint16_t *)H_256_768_22_H_cols
    },

    /* used for 4FSK/LLR experiments */
    {
        "H_256_512_4",
        H_256_512_4_MAX_ITER,
        0,
        1,
        1,
        H_256_512_4_CODELENGTH,
        H_256_512_4_NUMBERPARITYBITS,
        H_256_512_4_NUMBERROWSHCOLS,
        H_256_512_4_MAX_ROW_WEIGHT,
        H_256_512_4_MAX_COL_WEIGHT,
        (uint16_t *)H_256_512_4_H_rows,
        (uint16_t *)H_256_512_4_H_cols
    },
    
    /* used for 4FSK/LLR experiments */
    {
        "HRAa_1536_512",
        HRAa_1536_512_MAX_ITER,
        0,
        1,
        1,
        HRAa_1536_512_CODELENGTH,
        HRAa_1536_512_NUMBERPARITYBITS,
        HRAa_1536_512_NUMBERROWSHCOLS,
        HRAa_1536_512_MAX_ROW_WEIGHT,
        HRAa_1536_512_MAX_COL_WEIGHT,
        (uint16_t *)HRAa_1536_512_H_rows,
        (uint16_t *)HRAa_1536_512_H_cols
    }
};

int ldpc_codes_num(void) { return sizeof(ldpc_codes)/sizeof(struct LDPC); }

void ldpc_codes_list() {
    fprintf(stderr, "\n");
    for(int c=0; c<ldpc_codes_num(); c++) {
        int n =  ldpc_codes[c].NumberRowsHcols + ldpc_codes[c].NumberParityBits;
        int k = ldpc_codes[c].NumberRowsHcols;
        float rate = (float)k/n;
        fprintf(stderr, "%-20s rate %3.2f (%d,%d) \n", ldpc_codes[c].name, rate, n, k);
    }
    fprintf(stderr, "\n");
}

int ldpc_codes_find(char name[]) {
    int code_index = 0;
    for(int c=0; c<ldpc_codes_num(); c++)
        if (strcmp(ldpc_codes[c].name, name) == 0)
            code_index = c;
    return code_index;
}

