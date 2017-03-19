/* 
  FILE...: ldpc_dec.c
  AUTHOR.: Matthew C. Valenti, Rohit Iyer Seshadri, David Rowe
  CREATED: Sep 2016

  Command line C LDPC decoder derived from MpDecode.c in the CML
  library.  Allows us to run the same decoder in Octave and C.  The
  code is defined by the parameters and array stored in the include
  file below, which can be machine generated from the Octave function
  ldpc_fsk_lib.m:ldpc_decode()

  The include file also contains test input/output vectors for the LDPC
  decoder for testing this program.

  Build:

    $ gcc -O2 -o ldpc_dec ldpc_dec.c mpdecode_core.c -Wall -lm -g

  Note: -O2 option was required to get identical results to MpDecode,
  which is also compiled with -O2.  Without it the number of bit errors
  between C and Octave was different, especially when the code did
  not converge and hit max_iters.

  TODO:
  [ ] C cmd line encoder
      [ ] SD output option
      [ ] Es/No option for testing
  [ ] decoder
      [X] test mode or file I/O (incl stdin/stdout)
      [X] Octave code to generate include file
          + MAX_ITER as well
      [X] check into SVN
      [ ] enc/dec running on cmd line
      [ ] fsk_demod modified for soft decisions
      [ ] drs232 modified for SD
          + use UW syn cin this program to check BER with coding
      [ ] revisit CML support, maybe blog post
*/

#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "mpdecode_core.h"

/* Machine generated consts, H_rows, H_cols, test input/output data to
   change LDPC code regenerate this file. */

#ifdef HRA_112_112
#include "HRA_112_112.h"  
#else
#include "H2064_516_sparse.h"  
#endif

int opt_exists(char *argv[], int argc, char opt[]) {
    int i;
    for (i=0; i<argc; i++) {
        if (strcmp(argv[i], opt) == 0) {
            return 1;
        }
    }
    return 0;
}

void extract_output(char out_char[], int DecodedBits[], int ParityCheckCount[], int max_iter, int CodeLength, int NumberParityBits);

int main(int argc, char *argv[])
{    
    int         CodeLength, NumberParityBits;
    int         i, r, num_ok, num_runs;
    char        out_char[CODELENGTH];
    struct LDPC ldpc;

    /* derive some parameters */

    CodeLength = CODELENGTH;                    /* length of entire codeword */
    NumberParityBits = NUMBERPARITYBITS;
	
    if (argc < 2) {
        fprintf(stderr, "\n");
        fprintf(stderr, "usage: %s --test\n\n", argv[0]);
        fprintf(stderr, "  Run internal self test and print code parameters.\n\n");
        fprintf(stderr, "usage: %s InOneSymbolPerDouble OutOneBitPerByte [--sd] [--half]\n\n", argv[0]);
        fprintf(stderr, "   InOneSymbolPerDouble    Input file of double LLRs, use - for the \n");        
        fprintf(stderr, "                           file names to use stdin/stdout\n");
        fprintf(stderr, "   --sd                    Treat input file samples as Soft Decision\n");
        fprintf(stderr, "                           demod outputs rather than LLRs\n");
        fprintf(stderr, "   --half                  Load framesize/2 input samples for each decode\n");
        fprintf(stderr, "                           attempt, only output decoded bits if decoder\n");
        fprintf(stderr, "                           converges.  Form of frame sync.\n");
        fprintf(stderr, "\n");
        exit(0);
    }

    /* set up LDPC code from include file constants */

    ldpc.max_iter = MAX_ITER;
    ldpc.dec_type = 0;
    ldpc.q_scale_factor = 1;
    ldpc.r_scale_factor = 1;
    ldpc.CodeLength = CODELENGTH;
    ldpc.NumberParityBits = NUMBERPARITYBITS;
    ldpc.NumberRowsHcols = NUMBERROWSHCOLS;
    ldpc.max_row_weight = MAX_ROW_WEIGHT;
    ldpc.max_col_weight = MAX_COL_WEIGHT;
    ldpc.H_rows = H_rows;
    ldpc.H_cols = H_cols;

    if (!strcmp(argv[1],"--test")) {

        /* test mode --------------------------------------------------------*/

        fprintf(stderr, "Starting test using pre-compiled test data .....\n");
        fprintf(stderr, "Codeword length: %d\n",  CodeLength);
        fprintf(stderr, "Parity Bits....: %d\n",  NumberParityBits);

        num_runs = 100; num_ok = 0;

        for(r=0; r<num_runs; r++) {

            run_ldpc_decoder(&ldpc, out_char, input);

            int ok = 0;
            for (i=0; i<CodeLength; i++) {
                if (out_char[i] == detected_data[i])                    
                    ok++;
            }

            if (ok == CodeLength)
                num_ok++;            
        }

        fprintf(stderr, "test runs......: %d\n",  num_runs);
        fprintf(stderr, "test runs OK...: %d\n",  num_ok);
        if (num_runs == num_ok)
            fprintf(stderr, "test runs OK...: PASS\n");
        else
            fprintf(stderr, "test runs OK...: FAIL\n");
    }
    else {
        FILE *fin, *fout;
        int   sdinput, readhalfframe, nread, offset, iter;

        /* File I/O mode ------------------------------------------------*/

        if (strcmp(argv[1], "-")  == 0) fin = stdin;
        else if ( (fin = fopen(argv[1],"rb")) == NULL ) {
            fprintf(stderr, "Error opening input SD file: %s: %s.\n",
                    argv[1], strerror(errno));
            exit(1);
        }
        
        if (strcmp(argv[2], "-") == 0) fout = stdout;
        else if ( (fout = fopen(argv[2],"wb")) == NULL ) {
            fprintf(stderr, "Error opening output bit file: %s: %s.\n",
                    argv[2], strerror(errno));
            exit(1);
        }

        sdinput = 0;
        readhalfframe = 0;
        if (opt_exists(argv, argc, "--sd")) {
            sdinput = 1;
        }
        if (opt_exists(argv, argc, "--half")) {
            readhalfframe = 1;
        }

        double *input_double = calloc(CodeLength, sizeof(double));

        nread = CodeLength;
        offset = 0;
        if (readhalfframe) {
            nread = CodeLength/2;
            offset = CodeLength/2;
            for(i=0; i<offset; i++) {
                input_double[i] = 0.0;
            }
        }

        while(fread(&input_double[offset], sizeof(double), nread, fin) == nread) {
            if (sdinput) {
                sd_to_llr(input_double, input_double, CodeLength);
            }

            iter = run_ldpc_decoder(&ldpc, out_char, input_double);
            fprintf(stderr, "%4d ", iter);

            // output data bits if decoder converged

            if (iter != MAX_ITER) {
              fwrite(out_char, sizeof(char), NUMBERROWSHCOLS, fout);
            }

            for(i=0; i<offset; i++) {
                input_double[i] = input_double[i+offset];
            }
        }

        free(input_double);
    }

    return 0;
}


