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

/* TODO: Better octave/C support for multuple codes */

#include "H2064_516_sparse.h"  
#include "HRA_112_112.h"  

int opt_exists(char *argv[], int argc, char opt[]) {
    int i;
    for (i=0; i<argc; i++) {
        if (strcmp(argv[i], opt) == 0) {
            return i;
        }
    }
    return 0;
}

void extract_output(char out_char[], int DecodedBits[], int ParityCheckCount[], int max_iter, int CodeLength, int NumberParityBits);

int main(int argc, char *argv[])
{    
    int         CodeLength, NumberParityBits;
    int         i, r, num_ok, num_runs, codename, parityCheckCount, mute, state, next_state, frame;
    char        out_char[CODELENGTH], *adetected_data;
    struct LDPC ldpc;
    double     *ainput;
	
    if (argc < 2) {
        fprintf(stderr, "\n");
        fprintf(stderr, "usage: %s --test [--code CodeName]\n\n", argv[0]);
        fprintf(stderr, "  Run internal self test and print code parameters.\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "usage: %s --listcodes\n\n", argv[0]);
        fprintf(stderr, "  List supported codes (more can be added via using Octave ldpc scripts)\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "usage: %s InOneSymbolPerDouble OutOneBitPerByte [--sd] [--half] [--code CodeName]\n\n", argv[0]);
        fprintf(stderr, "   InOneSymbolPerDouble    Input file of double LLRs, use - for the \n");        
        fprintf(stderr, "                           file names to use stdin/stdout\n");
        fprintf(stderr, "   --code                  Use LDPC code CodeName\n");
        fprintf(stderr, "   --listcodes             List available LDPC codes\n");
        fprintf(stderr, "   --sd                    Treat input file samples as Soft Decision\n");
        fprintf(stderr, "                           demod outputs rather than LLRs\n");
        fprintf(stderr, "   --half                  Load framesize/2 input samples for each decode\n");
        fprintf(stderr, "                           attempt, only output decoded bits if decoder\n");
        fprintf(stderr, "                           converges.  Form of frame sync.\n");
        fprintf(stderr, "   --mute                  Only output frames with < 10%% parity check fails\n");
        fprintf(stderr, "\n");
        exit(0);
    }


    if ((codename = opt_exists(argv, argc, "--listcodes")) != 0) {
        fprintf(stderr,"\n");
        fprintf(stderr,"H2064_516_sparse\n");
        fprintf(stderr,"HRA_112_112\n");
        fprintf(stderr,"\n");
        exit(0);
    }

    /* default Wenet High Alitiude Balloon rate 0.8 code */

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
    ainput = input;
    adetected_data = detected_data;

    if ((codename = opt_exists(argv, argc, "--code")) != 0) {

        /* short rate 1/2 code for FreeDV HF digital voice */
        
        if (strcmp(argv[codename+1], "HRA_112_112") == 0) {
            fprintf(stderr, "code: %s\n", argv[codename+1]);
            ldpc.max_iter = HRA_112_112_MAX_ITER;
            ldpc.dec_type = 0;
            ldpc.q_scale_factor = 1;
            ldpc.r_scale_factor = 1;
            ldpc.CodeLength = HRA_112_112_CODELENGTH;
            ldpc.NumberParityBits = HRA_112_112_NUMBERPARITYBITS;
            ldpc.NumberRowsHcols = HRA_112_112_NUMBERROWSHCOLS;
            ldpc.max_row_weight = HRA_112_112_MAX_ROW_WEIGHT;
            ldpc.max_col_weight = HRA_112_112_MAX_COL_WEIGHT;
            ldpc.H_rows = HRA_112_112_H_rows;
            ldpc.H_cols = HRA_112_112_H_cols;
            ainput = HRA_112_112_input;
            adetected_data = HRA_112_112_detected_data;
        }
    }

    CodeLength = ldpc.CodeLength;                    /* length of entire codeword */
    NumberParityBits = ldpc.NumberParityBits;

    if (!strcmp(argv[1],"--test")) {

        /* test mode --------------------------------------------------------*/

        fprintf(stderr, "Starting test using pre-compiled test data .....\n");
        fprintf(stderr, "Codeword length: %d\n",  CodeLength);
        fprintf(stderr, "Parity Bits....: %d\n",  NumberParityBits);

        num_runs = 100; num_ok = 0;

        for(r=0; r<num_runs; r++) {

            run_ldpc_decoder(&ldpc, out_char, ainput, &parityCheckCount);

            int ok = 0;
            for (i=0; i<CodeLength; i++) {
                if (out_char[i] == adetected_data[i])                    
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
        mute = 0; state = 0; frame = 0;
        if (opt_exists(argv, argc, "--sd")) {
            sdinput = 1;
        }
        if (opt_exists(argv, argc, "--half")) {
            readhalfframe = 1;
        }
        if (opt_exists(argv, argc, "--mute")) {
            mute = 1;
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
        fprintf(stderr, "CodeLength: %d offset: %d\n", CodeLength, offset);

        while(fread(&input_double[offset], sizeof(double), nread, fin) == nread) {
            if (sdinput) {
                sd_to_llr(input_double, input_double, CodeLength);
            }

            iter = run_ldpc_decoder(&ldpc, out_char, input_double, &parityCheckCount);

            if (mute) {

                // Output data bits if decoder converged, or was
                // within 10% of all parity checks converging (10% est
                // BER).  useful for real world operation as it can
                // resync and won't send crappy packets to the decoder
                
                float ber_est = (float)(ldpc.NumberParityBits - parityCheckCount)/ldpc.NumberParityBits;
                //fprintf(stderr, "iter: %4d parityCheckErrors: %4d ber: %3.2f\n", iter, ldpc.NumberParityBits - parityCheckCount, ber_est);
                if (ber_est < 0.1) {
                    fwrite(out_char, sizeof(char), ldpc.NumberRowsHcols, fout);
                }

            } else {
                
                // Output all data packets, based on initial FEC sync
                // estimate.  Useful for testing with cohpsk_put_bits,
                // as it maintains sync with test bits state machine.
                
                next_state = state;
                switch(state) {
                case 0:
                    if (iter < MAX_ITER) {
                        /* OK we've found which frame to sync on */
                        next_state = 1;
                        frame = 0;
                    }
                    break;
                case 1:
                    frame++;
                    if ((frame % 2) == 0) {
                        /* write decoded packets every second input frame */
                        fwrite(out_char, sizeof(char), ldpc.NumberRowsHcols, fout);
                    }
                    break;
                }
                state = next_state;
                fprintf(stderr, "state: %d iter: %d\n", state, iter);
            }

            for(i=0; i<offset; i++) {
                input_double[i] = input_double[i+offset];
            }
        }

        free(input_double);
        if (fin  != NULL) fclose(fin);
        if (fout != NULL) fclose(fout);
    }

    return 0;
}


