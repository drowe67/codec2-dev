/* 
  FILE...: ldpc_dec.c
  AUTHOR.: Matthew C. Valenti, Rohit Iyer Seshadri, David Rowe
  CREATED: Sep 2016

  Command line C LDPC decoder derived from MpDecode.c in the CML
  library.  Allows us to run the same decoder in Octave and C.
*/

#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>

#include "mpdecode_core.h"
#include "ofdm_internal.h"

/* Machine generated consts, H_rows, H_cols, test input/output data to
   change LDPC code regenerate this file. */

/* TODO: Better octave/C support for multuple codes */

#include "H2064_516_sparse.h"  
#include "HRA_112_112.h"  
#include "HRAb_396_504.h"
#include "H_256_768_22.h"

int opt_exists(char *argv[], int argc, char opt[]) {
    int i;
    for (i=0; i<argc; i++) {
        if (strcmp(argv[i], opt) == 0) {
            return i;
        }
    }
    return 0;
}


int main(int argc, char *argv[])
{    
    int         CodeLength, NumberParityBits;
    int         i, r, num_ok, num_runs, codename, parityCheckCount, mute, testframes;
    int         data_bits_per_frame;
    char        *adetected_data;
    struct LDPC ldpc;
    float      *ainput;
    int         iter, total_iters;
    int         Tbits, Terrs, Tbits_raw, Terrs_raw, Tpackets, Tpacketerrs;

    int unused_data_bits = 84;

    if (argc < 2) {
        fprintf(stderr, "\n");
        fprintf(stderr, "usage: %s --test [--code CodeName]\n\n", argv[0]);
        fprintf(stderr, "  Run internal self test and print code parameters.\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "usage: %s --listcodes\n\n", argv[0]);
        fprintf(stderr, "  List supported codes (more can be added via using Octave ldpc scripts)\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "usage: %s InOneSymbolPerFloat OutOneBitPerByte [--sd] [--half] [--code CodeName] [--testframes]", argv[0]);
        fprintf(stderr, " [--unused numUnusedDataBits]\n\n");
        fprintf(stderr, "   InOneSymbolPerFloat     Input file of float LLRs, use - for the \n");        
        fprintf(stderr, "                           file names to use stdin/stdout\n");
        fprintf(stderr, "   --code                  Use LDPC code CodeName\n");
        fprintf(stderr, "   --listcodes             List available LDPC codes\n");
        fprintf(stderr, "   --sd                    Treat input file samples as BPSK Soft Decision\n");
        fprintf(stderr, "                           demod outputs rather than LLRs\n");
        fprintf(stderr, "   --mute                  Only output frames with < 10%% parity check fails\n");
        fprintf(stderr, "   --testframes            built in test frame modem, requires --testframes at encoder\n");
        fprintf(stderr, "    --unused               number of unused data bits, which are set to 1's at enc and dec\n");
        fprintf(stderr, "\n");

        fprintf(stderr, "Example in testframe mode:\n\n");
        fprintf(stderr, " $ ./ldpc_enc /dev/zero - --sd --code HRA_112_112 --testframes 10 |\n");
        fprintf(stderr, "   ./ldpc_dec - /dev/null --code HRA_112_112 --sd --testframes\n");
        exit(0);
    }


    if ((codename = opt_exists(argv, argc, "--listcodes")) != 0) {
        fprintf(stderr,"\n");
        fprintf(stderr,"H2064_516_sparse\n");
        fprintf(stderr,"HRA_112_112\n");
        fprintf(stderr,"HRAb_396_504\n");
        fprintf(stderr,"H_256_768\n");
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
            ldpc.H_rows = (uint16_t *)HRA_112_112_H_rows;
            ldpc.H_cols = (uint16_t *)HRA_112_112_H_cols;
            ainput = (float *)HRA_112_112_input;
            adetected_data = (char *)HRA_112_112_detected_data;
        }
        else if (strcmp(argv[codename+1], "HRAb_396_504") == 0) {
            fprintf(stderr, "code: %s\n", argv[codename+1]);
            ldpc.max_iter = HRAb_396_504_MAX_ITER;
            ldpc.dec_type = 0;
            ldpc.q_scale_factor = 1;
            ldpc.r_scale_factor = 1;
            ldpc.CodeLength = HRAb_396_504_CODELENGTH;
            ldpc.NumberParityBits = HRAb_396_504_NUMBERPARITYBITS;
            ldpc.NumberRowsHcols = HRAb_396_504_NUMBERROWSHCOLS;
            ldpc.max_row_weight = HRAb_396_504_MAX_ROW_WEIGHT;
            ldpc.max_col_weight = HRAb_396_504_MAX_COL_WEIGHT;
            ldpc.H_rows = (uint16_t *)HRAb_396_504_H_rows;
            ldpc.H_cols = (uint16_t *)HRAb_396_504_H_cols;
            ainput = (float *)HRAb_396_504_input;
            adetected_data = (char *)HRAb_396_504_detected_data;
        }
        else if (strcmp(argv[codename+1], "H_256_768") == 0) {
            fprintf(stderr, "code: %s\n", argv[codename+1]);
            ldpc.max_iter = H_256_768_22_MAX_ITER;
            ldpc.dec_type = 0;
            ldpc.q_scale_factor = 1;
            ldpc.r_scale_factor = 1;
            ldpc.CodeLength = H_256_768_22_CODELENGTH;
            ldpc.NumberParityBits = H_256_768_22_NUMBERPARITYBITS;
            ldpc.NumberRowsHcols = H_256_768_22_NUMBERROWSHCOLS;
            ldpc.max_row_weight = H_256_768_22_MAX_ROW_WEIGHT;
            ldpc.max_col_weight = H_256_768_22_MAX_COL_WEIGHT;
            ldpc.H_rows = (uint16_t *)H_256_768_22_H_rows;
            ldpc.H_cols = (uint16_t *)H_256_768_22_H_cols;
            ainput = (float *)H_256_768_22_input;
            adetected_data = (char *)H_256_768_22_detected_data;
        }
        else {
            fprintf(stderr, "Unknown code: %s, defaulting to H2064_516_sparse\n", argv[codename+1]);
        }
    }

    CodeLength = ldpc.CodeLength;                    /* length of entire codeword */
    NumberParityBits = ldpc.NumberParityBits;
    data_bits_per_frame = ldpc.NumberRowsHcols;
    unsigned char ibits[data_bits_per_frame];
    unsigned char pbits[NumberParityBits];
    uint8_t out_char[CodeLength];

    testframes = 0;
    total_iters = 0;
    Tbits = Terrs = Tbits_raw = Terrs_raw = Tpackets = Tpacketerrs = 0;
    
    if (!strcmp(argv[1],"--test")) {

        /* test mode --------------------------------------------------------*/

        fprintf(stderr, "Starting test using pre-compiled test data .....\n");
        fprintf(stderr, "Codeword length: %d\n",  CodeLength);
        fprintf(stderr, "Parity Bits....: %d\n",  NumberParityBits);

        num_runs = 100; num_ok = 0;

        for(r=0; r<num_runs; r++) {

            iter = run_ldpc_decoder(&ldpc, out_char, ainput, &parityCheckCount);
            //fprintf(stderr, "iter: %d\n", iter);
	    total_iters += iter;

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
        int   sdinput, nread, offset;

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

        sdinput = 0; mute = 0; testframes = 0;
        if (opt_exists(argv, argc, "--sd")) {
            sdinput = 1;
        }
        if (opt_exists(argv, argc, "--mute")) {
            mute = 1;
        }
        unused_data_bits = 0; int arg;
        if ((arg = opt_exists(argv, argc, "--unused"))) {
            unused_data_bits = atoi(argv[arg+1]);
        }
        if (opt_exists(argv, argc, "--testframes")) {
            testframes = 1;
            uint16_t r[data_bits_per_frame];
            ofdm_rand(r, data_bits_per_frame);

            for(i=0; i<data_bits_per_frame-unused_data_bits; i++) {
                ibits[i] = r[i] > 16384;
            }
            for(i=data_bits_per_frame-unused_data_bits; i<data_bits_per_frame; i++) {
                ibits[i] = 1;
            }
            encode(&ldpc, ibits, pbits);  
        }

        float  *input_float  = calloc(CodeLength, sizeof(float));

        nread = CodeLength - unused_data_bits;
        fprintf(stderr, "CodeLength: %d offset: %d\n", CodeLength, offset);

        while(fread(input_float, sizeof(float), nread, fin) == nread) {
            if (testframes) {
                /* BPSK SD and bit LLRs get mapped roughly the same way so this just happens to work for both */
                char in_char;
                for (i=0; i<data_bits_per_frame-unused_data_bits; i++) {
                    in_char = input_float[i] < 0;
                    if (in_char != ibits[i]) {
                        Terrs_raw++;
                    }
                    Tbits_raw++;
                }
                for (i=0; i<NumberParityBits; i++) {
                    in_char = input_float[i+data_bits_per_frame-unused_data_bits] < 0;
                    if (in_char != pbits[i]) {
                        Terrs_raw++;
                    }
                    Tbits_raw++;
                }
            }

            if (sdinput) {
                /* map BPSK SDs to bit LLRs */
                float llr[CodeLength-unused_data_bits];
                sd_to_llr(llr, input_float, CodeLength-unused_data_bits);

                /* insert unused data LLRs */

                float llr_tmp[CodeLength];
                for(i=0; i<data_bits_per_frame-unused_data_bits; i++)
                    llr_tmp[i] = llr[i];  // rx data bits
                for(i=data_bits_per_frame-unused_data_bits; i<data_bits_per_frame; i++)
                    llr_tmp[i] = -10.0;           // known data bits high likelhood
                for(i=data_bits_per_frame; i<CodeLength; i++)
                    llr_tmp[i] = llr[i-unused_data_bits];  // rx parity bits
                memcpy(input_float, llr_tmp, sizeof(float)*CodeLength);                
            }

            iter = run_ldpc_decoder(&ldpc, out_char, input_float, &parityCheckCount);
            //fprintf(stderr, "iter: %d\n", iter);
	    total_iters += iter;
            
            if (mute) {

                // Output data bits only if decoder converged, or was
                // within 10% of all parity checks converging (10% est
                // BER).  Useful for real world operation as it can
                // resync and won't send crappy packets to the decoder
                
                float ber_est = (float)(ldpc.NumberParityBits - parityCheckCount)/ldpc.NumberParityBits;
                //fprintf(stderr, "iter: %4d parityCheckErrors: %4d ber: %3.2f\n", iter, ldpc.NumberParityBits - parityCheckCount, ber_est);
                if (ber_est < 0.1) {
                    fwrite(out_char, sizeof(char), ldpc.NumberRowsHcols, fout);
                }

            }

            fwrite(out_char, sizeof(char), data_bits_per_frame, fout);

            if (testframes) {
                int perr = 0;
                for (i=0; i<data_bits_per_frame; i++) {
                    if (out_char[i] != ibits[i]) {
                        Terrs++;
                        perr = 1;
                    }
                    Tbits++;
                }
                Tpackets++; if (perr) { Tpacketerrs++; fprintf(stderr,"x"); } else fprintf(stderr,".");
            }
            //fprintf(stderr, "Terrs_raw: %d  Tbits_raw: %d Terr: %d Tbits: %d\n", Terrs_raw, Tbits_raw, Terrs, Tbits);
        }

        free(input_float);
        if (fin  != NULL) fclose(fin);
        if (fout != NULL) fclose(fout);
    }

    fprintf(stderr, "total iters %d\n", total_iters);

    if (testframes) {
        fprintf(stderr, "Raw   Tbits: %6d Terr: %6d BER: %4.3f\n", Tbits_raw, Terrs_raw,
                (float)Terrs_raw/(Tbits_raw+1E-12));
        float coded_ber = (float)Terrs/(Tbits+1E-12);
        fprintf(stderr, "Coded Tbits: %6d Terr: %6d BER: %4.3f\n", Tbits, Terrs, coded_ber);
        fprintf(stderr, "      Tpkts: %6d Tper: %6d PER: %4.3f\n", Tpackets, Tpacketerrs, Tpacketerrs/(Tpackets+1E-12));

        /* set return code for Ctest */
        if (coded_ber < 0.01)
            return 0;
        else
            return 1;
    }
    
    return 0;
}


