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

#include "ldpc_code.h"  

void run_ldpc_decoder(int DecodedBits[], int ParityCheckCount[], double input[]);

int main(int argc, char *argv[])
{    
    int         CodeLength, NumberParityBits, max_iter;
    int         i, j, r, num_ok, num_runs;
    char        out_char[CODELENGTH];

    /* derive some parameters */

    max_iter   = MAX_ITER;
    CodeLength = CODELENGTH;                    /* length of entire codeword */
    NumberParityBits = NUMBERPARITYBITS;
	
    if (argc < 2) {
        fprintf(stderr, "usage: %s --test\n", argv[0]);
        fprintf(stderr, "  Run internal self test and print code parameters.\n\n");
        fprintf(stderr, "usage: %s InOneSymbolPerDouble OutOneBitPerByte [--sdinput]\n", argv[0]);
        fprintf(stderr, "  InOneSymbolPerDouble is a file of double LLRs.  If the\n");
        fprintf(stderr, "  --sd flag is used the input file can be Soft Decision\n");
        fprintf(stderr, "  symbols, and LLRs will be calculated internally. Use -\n");
        fprintf(stderr, "  for the file names to use stdin/stdout.\n");
        exit(0);
    }

    int *DecodedBits = calloc( max_iter*CodeLength, sizeof( int ) );
    int *ParityCheckCount = calloc( max_iter, sizeof(int) );

    if (!strcmp(argv[1],"--test")) {

        /* test mode --------------------------------------------------------*/

        fprintf(stderr, "Starting test using pre-compiled test data .....\n");
        fprintf(stderr, "Codeword length: %d\n",  CodeLength);
        fprintf(stderr, "Parity Bits....: %d\n",  NumberParityBits);

        num_runs = 100; num_ok = 0;

        for(r=0; r<num_runs; r++) {

            run_ldpc_decoder(DecodedBits, ParityCheckCount, input);

            /* Check output by comparing every output iteration */

            int ok = 0;
            for (i=0;i<max_iter*CodeLength;i++) {
                if (output[i] == DecodedBits[i])                    
                            ok++;
            }
            if (ok == max_iter*CodeLength)
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
        int   sdinput;

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
        printf("argc: %d\n", argc);
        if (argc == 4)
            if (strcmp(argv[3], "--sdinput") == 0)
                sdinput = 1;

        double *input_double = calloc(CodeLength, sizeof(double));
        double sum, mean, sign, sumsq, estvar, estEsN0, x;

        while(fread(input_double, sizeof(double), CodeLength, fin) == CodeLength) {
            if (sdinput) {
                /* convert SD samples to LLRs -------------------------------*/

                sum = 0.0;
                for(i=0; i<CodeLength; i++)
                    sum += fabs(input_double[i]);
                mean = sum/CodeLength;
                
                /* scale by mean to map onto +/- 1 symbol position */

                for(i=0; i<CodeLength; i++) {
                    input_double[i] /= mean;
                }

                /* find variance from +/-1 symbol position */

                sum = sumsq = 0.0; 
                for(i=0; i<CodeLength; i++) {
                    sign = (input_double[i] > 0.0) - (input_double[i] < 0.0);
                    x = (input_double[i] - sign);
                    sum += x;
                    sumsq += x*x;
                }
                mean = sum/CodeLength;
                estvar = sumsq/CodeLength - mean*mean;

                estEsN0 = 1.0/(2.0 * estvar + 1E-3); 
                for(i=0; i<CodeLength; i++)
                    input_double[i] = 4.0 * estEsN0 * input_double[i];              
            }

            run_ldpc_decoder(DecodedBits, ParityCheckCount, input_double);
            
            /* extract output bits from ouput iteration that solved all parity equations, or failing that
               the last iteration. */

            int converged = 0;
            int iter = 0;
            for (i=0;i<max_iter;i++) {
                if (converged == 0)
                    iter++;
                if ((ParityCheckCount[i] == NumberParityBits)) {
                    for (j=0; j<CodeLength; j++) {
                        out_char[j] = DecodedBits[i+j*max_iter];
                    }
                    converged = 1;
                }               
            }
            if (converged == 0) {
                for (j=0; j<CodeLength; j++) {
                    out_char[j] = DecodedBits[max_iter-1+j*max_iter];
                }
            }
            //printf("%4d ", iter);
            fwrite(out_char, sizeof(char), CodeLength, fout);
        }

        free(input_double);
    }

    /* Clean up memory */

    free(ParityCheckCount);
    free(DecodedBits);

    return 0;
}


void run_ldpc_decoder(int DecodedBits[], int ParityCheckCount[], double input[]) {
    int		max_iter, dec_type;
    float       q_scale_factor, r_scale_factor;
    int		max_row_weight, max_col_weight;
    int         CodeLength, NumberParityBits, NumberRowsHcols, shift, H1;
    int         i;
    struct c_node *c_nodes;
    struct v_node *v_nodes;
    
    /* default values */

    max_iter  = MAX_ITER;
    dec_type  = 0;
    q_scale_factor = 1;
    r_scale_factor = 1;

    /* derive some parameters */

    CodeLength = CODELENGTH;                    /* length of entire codeword */
    NumberParityBits = NUMBERPARITYBITS;
    NumberRowsHcols = NUMBERROWSHCOLS;

    shift = (NumberParityBits + NumberRowsHcols) - CodeLength;
    if (NumberRowsHcols == CodeLength) {
        H1=0;
        shift=0;
    } else {
        H1=1;
    }
	
    max_row_weight = MAX_ROW_WEIGHT;
    max_col_weight = MAX_COL_WEIGHT;
	
    c_nodes = calloc( NumberParityBits, sizeof( struct c_node ) );
    v_nodes = calloc( CodeLength, sizeof( struct v_node));
	
    /* initialize c-node and v-node structures */

    c_nodes = calloc( NumberParityBits, sizeof( struct c_node ) );
    v_nodes = calloc( CodeLength, sizeof( struct v_node));
	
    init_c_v_nodes(c_nodes, shift, NumberParityBits, max_row_weight, H_rows, H1, CodeLength, 
                   v_nodes, NumberRowsHcols, H_cols, max_col_weight, dec_type, input);

    int DataLength = CodeLength - NumberParityBits;
    int *data_int = calloc( DataLength, sizeof(int) );
	
    /* need to clear these on each call */

    for(i=0; i<max_iter; i++)
        ParityCheckCount[i] = 0;
     for(i=0; i<max_iter*CodeLength; i++)
         DecodedBits[i] = 0;

    /* Call function to do the actual decoding */

    if ( dec_type == 1) {
        MinSum( ParityCheckCount, DecodedBits, c_nodes, v_nodes, CodeLength, 
                NumberParityBits, max_iter, r_scale_factor, q_scale_factor, data_int );
    } else if ( dec_type == 2) {
        fprintf(stderr, "dec_type = 2 not currently supported");
        /* ApproximateMinStar( BitErrors, DecodedBits, c_nodes, v_nodes, 
           CodeLength, NumberParityBits, max_iter, r_scale_factor, q_scale_factor );*/
    } else {
        SumProduct( ParityCheckCount, DecodedBits, c_nodes, v_nodes, CodeLength, 
                    NumberParityBits, max_iter, r_scale_factor, q_scale_factor, data_int ); 
    }

    /* Clean up memory */

    free( data_int );

    /*  Cleaning c-node elements */

    for (i=0;i<NumberParityBits;i++) {
        free( c_nodes[i].index );
        free( c_nodes[i].message );
        free( c_nodes[i].socket );
    }
	
    /* printf( "Cleaning c-nodes \n" ); */
    free( c_nodes );
	
    /* printf( "Cleaning v-node elements\n" ); */
    for (i=0;i<CodeLength;i++) {
        free( v_nodes[i].index);
        free( v_nodes[i].sign );
        free( v_nodes[i].message );
        free( v_nodes[i].socket );
    }
	
    /* printf( "Cleaning v-nodes \n" ); */
    free( v_nodes );
}
