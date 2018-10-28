/*
  FILE...: mpdecode_core.c
  AUTHOR.: Matthew C. Valenti, Rohit Iyer Seshadri, David Rowe
  CREATED: Sep 2016

  C-callable core functions moved from MpDecode.c, so they can be used for
  Octave and C programs.
*/

#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include "mpdecode_core.h"
#ifndef USE_ORIGINAL_PHI0
#include "phi0.h"
#endif

#ifdef __EMBEDDED__
#include "machdep.h"
#endif

#define QPSK_CONSTELLATION_SIZE 4
#define QPSK_BITS_PER_SYMBOL    2

#undef PRINT_PROGRESS
#undef PRINT_ALLOCS

/* QPSK constellation for symbol likelihood calculations */

static COMP S_matrix[] = {
    { 1.0f,  0.0f},
    { 0.0f,  1.0f},
    { 0.0f, -1.0f},
    {-1.0f,  0.0f}
};


void encode(struct LDPC *ldpc, unsigned char ibits[], unsigned char pbits[]) {
    unsigned int p, i, tmp, par, prev=0;
    int          ind;
    uint16_t     *H_rows = ldpc->H_rows;

    for (p=0; p<ldpc->NumberParityBits; p++) {
        par = 0;

        for (i=0; i<ldpc->max_row_weight; i++) {
            ind = H_rows[p + i*ldpc->NumberParityBits];
            par = par + ibits[ind-1];
        }

        tmp = par + prev;

        tmp &= 1;    // only retain the lsb
        prev = tmp;
        pbits[p] = tmp;
    }
}

#ifdef USE_ORIGINAL_PHI0
/* Phi function */
static float phi0(
                  float x )
{
  float z;

  if (x>10)
    return( 0 );
  else if (x< 9.08e-5 )
    return( 10 );
  else if (x > 9)
    return( 1.6881e-4 );
  /* return( 1.4970e-004 ); */
  else if (x > 8)
    return( 4.5887e-4 );
  /* return( 4.0694e-004 ); */
  else if (x > 7)
    return( 1.2473e-3 );
  /* return( 1.1062e-003 ); */
  else if (x > 6)
    return( 3.3906e-3 );
  /* return( 3.0069e-003 ); */
  else if (x > 5)
    return( 9.2168e-3 );
  /* return( 8.1736e-003 ); */
  else {
    z = (float) exp(x);
    return( (float) log( (z+1)/(z-1) ) );
  }
}
#endif


/* Values for linear approximation (DecoderType=5) */

#define AJIAN -0.24904163195436
#define TJIAN 2.50681740420944

/* The linear-log-MAP algorithm */

static float max_star0(
                       float delta1,
                       float delta2 )
{
    register float diff;

    diff = delta2 - delta1;

    if ( diff > TJIAN )
        return( delta2 );
    else if ( diff < -TJIAN )
        return( delta1 );
    else if ( diff > 0 )
        return( delta2 + AJIAN*(diff-TJIAN) );
    else
        return( delta1 - AJIAN*(diff+TJIAN) );
}

void init_c_v_nodes(struct c_node *c_nodes,
                    int     shift,
                    int     NumberParityBits,
                    int     max_row_weight,
                    uint16_t *H_rows,
                    int     H1,
                    int     CodeLength,
                    struct v_node *v_nodes,
                    int     NumberRowsHcols,
                    uint16_t *H_cols,
                    int     max_col_weight,
                    int     dec_type,
                    float  *input)
{
    int i, j, k, count, cnt, c_index, v_index;

    /* first determine the degree of each c-node */

    if (shift ==0){
        for (i=0;i<NumberParityBits;i++) {
            count = 0;
            for (j=0;j<max_row_weight;j++) {
                if ( H_rows[i+j*NumberParityBits] > 0 ) {
                    count++;
                }
            }
            c_nodes[i].degree = count;
            if (H1){
                if (i==0){
                    c_nodes[i].degree=count+1;
                }
                else{
                    c_nodes[i].degree=count+2;
                }
            }
        }
    }
    else{
        cnt=0;
        for (i=0;i<(NumberParityBits/shift);i++) {
            for (k=0;k<shift;k++){
                count = 0;
                for (j=0;j<max_row_weight;j++) {
                    if ( H_rows[cnt+j*NumberParityBits] > 0 ) {
                        count++;
                    }
                }
                c_nodes[cnt].degree = count;
                if ((i==0)||(i==(NumberParityBits/shift)-1)){
                    c_nodes[cnt].degree=count+1;
                }
                else{
                    c_nodes[cnt].degree=count+2;
                }
                cnt++;
            }
        }
    }

    if (H1){

        if (shift ==0){
            for (i=0;i<NumberParityBits;i++) {
                /* now that we know the size, we can dynamically allocate memory */
                c_nodes[i].index =  calloc( c_nodes[i].degree, sizeof( int ) );
                c_nodes[i].message =calloc( c_nodes[i].degree, sizeof( float ) );
                c_nodes[i].socket = calloc( c_nodes[i].degree, sizeof( int ) );

                for (j=0;j<c_nodes[i].degree-2;j++) {
                    c_nodes[i].index[j] = (int) (H_rows[i+j*NumberParityBits] - 1);
                }
                j=c_nodes[i].degree-2;

                if (i==0){
                    c_nodes[i].index[j] = (int) (H_rows[i+j*NumberParityBits] - 1);
                }
                else {
                    c_nodes[i].index[j] = (CodeLength-NumberParityBits)+i-1;
                }

                j=c_nodes[i].degree-1;
                c_nodes[i].index[j] = (CodeLength-NumberParityBits)+i;

            }
        }
        if (shift >0){
            cnt=0;
            for (i=0;i<(NumberParityBits/shift);i++){

                for (k =0;k<shift;k++){
                    c_nodes[cnt].index =  calloc( c_nodes[cnt].degree, sizeof( int ) );
                    c_nodes[cnt].message =calloc( c_nodes[cnt].degree, sizeof( float ) );
                    c_nodes[cnt].socket = calloc( c_nodes[cnt].degree, sizeof( int ) );

                    for (j=0;j<c_nodes[cnt].degree-2;j++) {
                        c_nodes[cnt].index[j] = (int) (H_rows[cnt+j*NumberParityBits] - 1);
                    }
                    j=c_nodes[cnt].degree-2;
                    if ((i ==0)||(i==(NumberParityBits/shift-1))){
                        c_nodes[cnt].index[j] = (int) (H_rows[cnt+j*NumberParityBits] - 1);
                    }
                    else{
                        c_nodes[cnt].index[j] = (CodeLength-NumberParityBits)+k+shift*(i);
                    }
                    j=c_nodes[cnt].degree-1;
                    c_nodes[cnt].index[j] = (CodeLength-NumberParityBits)+k+shift*(i+1);
                    if (i== (NumberParityBits/shift-1))
                        {
                            c_nodes[cnt].index[j] = (CodeLength-NumberParityBits)+k+shift*(i);
                        }
                    cnt++;
                }
            }
        }

    } else {
        for (i=0;i<NumberParityBits;i++) {
            /* now that we know the size, we can dynamically allocate memory */
            c_nodes[i].index =  calloc( c_nodes[i].degree, sizeof( int ) );
            c_nodes[i].message =calloc( c_nodes[i].degree, sizeof( float ) );
            c_nodes[i].socket = calloc( c_nodes[i].degree, sizeof( int ) );
            for (j=0;j<c_nodes[i].degree;j++){
                c_nodes[i].index[j] = (int) (H_rows[i+j*NumberParityBits] - 1);
            }
        }
    }


    /* determine degree of each v-node */

    for(i=0;i<(CodeLength-NumberParityBits+shift);i++){
        count=0;
        for (j=0;j<max_col_weight;j++) {
            if ( H_cols[i+j*NumberRowsHcols] > 0 ) {
                count++;
            }
        }
        v_nodes[i].degree = count;
    }

    for(i=CodeLength-NumberParityBits+shift;i<CodeLength;i++){
        count=0;
        if (H1){
            if(i!=CodeLength-1){
                v_nodes[i].degree=2;
            }  else{
                v_nodes[i].degree=1;
            }

        } else{
            for (j=0;j<max_col_weight;j++) {
                if ( H_cols[i+j*NumberRowsHcols] > 0 ) {
                    count++;
                }
            }
            v_nodes[i].degree = count;
        }
    }

    if (shift>0){
        v_nodes[CodeLength-1].degree =v_nodes[CodeLength-1].degree+1;
    }


    /* set up v_nodes */

    for (i=0;i<CodeLength;i++) {
        /* allocate memory according to the degree of the v-node */
        v_nodes[i].index = calloc( v_nodes[i].degree, sizeof( int ) );
        v_nodes[i].message = calloc( v_nodes[i].degree, sizeof( float ) );
        v_nodes[i].sign = calloc( v_nodes[i].degree, sizeof( int ) );
        v_nodes[i].socket = calloc( v_nodes[i].degree, sizeof( int ) );

        /* index tells which c-nodes this v-node is connected to */
        v_nodes[i].initial_value = input[i];
        count=0;

        for (j=0;j<v_nodes[i].degree;j++) {
            if ((H1)&& (i>=CodeLength-NumberParityBits+shift)){
                v_nodes[i].index[j]=i-(CodeLength-NumberParityBits+shift)+count;
                if (shift ==0){
                    count=count+1;
                }
                else{
                    count=count+shift;
                }
            } else  {
                v_nodes[i].index[j] = (int) (H_cols[i+j*NumberRowsHcols] - 1);
            }

            /* search the connected c-node for the proper message value */
            for (c_index=0;c_index<c_nodes[ v_nodes[i].index[j] ].degree;c_index++)
                if ( c_nodes[ v_nodes[i].index[j] ].index[c_index] == i ) {
                    v_nodes[i].socket[j] = c_index;
                    break;
                }
            /* initialize v-node with received LLR */
            if ( dec_type == 1)
                v_nodes[i].message[j] = fabs(input[i]);
            else
                v_nodes[i].message[j] = phi0( fabs(input[i]) );

            if (input[i] < 0)
                v_nodes[i].sign[j] = 1;
        }

    }



    /* now finish setting up the c_nodes */
    for (i=0;i<NumberParityBits;i++) {
        /* index tells which v-nodes this c-node is connected to */
        for (j=0;j<c_nodes[i].degree;j++) {
            /* search the connected v-node for the proper message value */
            for (v_index=0;v_index<v_nodes[ c_nodes[i].index[j] ].degree;v_index++)
                if (v_nodes[ c_nodes[i].index[j] ].index[v_index] == i ) {
                    c_nodes[i].socket[j] = v_index;
                    break;
                }
        }
    }

}


/* function for doing the MP decoding */
// Returns the iteration count
int SumProduct( int       *parityCheckCount,
                char     DecodedBits[],
                struct c_node c_nodes[],
                struct v_node v_nodes[],
                int       CodeLength,
                int       NumberParityBits,
                int       max_iter,
                float    r_scale_factor,
                float    q_scale_factor,
                int      data[] )
{
  int result;
  int bitErrors;
  int i,j, iter;
  float phi_sum;
  int sign;
  float temp_sum;
  float Qi;
  int   ssum;

  #ifdef PRINT_PROGRESS
  fprintf(stderr, "SumProduct\n");
  #endif

//#ifdef __EMBEDDED__
//PROFILE_VAR(ldpc_SP_iter, ldpc_SP_upr, ldpc_SP_upq, ldpc_SP_misc);
//#endif

  result = max_iter;
  for (iter=0;iter<max_iter;iter++) {
    #ifdef PRINT_PROGRESS
    fprintf(stderr, "  iter %d\n", iter);
    #endif

    for(i=0; i<CodeLength; i++) DecodedBits[i] = 0; // Clear each pass!
    bitErrors = 0;

    /* update r */
    ssum = 0;
    for (j=0;j<NumberParityBits;j++) {
      sign = v_nodes[ c_nodes[j].index[0] ].sign[ c_nodes[j].socket[0] ];
      phi_sum = v_nodes[ c_nodes[j].index[0] ].message[ c_nodes[j].socket[0] ];

      for (i=1;i<c_nodes[j].degree;i++) {
        phi_sum += v_nodes[ c_nodes[j].index[i] ].message[ c_nodes[j].socket[i] ];
        sign ^= v_nodes[ c_nodes[j].index[i] ].sign[ c_nodes[j].socket[i] ];
      }

      if (sign==0) ssum++;

      for (i=0;i<c_nodes[j].degree;i++) {
        if ( sign^v_nodes[ c_nodes[j].index[i] ].sign[ c_nodes[j].socket[i] ] ) {
          c_nodes[j].message[i] = -phi0( phi_sum - v_nodes[ c_nodes[j].index[i] ].message[ c_nodes[j].socket[i] ] )*r_scale_factor;
        } else
          c_nodes[j].message[i] = phi0( phi_sum - v_nodes[ c_nodes[j].index[i] ].message[ c_nodes[j].socket[i] ] )*r_scale_factor;
      }
    }

    /* update q */
    for (i=0;i<CodeLength;i++) {

      /* first compute the LLR */
      Qi = v_nodes[i].initial_value;
      for (j=0;j<v_nodes[i].degree;j++) {
        Qi += c_nodes[ v_nodes[i].index[j] ].message[ v_nodes[i].socket[j] ];
      }

      /* make hard decision */
      if (Qi < 0) {
            DecodedBits[i] = 1;
      }

      /* now subtract to get the extrinsic information */
      for (j=0;j<v_nodes[i].degree;j++) {
        temp_sum = Qi - c_nodes[ v_nodes[i].index[j] ].message[ v_nodes[i].socket[j] ];

        v_nodes[i].message[j] = phi0( fabs( temp_sum ) )*q_scale_factor;
        if (temp_sum > 0)
          v_nodes[i].sign[j] = 0;
        else
          v_nodes[i].sign[j] = 1;
      }
    }

    /* count data bit errors, assuming that it is systematic */
    for (i=0;i<CodeLength-NumberParityBits;i++)
      if ( DecodedBits[i] != data[i] )
            bitErrors++;

    #ifdef PRINT_PROGRESS
    fprintf(stderr, "    bitErrors is %d \n", bitErrors);
    #endif

    /* Halt if zero errors */
    if (bitErrors == 0) {
      #ifdef PRINT_PROGRESS
      fprintf(stderr, "    SumProducts: 0 errors\n");
      #endif
      result = iter + 1;
      break;
    }

    // count the number of PC satisfied and exit if all OK
    #ifdef PRINT_PROGRESS
    fprintf(stderr, "    ssum is %d \n", ssum);
    #endif
    *parityCheckCount = ssum;
    if (ssum==NumberParityBits)  {
      #ifdef PRINT_PROGRESS
      fprintf(stderr, "    SumProducts: ssum == NumParityBits\n");
      #endif
      result = iter + 1;
      break;
    }


  }

#ifdef PRINT_PROGRESS
fprintf(stderr, "SumProducts %d iterations\n", result);
#endif
return(result);
}


/* Convenience function to call LDPC decoder from C programs */

int run_ldpc_decoder(struct LDPC *ldpc, char out_char[], float input[], int *parityCheckCount) {
    int         max_iter, dec_type;
    float       q_scale_factor, r_scale_factor;
    int         max_row_weight, max_col_weight;
    int         CodeLength, NumberParityBits, NumberRowsHcols, shift, H1;
    int         i;
    struct c_node *c_nodes;
    struct v_node *v_nodes;

    /* default values */

    max_iter  = ldpc->max_iter;
    dec_type  = ldpc->dec_type;
    q_scale_factor = ldpc->q_scale_factor;
    r_scale_factor = ldpc->r_scale_factor;

    CodeLength = ldpc->CodeLength;                    /* length of entire codeword */
    NumberParityBits = ldpc->NumberParityBits;
    NumberRowsHcols = ldpc->NumberRowsHcols;

    char *DecodedBits = calloc( CodeLength, sizeof( char ) );

    /* derive some parameters */

    shift = (NumberParityBits + NumberRowsHcols) - CodeLength;
    if (NumberRowsHcols == CodeLength) {
        H1=0;
        shift=0;
    } else {
        H1=1;
    }

    max_row_weight = ldpc->max_row_weight;
    max_col_weight = ldpc->max_col_weight;
    /*
    c_nodes = calloc( NumberParityBits, sizeof( struct c_node ) );
    v_nodes = calloc( CodeLength, sizeof( struct v_node));
    */
    /* initialize c-node and v-node structures */

    c_nodes = calloc( NumberParityBits, sizeof( struct c_node ) );
    v_nodes = calloc( CodeLength, sizeof( struct v_node));

    init_c_v_nodes(c_nodes, shift, NumberParityBits, max_row_weight, ldpc->H_rows, H1, CodeLength,
                   v_nodes, NumberRowsHcols, ldpc->H_cols, max_col_weight, dec_type, input);

    int DataLength = CodeLength - NumberParityBits;
    int *data_int = calloc( DataLength, sizeof(int) );

    /* need to clear these on each call */

    for(i=0; i<CodeLength; i++) DecodedBits[i] = 0;

    /* Call function to do the actual decoding */
    int iter = SumProduct( parityCheckCount, DecodedBits, c_nodes, v_nodes, 
                           CodeLength, NumberParityBits, max_iter, 
                           r_scale_factor, q_scale_factor, data_int );

    for (i=0; i<CodeLength; i++) out_char[i] = DecodedBits[i];

    /* Clean up memory */

    free(DecodedBits);
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

    return iter;
}


void sd_to_llr(float llr[], double sd[], int n) {
    double sum, mean, sign, sumsq, estvar, estEsN0, x;
    int i;

    /* convert SD samples to LLRs -------------------------------*/

    sum = 0.0;
    for(i=0; i<n; i++)
        sum += fabs(sd[i]);
    mean = sum/n;

    /* find variance from +/-1 symbol position */

    sum = sumsq = 0.0;
    for(i=0; i<n; i++) {
        sign = (sd[i] > 0.0L) - (sd[i] < 0.0L);
        x = (sd[i]/mean - sign);
        sum += x;
        sumsq += x*x;
    }
    estvar = (n * sumsq - sum * sum) / (n * (n - 1));
    //fprintf(stderr, "mean: %f var: %f\n", mean, estvar);

    estEsN0 = 1.0/(2.0L * estvar + 1E-3);
    for(i=0; i<n; i++)
        llr[i] = 4.0L * estEsN0 * sd[i];
}


/*
   Determine symbol likelihood from received QPSK symbols.

   Notes:

   1) We assume fading[] is real, it is also possible to compute
      with complex fading, see CML library Demod2D.c source code.
   2) Using floats instead of doubles, for stm32.
      Testing shows good BERs with floats.
*/

void Demod2D(float   symbol_likelihood[],       /* output, M*number_symbols              */
             COMP    r[],                       /* received QPSK symbols, number_symbols */
             COMP    S_matrix[],                /* constellation of size M               */
             float   EsNo,
             float   fading[],                  /* real fading values, number_symbols    */
             float   mean_amp,
             int     number_symbols)
{
    int     M=QPSK_CONSTELLATION_SIZE;
    int     i,j;
    float  tempsr, tempsi, Er, Ei;

    /* determine output */

    for (i=0;i<number_symbols;i++) {                /* go through each received symbol */
        for (j=0;j<M;j++) {                         /* each postulated symbol          */
            tempsr = fading[i]*S_matrix[j].real/mean_amp;
            tempsi = fading[i]*S_matrix[j].imag/mean_amp;
            Er = r[i].real/mean_amp - tempsr;
            Ei = r[i].imag/mean_amp - tempsi;
            symbol_likelihood[i*M+j] = -EsNo*(Er*Er+Ei*Ei);
            //printf("symbol_likelihood[%d][%d] = %f\n", i,j,symbol_likelihood[i*M+j]);
        }
        //exit(0);
    }

}


void Somap(float  bit_likelihood[],      /* number_bits, bps*number_symbols */
           float  symbol_likelihood[],   /* M*number_symbols                */
           int     number_symbols)
{
    int    M=QPSK_CONSTELLATION_SIZE, bps = QPSK_BITS_PER_SYMBOL;
    int    n,i,j,k,mask;
    float num[bps], den[bps];
    float metric;

    for (n=0; n<number_symbols; n++) { /* loop over symbols */
        for (k=0;k<bps;k++) {
            /* initialize */
            num[k] = -1000000;
            den[k] = -1000000;
        }

        for (i=0;i<M;i++) {
            metric =  symbol_likelihood[n*M+i]; /* channel metric for this symbol */

            mask = 1 << (bps - 1);
            for (j=0;j<bps;j++) {
                mask = mask >> 1;
            }
            mask = 1 << (bps - 1);

            for (k=0;k<bps;k++) {       /* loop over bits */
                if (mask&i) {
                    /* this bit is a one */
                    num[k] = max_star0( num[k], metric );
                } else {
                    /* this bit is a zero */
                    den[k] = max_star0( den[k], metric );
                }
                mask = mask >> 1;
            }
        }
        for (k=0;k<bps;k++) {
            bit_likelihood[bps*n+k] = num[k] - den[k];
        }
    }
}


void symbols_to_llrs(float llr[], COMP rx_qpsk_symbols[], float rx_amps[], float EsNo, float mean_amp, int nsyms) {
    int i;

    float symbol_likelihood[nsyms*QPSK_CONSTELLATION_SIZE];
    float bit_likelihood[nsyms*QPSK_BITS_PER_SYMBOL];

    Demod2D(symbol_likelihood, rx_qpsk_symbols, S_matrix, EsNo, rx_amps, mean_amp, nsyms);
    Somap(bit_likelihood, symbol_likelihood, nsyms);
    for(i=0; i<nsyms*QPSK_BITS_PER_SYMBOL; i++) {
        llr[i] = -bit_likelihood[i];
    }
}

void ldpc_print_info(struct LDPC *ldpc) {
fprintf(stderr, "ldpc->max_iter = %d\n", ldpc->max_iter);
fprintf(stderr, "ldpc->dec_type = %d\n", ldpc->dec_type);
fprintf(stderr, "ldpc->q_scale_factor = %d\n", ldpc->q_scale_factor);
fprintf(stderr, "ldpc->r_scale_factor = %d\n", ldpc->r_scale_factor);
fprintf(stderr, "ldpc->CodeLength = %d\n", ldpc->CodeLength);
fprintf(stderr, "ldpc->NumberParityBits = %d\n", ldpc->NumberParityBits);
fprintf(stderr, "ldpc->NumberRowsHcols = %d\n", ldpc->NumberRowsHcols);
fprintf(stderr, "ldpc->max_row_weight = %d\n", ldpc->max_row_weight);
fprintf(stderr, "ldpc->max_col_weight = %d\n", ldpc->max_col_weight);
fprintf(stderr, "ldpc->data_bits_per_frame = %d\n", ldpc->data_bits_per_frame);
fprintf(stderr, "ldpc->coded_bits_per_frame = %d\n", ldpc->coded_bits_per_frame);
fprintf(stderr, "ldpc->coded_syms_per_frame = %d\n", ldpc->coded_syms_per_frame);
}

/* vi:set ts=4 et sts=4: */
