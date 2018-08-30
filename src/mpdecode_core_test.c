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
#include "mpdecode_core_test.h"
#include "phi0.h"

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

// c_nodes will be an array of NumberParityBits of struct c_node
// Each c_node contains an array of <degree> c_sub_node elements
// This structure reduces the indexing caluclations in SumProduct()

struct c_sub_node { // Order is important here to keep total size small.
  uint16_t index;   // Values from H_rows (except last 2 entries)
  uint16_t socket;  // The socket number at the v_node
  float    message; // modified during operation!
};

struct c_node {
  int degree;       // A count of elements in the following arrays
  struct c_sub_node *subs;
};

// v_nodes will be an array of CodeLength of struct v_node

struct v_sub_node {
  uint16_t index;  //    the index of a c_node it is connected to
                   //    Filled with values from H_cols (except last 2 entries)
  uint16_t socket; //    socket number at the c_node
  float message;   //    Loaded with input data
                   //    modified during operation!
  uint8_t sign;    //    1 if input is negative
                   //    modified during operation!
};

struct v_node {
  int degree;       // A count of ???
  float initial_value;
  struct v_sub_node *subs;
};

// pointers to allocated memory
static struct c_node *c_nodes;
static struct v_node *v_nodes;
static char *DecodedBits;
static int *data_int;

// Variables common to several routines (there are better ways but this works).
static int shift;
static int H1;

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


///////////////////////////////////////
void count_c_v_nodes(
                    struct c_node *c_nodes,
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
                    // Outputs   
                    int     *count_cnodes,
                    int     *count_vnodes) {

    int i, j, k, cnt, count;

    #ifdef PRINT_PROGRESS
    fprintf(stderr, "count_c_v_nodes(*, %d, %d, %d, *, %d, %d, *, %d, *, %d, %d)\n",
                    shift,
                    NumberParityBits,
                    max_row_weight,
                    H1,
                    CodeLength,
                    NumberRowsHcols,
                    max_col_weight,
                    dec_type);
    #endif

    /* first determine the degree of each c-node */

    *count_cnodes = 0;

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
        *count_cnodes += c_nodes[i].degree;
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
        *count_cnodes += c_nodes[i].degree;
        }
    }

    /* determine degree of each v-node */

    *count_vnodes = 0;

    for(i=0;i<(CodeLength-NumberParityBits+shift);i++) {
        count=0;
        for (j=0;j<max_col_weight;j++) {
            if ( H_cols[i+j*NumberRowsHcols] > 0 ) {
                count++;
            }
        }
        v_nodes[i].degree = count;
        *count_vnodes += v_nodes[i].degree;
    }

    for(i=CodeLength-NumberParityBits+shift;i<CodeLength;i++) {
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
        *count_vnodes += v_nodes[i].degree;
    }

    if (shift>0) {
        v_nodes[CodeLength-1].degree =v_nodes[CodeLength-1].degree+1;
        *count_vnodes += 1;
    }
}


///////////////////////////////////////
extern void ldpc_init(struct LDPC *ldpc, int *size_common) {

    int count_cnodes, count_vnodes;

    shift = (ldpc->NumberParityBits + ldpc->NumberRowsHcols) - ldpc->CodeLength;
    if (ldpc->NumberRowsHcols == ldpc->CodeLength) {
        H1=0;
        shift=0;
    } else {
        H1=1;
    }

    #ifdef PRINT_ALLOCS
    fprintf(stderr, "c_nodes = calloc(%d)\n", 
        (int)(ldpc->NumberParityBits * sizeof(struct c_node)));
    #endif
    c_nodes = calloc(ldpc->NumberParityBits, sizeof( struct c_node ) );
    #ifdef PRINT_ALLOCS
    fprintf(stderr, "v_nodes = calloc(%d)\n", 
        (int)(ldpc->CodeLength * sizeof(struct v_node)));
    #endif
    v_nodes = calloc(ldpc->CodeLength, sizeof( struct v_node));

    #ifdef PRINT_ALLOCS
    fprintf(stderr, "DecodedBits = calloc(%d)\n",  
        (int)(ldpc->CodeLength * sizeof(char)));
    #endif
    DecodedBits = calloc(ldpc->CodeLength, sizeof( char ) );

    int DataLength = ldpc->CodeLength - ldpc->NumberParityBits;
    #ifdef PRINT_ALLOCS
    fprintf(stderr, "data_int = calloc(%d)\n", 
        (int)(DataLength * sizeof(int)));
    #endif
    data_int = calloc(DataLength, sizeof(int) );

    count_c_v_nodes(c_nodes, shift,
        ldpc->NumberParityBits, ldpc->max_row_weight, ldpc->H_rows, H1, ldpc->CodeLength,
        v_nodes, ldpc->NumberRowsHcols, ldpc->H_cols, ldpc->max_col_weight, 
        ldpc->dec_type,
        &count_cnodes, &count_vnodes);

    int bytes_c_sub_nodes = 4 * ceil((sizeof(struct c_sub_node) / 4.0));
    int bytes_v_sub_nodes = 4 * ceil((sizeof(struct v_sub_node) / 4.0));

    *size_common = count_cnodes * bytes_c_sub_nodes +
                   count_vnodes * bytes_v_sub_nodes;

    #ifdef PRINT_ALLOCS
    fprintf(stderr, "common memory = %d\n", *size_common);
    #endif
    }


///////////////////////////////////////
void init_c_v_nodes(
                    struct c_node *c_nodes,
                    int     shift,
                    int     NumberParityBits,
                    int     max_row_weight,
                    uint16_t *H_rows,
                    int     H1,
                    int     CodeLength,
                    struct  v_node *v_nodes,
                    int     NumberRowsHcols,
                    uint16_t *H_cols,
                    int     max_col_weight,
                    int     dec_type,
                    void    *mem_common,
                    float   *input) {

    int cnt, count, c_index, v_index;
    int i, j, k;

    struct c_sub_node *cur_c_sub_node = mem_common;
    // TODO: align to word boundary!  Ideally align each one!

    if (H1){

        if (shift ==0){
            for (i=0;i<NumberParityBits;i++) {

                // Assign space from common array
                c_nodes[i].subs = cur_c_sub_node;
                cur_c_sub_node += c_nodes[i].degree;

                for (j=0;j<c_nodes[i].degree-2;j++) {
                    c_nodes[i].subs[j].index = (H_rows[i+j*NumberParityBits] - 1);
                }
                j=c_nodes[i].degree-2;

                if (i==0){
                    c_nodes[i].subs[j].index = (H_rows[i+j*NumberParityBits] - 1);
                }
                else {
                    c_nodes[i].subs[j].index = (CodeLength-NumberParityBits)+i-1;
                }

                j=c_nodes[i].degree-1;
                c_nodes[i].subs[j].index = (CodeLength-NumberParityBits)+i;

            }
        }
        if (shift >0){
            cnt=0;
            for (i=0;i<(NumberParityBits/shift);i++){

                // Assign space from common array
                c_nodes[i].subs = cur_c_sub_node;
                cur_c_sub_node += c_nodes[i].degree;

                for (k =0;k<shift;k++){

                    for (j=0;j<c_nodes[cnt].degree-2;j++) {
                        c_nodes[cnt].subs[j].index = (H_rows[cnt+j*NumberParityBits] - 1);
                    }
                    j=c_nodes[cnt].degree-2;
                    if ((i ==0)||(i==(NumberParityBits/shift-1))){
                        c_nodes[cnt].subs[j].index = (H_rows[cnt+j*NumberParityBits] - 1);
                    }
                    else{
                        c_nodes[cnt].subs[j].index = (CodeLength-NumberParityBits)+k+shift*(i);
                    }
                    j=c_nodes[cnt].degree-1;
                    c_nodes[cnt].subs[j].index = (CodeLength-NumberParityBits)+k+shift*(i+1);
                    if (i== (NumberParityBits/shift-1))
                        {
                            c_nodes[cnt].subs[j].index = (CodeLength-NumberParityBits)+k+shift*(i);
                        }
                    cnt++;
                }
            }
        }

    } else {
        for (i=0;i<NumberParityBits;i++) {

            // Assign space from common array
            c_nodes[i].subs = cur_c_sub_node;
            cur_c_sub_node += c_nodes[i].degree;

            for (j=0;j<c_nodes[i].degree;j++){
                c_nodes[i].subs[j].index = (H_rows[i+j*NumberParityBits] - 1);
            }
        }
    }


    /* set up v_nodes */
    struct v_sub_node *cur_v_sub_node = (struct v_sub_node *)cur_c_sub_node;
    // TODO: align to word boundary!  Ideally align each one!

    for (i=0;i<CodeLength;i++) {

        // Assign space from common array
        v_nodes[i].subs = cur_v_sub_node;
        cur_v_sub_node += v_nodes[i].degree;

        v_nodes[i].initial_value = input[i];

        count=0;

        for (j=0;j<v_nodes[i].degree;j++) {
            if ((H1) && (i>=CodeLength-NumberParityBits+shift)) {
                v_nodes[i].subs[j].index=i-(CodeLength-NumberParityBits+shift)+count;
                if (shift ==0){
                    count=count+1;
                }
                else{
                    count=count+shift;
                }
            } else  {
                v_nodes[i].subs[j].index = (H_cols[i+j*NumberRowsHcols] - 1);
            }

            /* search the connected c-node for the proper message value */
            for (c_index=0;c_index<c_nodes[ v_nodes[i].subs[j].index ].degree;c_index++)
                if ( c_nodes[ v_nodes[i].subs[j].index ].subs[c_index].index == i ) {
                    v_nodes[i].subs[j].socket = c_index;
                    break;
                }

            /* initialize v-node with received LLR */
            if ( dec_type == 1)
                v_nodes[i].subs[j].message = fabs(input[i]);
            else
                v_nodes[i].subs[j].message = phi0( fabs(input[i]) );

            if (input[i] < 0)
                v_nodes[i].subs[j].sign = 1;
        }
    }

    /* now finish setting up the c_nodes */
    for (i=0;i<NumberParityBits;i++) {
        /* index tells which v-nodes this c-node is connected to */
        for (j=0;j<c_nodes[i].degree;j++) {
            /* search the connected v-node for the proper message value */
            for (v_index=0;v_index<v_nodes[ c_nodes[i].subs[j].index ].degree;v_index++)
                if (v_nodes[ c_nodes[i].subs[j].index ].subs[v_index].index == i ) {
                    c_nodes[i].subs[j].socket = v_index;
                    break;
                }
        }
    }

}


/* function for doing the MP decoding */
// Returns the iteration count
int SumProduct(int	  *parityCheckCount,
			 char     DecodedBits[],
			 struct c_node c_nodes[],
			 struct v_node v_nodes[],
			 int	  CodeLength,
			 int	  NumberParityBits,
			 int	  max_iter,
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
//#ifdef __EMBEDDED__
//PROFILE_SAMPLE(ldpc_SP_iter);
//#endif
    for(i=0; i<CodeLength; i++) DecodedBits[i] = 0; // Clear each pass!
    bitErrors = 0;
    /* update r */
//#ifdef __EMBEDDED__
//PROFILE_SAMPLE(ldpc_SP_upr);
//#endif
    ssum = 0;
    for (j=0;j<NumberParityBits;j++) {
      sign = v_nodes[ c_nodes[j].subs[0].index ].subs[ c_nodes[j].subs[0].socket ].sign;
      phi_sum = v_nodes[ c_nodes[j].subs[0].index ].subs[ c_nodes[j].subs[0].socket ].message;

      for (i=1;i<c_nodes[j].degree;i++) {
        // Compiler should optomize this but write the best we can to start from.
        struct c_sub_node *cp = &c_nodes[j].subs[i];
        struct v_sub_node *vp = &v_nodes[ cp->index ].subs[ cp->socket ];
	    phi_sum += vp->message;
	    sign ^= vp->sign;
      }

      if (sign==0) ssum++;

//fprintf(stderr, "  up-r: %d: sign=%d, phi_sum=%f\n", j, sign, (double)phi_sum);

      for (i=0;i<c_nodes[j].degree;i++) {
        struct c_sub_node *cp = &c_nodes[j].subs[i];
        struct v_sub_node *vp = &v_nodes[ cp->index ].subs[ cp->socket ];
	    if ( sign ^ vp->sign ) {
	      cp->message = -phi0( phi_sum - vp->message ); // *r_scale_factor;
	    } else
	      cp->message =  phi0( phi_sum - vp->message ); // *r_scale_factor;
      }
    }

    /* update q */
//#ifdef __EMBEDDED__
//PROFILE_SAMPLE_AND_LOG(ldpc_SP_upq, ldpc_SP_upr, "ldpc_SP_update_r");
//#endif
    for (i=0;i<CodeLength;i++) {

      /* first compute the LLR */
      Qi = v_nodes[i].initial_value;
      for (j=0;j<v_nodes[i].degree;j++) {
        struct v_sub_node *vp = &v_nodes[i].subs[j];
	    Qi += c_nodes[ vp->index ].subs[ vp->socket ].message;
      }

      /* make hard decision */
      if (Qi < 0) {
	    DecodedBits[i] = 1;
      }

      /* now subtract to get the extrinsic information */
      for (j=0;j<v_nodes[i].degree;j++) {
        struct v_sub_node *vp = &v_nodes[i].subs[j];
	    temp_sum = Qi - c_nodes[ vp->index ].subs[ vp->socket ].message;

	    vp->message = phi0( fabs( temp_sum ) ); // *q_scale_factor;
	    if (temp_sum > 0)
	      vp->sign = 0;
	    else
	      vp->sign = 1;
      }
    }
//#ifdef __EMBEDDED__
//PROFILE_SAMPLE_AND_LOG(ldpc_SP_misc, ldpc_SP_upq, "ldpc_SP_update_q");
//#endif

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

//#ifdef __EMBEDDED__
//PROFILE_SAMPLE_AND_LOG2(ldpc_SP_misc, "ldpc_SP_misc");
//PROFILE_SAMPLE_AND_LOG2(ldpc_SP_iter, "ldpc_SP_iter");
//#endif
  }

#ifdef PRINT_PROGRESS
fprintf(stderr, "SumProducts %d iterations\n", result);
#endif
return(result);
}


void ldpc_free_mem(struct LDPC *ldpc) {
    free(c_nodes);
    free(v_nodes);
    free(DecodedBits);
    free(data_int);
}


/* Convenience function to call LDPC decoder from C programs */

int run_ldpc_decoder(struct LDPC *ldpc, void *mem_common, char out_char[], 
                       float input[], int *parityCheckCount) {

    int		    max_iter;
    float       q_scale_factor, r_scale_factor;
    int         CodeLength, NumberParityBits;
    int         i;

    #ifdef PRINT_PROGRESS
    fprintf(stderr, "run_ldpc_decoder()\n");
    #endif

    #ifdef __EMBEDDED__
    PROFILE_VAR(ldpc_init, ldpc_SP);
    PROFILE_SAMPLE(ldpc_init);
    #endif

    /* default values */
    max_iter  = ldpc->max_iter;
    q_scale_factor = ldpc->q_scale_factor;
    r_scale_factor = ldpc->r_scale_factor;

    CodeLength = ldpc->CodeLength;         /* length of entire codeword */
    NumberParityBits = ldpc->NumberParityBits;

    /* initialize c-node and v-node structures */
    init_c_v_nodes( c_nodes, shift,
        ldpc->NumberParityBits, ldpc->max_row_weight, ldpc->H_rows, H1, ldpc->CodeLength,
        v_nodes, ldpc->NumberRowsHcols, ldpc->H_cols, ldpc->max_col_weight, ldpc->dec_type,
        mem_common, input);

    #ifdef __EMBEDDED__
    PROFILE_SAMPLE_AND_LOG(ldpc_SP, ldpc_init, "ldpc_init");
    #endif

    /* Call function to do the actual decoding */

    int iter = SumProduct(parityCheckCount, DecodedBits, c_nodes, v_nodes, 
                CodeLength, NumberParityBits, max_iter, 
                r_scale_factor, q_scale_factor, data_int);

    for (i=0; i<CodeLength; i++) out_char[i] = DecodedBits[i];

    #ifdef __EMBEDDED__
    PROFILE_SAMPLE_AND_LOG2(ldpc_SP, "ldpc_SP");
    #endif

    #ifdef PRINT_PROGRESS
    fprintf(stderr, "parityCheckCount = %d\n", *parityCheckCount);
    fprintf(stderr, "iter = %d\n", iter);
    #endif

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
   2) Using doubles, as experience with FSK in drs232_ldpc.c showed
      doubles were required to obtain the same answers as Octave.
*/

void Demod2D(double  symbol_likelihood[],       /* output, M*number_symbols              */
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


void Somap(double  bit_likelihood[],      /* number_bits, bps*number_symbols */
           double  symbol_likelihood[],   /* M*number_symbols                */
           int     number_symbols)
{
    int    M=QPSK_CONSTELLATION_SIZE, bps = QPSK_BITS_PER_SYMBOL;
    int    n,i,j,k,mask;
    double num[bps], den[bps];
    double metric;

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

            for (k=0;k<bps;k++) {	/* loop over bits */
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

void symbols_to_llrs(double llr[], COMP rx_qpsk_symbols[], float rx_amps[], float EsNo, float mean_amp, int nsyms) {
    int i;

    double symbol_likelihood[nsyms*QPSK_CONSTELLATION_SIZE];
    double bit_likelihood[nsyms*QPSK_BITS_PER_SYMBOL];

    Demod2D(symbol_likelihood, rx_qpsk_symbols, S_matrix, EsNo, rx_amps, mean_amp, nsyms);
    Somap(bit_likelihood, symbol_likelihood, nsyms);
    for(i=0; i<nsyms*QPSK_BITS_PER_SYMBOL; i++) {
        llr[i] = -bit_likelihood[i];
    }
}
