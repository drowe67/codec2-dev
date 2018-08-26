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

struct c_sub_node {
  uint16_t index;   // Values from H_rows (except last 2 entries)
  float    message; // modified during operation!
  uint16_t socket;  // The socket number at the v_node
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

extern void ldpc_init(void) {
    }

/*******************************************************
 * Macros for a signed scaled integer 
 * 16 bits of integer and 16 bits of fraction in 32 bits.
#define SI_16p16(f) ((int)(f * (1 << 16)))
 */

/* Phi function */
static float phi0( float x ) {

  //int x_si = SI_16p16(x);

  if (x > 10.0f) return( 0.0f ); else
    if (x > 2.375000f) {
      if (x > 3.812500f) {
        if (x > 4.562500f) {
          if (x > 4.937500f) {
            if (x > 7.000000f) {
              if (x > 8.000000f) {
                if (x > 9.000000f) {
                  return(0.00016879f);
                } else {
                  return(0.00045882f);
                }
              } else {
                return(0.0012472f);
              }
            } else {
              if (x > 5.000000f) {
                if (x > 6.000000f) {
                  return(0.00339025f);
                } else {
                  return(0.0131598f);
                }
              } else {
                return(0.0140086f);
              }
            }
          } else {
            if (x > 4.750000f) {
              if (x > 4.812500f) {
                if (x > 4.875000f) {
                  return(0.0149121f);
                } else {
                  return(0.0158739f);
                }
              } else {
                return(0.0168977f);
              }
            } else {
              if (x > 4.625000f) {
                if (x > 4.687500f) {
                  return(0.0179875f);
                } else {
                  return(0.0191477f);
                }
              } else {
                return(0.0203827f);
              }
            }
          }
        } else {
          if (x > 4.187500f) {
            if (x > 4.375000f) {
              if (x > 4.437500f) {
                if (x > 4.500000f) {
                  return(0.0216974f);
                } else {
                  return(0.0230969f);
                }
              } else {
                return(0.0245866f);
              }
            } else {
              if (x > 4.250000f) {
                if (x > 4.312500f) {
                  return(0.0261725f);
                } else {
                  return(0.0278607f);
                }
              } else {
                return(0.0296578f);
              }
            }
          } else {
            if (x > 4.000000f) {
              if (x > 4.062500f) {
                if (x > 4.125000f) {
                  return(0.0315709f);
                } else {
                  return(0.0336074f);
                }
              } else {
                return(0.0357753f);
              }
            } else {
              if (x > 3.875000f) {
                if (x > 3.937500f) {
                  return(0.0380832f);
                } else {
                  return(0.04054f);
                }
              } else {
                return(0.0431554f);
              }
            }
          }
        }
      } else {
        if (x > 3.062500f) {
          if (x > 3.437500f) {
            if (x > 3.625000f) {
              if (x > 3.687500f) {
                if (x > 3.750000f) {
                  return(0.0459396f);
                } else {
                  return(0.0489036f);
                }
              } else {
                return(0.052059f);
              }
            } else {
              if (x > 3.500000f) {
                if (x > 3.562500f) {
                  return(0.0554182f);
                } else {
                  return(0.0589944f);
                }
              } else {
                return(0.0628016f);
              }
            }
          } else {
            if (x > 3.250000f) {
              if (x > 3.312500f) {
                if (x > 3.375000f) {
                  return(0.0668549f);
                } else {
                  return(0.0711702f);
                }
              } else {
                return(0.0757646f);
              }
            } else {
              if (x > 3.125000f) {
                if (x > 3.187500f) {
                  return(0.0806561f);
                } else {
                  return(0.0858642f);
                }
              } else {
                return(0.0914094f);
              }
            }
          }
        } else {
          if (x > 2.687500f) {
            if (x > 2.875000f) {
              if (x > 2.937500f) {
                if (x > 3.000000f) {
                  return(0.0973138f);
                } else {
                  return(0.103601f);
                }
              } else {
                return(0.110296f);
              }
            } else {
              if (x > 2.750000f) {
                if (x > 2.812500f) {
                  return(0.117425f);
                } else {
                  return(0.125017f);
                }
              } else {
                return(0.133104f);
              }
            }
          } else {
            if (x > 2.500000f) {
              if (x > 2.562500f) {
                if (x > 2.625000f) {
                  return(0.141716f);
                } else {
                  return(0.150889f);
                }
              } else {
                return(0.160662f);
              }
            } else {
              if (x > 2.437500f) {
                return(0.171072f);
              } else {
                return(0.182165f);
              }
            }
          }
        }
      }
    } else {
      if (x > 0.875000f) {
        if (x > 1.625000f) {
          if (x > 2.000000f) {
            if (x > 2.187500f) {
              if (x > 2.250000f) {
                if (x > 2.312500f) {
                  return(0.193985f);
                } else {
                  return(0.206583f);
                }
              } else {
                return(0.22001f);
              }
            } else {
              if (x > 2.062500f) {
                if (x > 2.125000f) {
                  return(0.234326f);
                } else {
                  return(0.249591f);
                }
              } else {
                return(0.265873f);
              }
            }
          } else {
            if (x > 1.812500f) {
              if (x > 1.875000f) {
                if (x > 1.937500f) {
                  return(0.283243f);
                } else {
                  return(0.301781f);
                }
              } else {
                return(0.321571f);
              }
            } else {
              if (x > 1.687500f) {
                if (x > 1.750000f) {
                  return(0.342706f);
                } else {
                  return(0.365288f);
                }
              } else {
                return(0.389428f);
              }
            }
          }
        } else {
          if (x > 1.250000f) {
            if (x > 1.437500f) {
              if (x > 1.500000f) {
                if (x > 1.562500f) {
                  return(0.415249f);
                } else {
                  return(0.442887f);
                }
              } else {
                return(0.472491f);
              }
            } else {
              if (x > 1.312500f) {
                if (x > 1.375000f) {
                  return(0.504229f);
                } else {
                  return(0.538291f);
                }
              } else {
                return(0.574888f);
              }
            }
          } else {
            if (x > 1.062500f) {
              if (x > 1.125000f) {
                if (x > 1.187500f) {
                  return(0.614262f);
                } else {
                  return(0.65669f);
                }
              } else {
                return(0.702491f);
              }
            } else {
              if (x > 0.937500f) {
                if (x > 1.000000f) {
                  return(0.752038f);
                } else {
                  return(0.80577f);
                }
              } else {
                return(0.864208f);
              }
            }
          }
        }
      } else {
        if (x > 0.125000f) {
          if (x > 0.500000f) {
            if (x > 0.687500f) {
              if (x > 0.750000f) {
                if (x > 0.812500f) {
                  return(0.927985f);
                } else {
                  return(0.997872f);
                }
              } else {
                return(1.07483f);
              }
            } else {
              if (x > 0.562500f) {
                if (x > 0.625000f) {
                  return(1.16009f);
                } else {
                  return(1.25524f);
                }
              } else {
                return(1.36239f);
              }
            }
          } else {
            if (x > 0.312500f) {
              if (x > 0.375000f) {
                if (x > 0.437500f) {
                  return(1.48447f);
                } else {
                  return(1.6257f);
                }
              } else {
                return(1.79241f);
              }
            } else {
              if (x > 0.187500f) {
                if (x > 0.250000f) {
                  return(1.9949f);
                } else {
                  return(2.25157f);
                }
              } else {
                return(2.60048f);
              }
            }
          }
        } else {
          if (x > 0.001953f) {
            if (x > 0.015625f) {
              if (x > 0.031250f) {
                if (x > 0.062500f) {
                  return(3.14427f);
                } else {
                  return(3.83695f);
                }
              } else {
                return(4.52999f);
              }
            } else {
              if (x > 0.003906f) {
                if (x > 0.007812f) {
                  return(5.2231f);
                } else {
                  return(5.91624f);
                }
              } else {
                return(6.60939f);
              }
            }
          } else {
            if (x > 0.000244f) {
              if (x > 0.000488f) {
                if (x > 0.000977f) {
                  return(7.30254f);
                } else {
                  return(7.99568f);
                }
              } else {
                return(8.68883f);
              }
            } else {
              if (x > 0.000122f) {
                return(9.38198f);
              } else {
                return(10.0f);
              }
            }
          }
        }
      }
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

void alloc_c_v_nodes(struct c_node *c_nodes,
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
                    int     dec_type) {
    int i, j, k, count, cnt, c_index, v_index;

    #ifdef PRINT_PROGRESS
    fprintf(stderr, "alloc_c_v_nodes(*, %d, %d, %d, *, %d, %d, *, %d, *, %d, %d)\n",
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
                #ifdef PRINT_ALLOCS
                if (i<4) fprintf(stderr, 
                  "c_node.subs[%d].* = calloc(%d)\n", i,
                  (int)(c_nodes[i].degree * sizeof(struct c_sub_node)));
                if (i==4) fprintf(stderr, "...\n");
                #endif
                c_nodes[i].subs = calloc(c_nodes[i].degree, 
                                         sizeof(struct c_sub_node));

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

                for (k =0;k<shift;k++){
                #ifdef PRINT_ALLOCS
                fprintf(stderr, 
                  "cnodes[%d].* = calloc(%d)\n", cnt,
                  (int)(c_nodes[cnt].degree * sizeof(struct c_sub_node)));
                #endif
                    c_nodes[cnt].subs =  calloc( c_nodes[cnt].degree, sizeof( struct c_sub_node ) );

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
            /* now that we know the size, we can dynamically allocate memory */
            #ifdef PRINT_ALLOCS
            fprintf(stderr, 
              "cnodes[%d].* = calloc(%d)\n", cnt,
              (int)(c_nodes[cnt].degree * sizeof(struct c_sub_node)));
            #endif
            c_nodes[i].subs = calloc(c_nodes[i].degree, 
                                     sizeof( struct c_sub_node ) );
            for (j=0;j<c_nodes[i].degree;j++){
                c_nodes[i].subs[j].index = (H_rows[i+j*NumberParityBits] - 1);
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
        #ifdef PRINT_ALLOCS
        if (i<4) fprintf(stderr, 
          "cnodes[%d].* = calloc(%d)\n", i,
          (int)(v_nodes[i].degree * sizeof(struct v_sub_node)));
        if (i==4) fprintf(stderr, "...\n");
        #endif
        v_nodes[i].subs = calloc( v_nodes[i].degree, sizeof( struct v_sub_node ) );

        /* index tells which c-nodes this v-node is connected to */
        count=0;

        for (j=0;j<v_nodes[i].degree;j++) {
            if ((H1)&& (i>=CodeLength-NumberParityBits+shift)){
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

/*
//  Optional dump for development and debugging
void ldpc_dump_nodes(struct LDPC *ldpc) {

    int i, j;

    int NumberParityBits = ldpc->NumberParityBits;
    int CodeLength = ldpc->CodeLength;

    FILE *f = fopen("ldpc_node_dump.txt", "w");

    fprintf(f, "Cnodes:\n");
    fprintf(f, "  size: %d\n", NumberParityBits);
    fprintf(f, "  array:\n");
    for (i=0;i<NumberParityBits;i++) {
        fprintf(f, "    -\n");
        fprintf(f, "      degree: %d\n", c_nodes[i].degree);
        fprintf(f, "      index:\n");
        for (j=0;j<c_nodes[i].degree;j++)
            fprintf(f, "        - %d\n", c_nodes[i].index[j]);
        fprintf(f, "      socket:\n");
        for (j=0;j<c_nodes[i].degree;j++)
            fprintf(f, "        - %d\n", c_nodes[i].socket[j]);
    }

    fprintf(f, "\n");
    fprintf(f, "Vnodes:\n");
    fprintf(f, "  size: %d\n", CodeLength);
    fprintf(f, "  array:\n");
    for (i=0;i<CodeLength;i++) {
        fprintf(f, "    -\n");
        fprintf(f, "      degree: %d\n", v_nodes[i].degree);
        fprintf(f, "      index:\n");
        for (j=0;j<v_nodes[i].degree;j++)
            fprintf(f, "        - %d\n", v_nodes[i].index[j]);
        fprintf(f, "      socket:\n");
        for (j=0;j<v_nodes[i].degree;j++)
            fprintf(f, "        - %d\n", v_nodes[i].socket[j]);
    }

    fclose(f);
}
*/

void init_v_nodes(int    CodeLength,
                  struct v_node *v_nodes,
                  int    dec_type,
                  float  *input) {

    int i, j;

    /* set up v_nodes */
    for (i=0;i<CodeLength;i++) {

        v_nodes[i].initial_value = input[i];

        for (j=0;j<v_nodes[i].degree;j++) {
            /* initialize v-node with received LLR */
            if ( dec_type == 1)
                v_nodes[i].subs[j].message = fabs(input[i]);
            else
                v_nodes[i].subs[j].message = phi0( fabs(input[i]) );

            if (input[i] < 0)
                v_nodes[i].subs[j].sign = 1;
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

// Memory initialization, call when modem is setup
void ldpc_alloc_mem(struct LDPC *ldpc) {

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

    alloc_c_v_nodes(c_nodes, shift,
        ldpc->NumberParityBits, ldpc->max_row_weight, ldpc->H_rows, H1, ldpc->CodeLength,
        v_nodes, ldpc->NumberRowsHcols, ldpc->H_cols, ldpc->max_col_weight, 
        ldpc->dec_type);
}

void ldpc_free_mem(struct LDPC *ldpc) {
    int i;

    /*  Cleaning c-node elements */
    for (i=0;i<ldpc->NumberParityBits;i++) {
        free( c_nodes[i].subs );
    }
    free(c_nodes);

    /* Cleaning v-node elements */
    for (i=0;i<ldpc->CodeLength;i++) {
        free( v_nodes[i].subs);
    }
    free(v_nodes);

    free(DecodedBits);
    free(data_int);
}


/* Convenience function to call LDPC decoder from C programs */

int run_ldpc_decoder(struct LDPC *ldpc, char out_char[], float input[], int *parityCheckCount) {
    int		max_iter, dec_type;
    float       q_scale_factor, r_scale_factor;
    int         CodeLength, NumberParityBits;
    int         i;

    #ifdef PRINT_PROGRESS
    fprintf(stderr, "run_ldpc_decoder()\n");
    #endif

#ifdef __EMBEDDED__
PROFILE_VAR(ldpc_init, ldpc_SP);
#endif

    /* default values */
    max_iter  = ldpc->max_iter;
    dec_type  = ldpc->dec_type;
    q_scale_factor = ldpc->q_scale_factor;
    r_scale_factor = ldpc->r_scale_factor;

    CodeLength = ldpc->CodeLength;         /* length of entire codeword */
    NumberParityBits = ldpc->NumberParityBits;

#ifdef __EMBEDDED__
PROFILE_SAMPLE(ldpc_init);
#endif

    /* initialize c-node and v-node structures */
    init_v_nodes(CodeLength, v_nodes, dec_type, input);

#ifdef __EMBEDDED__
PROFILE_SAMPLE_AND_LOG(ldpc_SP, ldpc_init, "ldpc_init");
#endif

    /* Call function to do the actual decoding */

    int iter = SumProduct(parityCheckCount, DecodedBits, c_nodes, v_nodes, 
                CodeLength, NumberParityBits, max_iter, 
                r_scale_factor, q_scale_factor, data_int);

    for (i=0; i<CodeLength; i++) out_char[i] = DecodedBits[i];

    #ifdef PRINT_PROGRESS
    fprintf(stderr, "parityCheckCount = %d\n", *parityCheckCount);
    fprintf(stderr, "iter = %d\n", iter);
    #endif

#ifdef __EMBEDDED__
PROFILE_SAMPLE_AND_LOG2(ldpc_SP, "ldpc_SP");
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
