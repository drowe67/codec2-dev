/*
  FILE...: mpdecode_core.c
  AUTHOR.: Matthew C. Valenti, Rohit Iyer Seshadri, David Rowe
  CREATED: Sep 2016

  C-callable core functions moved from MpDecode.c, so they can be used for
  Octave and C programs.
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mpdecode_core.h"

#define QPSK_CONSTELLATION_SIZE 4
#define QPSK_BITS_PER_SYMBOL    2

/* QPSK constellation for symbol likelihood calculations */

static COMP S_matrix[] = {
    { 1.0f,  0.0f},
    { 0.0f,  1.0f},
    { 0.0f, -1.0f},
    {-1.0f,  0.0f}
};
         

int extract_output(char out_char[], int DecodedBits[], int ParityCheckCount[], 
                    int max_iter, int CodeLength, int NumberParityBits);

void encode(struct LDPC *ldpc, unsigned char ibits[], unsigned char pbits[]) {
    unsigned int p, i, tmp, par, prev=0;
    int          ind;
    double      *H_rows = ldpc->H_rows;

    for (p=0; p<ldpc->NumberParityBits; p++) {
        par = 0; 

        for (i=0; i<ldpc->max_row_weight; i++) {
            ind = (int)H_rows[p + i*ldpc->NumberParityBits];
            par = par + ibits[ind-1];
        }

        tmp = par + prev;

        tmp &= 1;    // only retain the lsb 
        prev = tmp; 
        pbits[p] = tmp; 
    }
}

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

static float correction(float xinput )
{
  if (xinput > 2.625 )
    return( 0 );
  else if (xinput < 1 )
    return( -0.375*xinput + 0.6825 );
  else 
    return( -0.1875*xinput + 0.5 );

}

static float LambdaAPPstar(	float mag1,
				float mag2 )
{
  if (mag1 > mag2)
    return( fabs( mag2 + correction( mag1 + mag2 ) - correction( mag1 - mag2 ) ) );
  else
    return( fabs( mag1 + correction( mag1 + mag2 ) - correction( mag2 - mag1 ) ) );
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

void init_c_v_nodes(struct c_node *c_nodes, 
                    int     shift, 
                    int     NumberParityBits, 
                    int     max_row_weight,
                    double *H_rows,
                    int     H1,
                    int     CodeLength,
                    struct v_node *v_nodes, 
                    int     NumberRowsHcols, 
                    double *H_cols,
                    int     max_col_weight,
                    int     dec_type,
                    double *input)
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
void ApproximateMinStar(	 int	  BitErrors[],
				 int      DecodedBits[],
				 struct c_node c_nodes[],
				 struct v_node v_nodes[],
				 int	  CodeLength,
				 int	  NumberParityBits,
				 int	  max_iter )
{
  int i,j, iter;
  int sign;
  float temp_sum;
  float Qi;

  float delta, minval, deltaAPP;
  int mink;

  for (iter=0;iter<max_iter;iter++) {
    /* update r */
    for (j=0;j<NumberParityBits;j++) {	
      /* start new code for approximate-min-star */
      mink = 0;
      sign = v_nodes[ c_nodes[j].index[0] ].sign[ c_nodes[j].socket[0] ];
      minval = v_nodes[ c_nodes[j].index[0] ].message[ c_nodes[j].socket[0] ];
		
      for (i=1;i<c_nodes[j].degree;i++) {
	/* first find the minimum magnitude input message */
	if ( v_nodes[ c_nodes[j].index[i] ].message[ c_nodes[j].socket[i] ] < minval ) {
	  mink = i;
	  minval = v_nodes[ c_nodes[j].index[i] ].message[ c_nodes[j].socket[i] ];							
	}
	/* update the aggregate sign */
	sign ^= v_nodes[ c_nodes[j].index[i] ].sign[ c_nodes[j].socket[i] ];
      }

      /* find the magnitude to send out the minimum input magnitude branch */
      if ( mink == 0 ) {
	delta = v_nodes[ c_nodes[j].index[1] ].message[ c_nodes[j].socket[1] ];
	for (i=2;i<c_nodes[j].degree;i++) {
	  delta = LambdaAPPstar( delta, v_nodes[ c_nodes[j].index[i] ].message[ c_nodes[j].socket[i] ] );
	}
      } else {
	delta = v_nodes[ c_nodes[j].index[0] ].message[ c_nodes[j].socket[0] ];
	for (i=1;i<c_nodes[j].degree;i++) {
	  if ( i != mink )
	    delta = LambdaAPPstar( delta, v_nodes[ c_nodes[j].index[i] ].message[ c_nodes[j].socket[i] ] );
	}
      }

      deltaAPP = LambdaAPPstar( delta, v_nodes[ c_nodes[j].index[mink] ].message[ c_nodes[j].socket[mink] ] );

      /* compute outgoing messages */
      for (i=0;i<c_nodes[j].degree;i++) {
	if ( i == mink ) {
	  if ( sign^v_nodes[ c_nodes[j].index[i] ].sign[ c_nodes[j].socket[i] ] )
	    c_nodes[j].message[i] = - delta;
	  else
	    c_nodes[j].message[i] = delta;
	} else {
	  if ( sign^v_nodes[ c_nodes[j].index[i] ].sign[ c_nodes[j].socket[i] ] )
	    c_nodes[j].message[i] = - deltaAPP;
	  else
	    c_nodes[j].message[i] = deltaAPP;
	}
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
	DecodedBits[iter+max_iter*i] = 1;
	BitErrors[iter]++;
      }

      /* now subtract to get the extrinsic information */
      for (j=0;j<v_nodes[i].degree;j++) {
	temp_sum = Qi - c_nodes[ v_nodes[i].index[j] ].message[ v_nodes[i].socket[j] ];
				
	v_nodes[i].message[j] = fabs( temp_sum );
	if (temp_sum > 0)
	  v_nodes[i].sign[j] = 0;
	else
	  v_nodes[i].sign[j] = 1;
      }
    }

    /* detect errors */
    if (BitErrors[iter] == 0)
      break; 
  }
}


/* function for doing the MP decoding */
void MinSum(		 int	  BitErrors[],
				 int      DecodedBits[],
				 struct c_node c_nodes[],
				 struct v_node v_nodes[],
				 int	  CodeLength,
				 int	  NumberParityBits,
				 int	  max_iter, 
				 float    r_scale_factor,
				 float    q_scale_factor, 
				 int      data[] )
{
  int i,j, iter, i_prime;
  float min_beta;
  int sign;
  float temp_sum;
  float Qi;

  for (iter=0;iter<max_iter;iter++) {

    /* update r */
    for (j=0;j<NumberParityBits;j++) {
      sign = 0;
      for (i=0;i<c_nodes[j].degree;i++) 
	sign ^= v_nodes[ c_nodes[j].index[i] ].sign[ c_nodes[j].socket[i] ];

      for (i=0;i<c_nodes[j].degree;i++) {
	min_beta = 1000;		
								
	for (i_prime=0;i_prime<c_nodes[j].degree;i_prime++) 
	  if ( ( v_nodes[ c_nodes[j].index[i_prime] ].message[c_nodes[j].socket[i_prime]] < min_beta )&&(i_prime != i) )
	    min_beta = v_nodes[ c_nodes[j].index[i_prime] ].message[c_nodes[j].socket[i_prime]];

	if ( sign^v_nodes[ c_nodes[j].index[i] ].sign[ c_nodes[j].socket[i] ] )
	  c_nodes[j].message[i] = -min_beta*r_scale_factor;
	else
	  c_nodes[j].message[i] = min_beta*r_scale_factor;
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
	DecodedBits[iter+max_iter*i] = 1;
      }

      /* now subtract to get the extrinsic information */
      for (j=0;j<v_nodes[i].degree;j++) {
	temp_sum = Qi - c_nodes[ v_nodes[i].index[j] ].message[ v_nodes[i].socket[j] ];
				
	v_nodes[i].message[j] = fabs( temp_sum )*q_scale_factor;
	if (temp_sum > 0)
	  v_nodes[i].sign[j] = 0;
	else
	  v_nodes[i].sign[j] = 1;
      }
    }

    /* count data bit errors, assuming that it is systematic */
    for (i=0;i<CodeLength-NumberParityBits;i++)
      if ( DecodedBits[iter+max_iter*i] != data[i] )
	BitErrors[iter]++;

    /* detect errors */
    if (BitErrors[iter] == 0)
      break; 
  }
}


/* function for doing the MP decoding */
void SumProduct(	 int	  BitErrors[],
			 int      DecodedBits[],
			 struct c_node c_nodes[],
			 struct v_node v_nodes[],
			 int	  CodeLength,
			 int	  NumberParityBits,
			 int	  max_iter,
			 float    r_scale_factor,
			 float    q_scale_factor, 
			 int      data[] )
{
  int i,j, iter;
  float phi_sum;
  int sign;
  float temp_sum;
  float Qi;
  int   ssum; 

  for (iter=0;iter<max_iter;iter++) {
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
	DecodedBits[iter+max_iter*i] = 1;
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
      if ( DecodedBits[iter+max_iter*i] != data[i] )
	BitErrors[iter]++;

    /* Halt if zero errors */
    if (BitErrors[iter] == 0)
      break; 

    // added by Bill -- reuse the BitErrors array to count PCs
    // count the number of PC satisfied and exit if all OK
    BitErrors[iter] = ssum;
    if (ssum==NumberParityBits) break;


  }
   
  // printf(" ssum is %d \n",   ssum); 
}


/* Convenience function to call LDPC decoder from C programs */

int run_ldpc_decoder(struct LDPC *ldpc, char out_char[], double input[], int *parityCheckCount) {
    int		max_iter, dec_type;
    float       q_scale_factor, r_scale_factor;
    int		max_row_weight, max_col_weight;
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

    int *DecodedBits = calloc( max_iter*CodeLength, sizeof( int ) );
    int *ParityCheckCount = calloc( max_iter, sizeof(int) );

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

    int iter = extract_output(out_char, DecodedBits, ParityCheckCount, max_iter, CodeLength, NumberParityBits);

    *parityCheckCount = ParityCheckCount[iter-1];

    /* Clean up memory */

    free(ParityCheckCount);
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


void sd_to_llr(double llr[], double sd[], int n) {
    double sum, mean, sign, sumsq, estvar, estEsN0, x;
    int i;

    /* convert SD samples to LLRs -------------------------------*/

    sum = 0.0;
    for(i=0; i<n; i++)
        sum += fabs(sd[i]);
    mean = sum/n;
                
    /* scale by mean to map onto +/- 1 symbol position */

    for(i=0; i<n; i++) {
        sd[i] /= mean;
    }

    /* find variance from +/-1 symbol position */

    sum = sumsq = 0.0; 
    for(i=0; i<n; i++) {
        sign = (sd[i] > 0.0) - (sd[i] < 0.0);
        x = (sd[i] - sign);
        sum += x;
        sumsq += x*x;
    }
    mean = sum/n;
    estvar = sumsq/n - mean*mean;

    estEsN0 = 1.0/(2.0 * estvar + 1E-3); 
    for(i=0; i<n; i++)
        llr[i] = 4.0 * estEsN0 * sd[i];              
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
    double  tempsr, tempsi, Er, Ei;

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


int extract_output(char out_char[], int DecodedBits[], int ParityCheckCount[], int max_iter, int CodeLength, int NumberParityBits) {
    int i, j;

    /* extract output bits from iteration that solved all parity
       equations, or failing that the last iteration. */

    int converged = 0;
    int iter = 0;
    for (i=0;i<max_iter;i++) {
        if (converged == 0)
            iter++;
        if (ParityCheckCount[i] == NumberParityBits) {
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
    //fprintf(stderr, "iter: %d\n", iter);
    return iter;
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
