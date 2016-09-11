/*
  FILE...: mpdecode_core.c
  AUTHOR.: Matthew C. Valenti, Rohit Iyer Seshadri, David Rowe
  CREATED: Sep 2016

  C-callable core functions moved from MpDecode.c, so they can be used for
  Octave and C programs.
*/

#include <math.h>
#include <stdlib.h>
#include "mpdecode_core.h"

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

static float correction(
			float xinput )
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


