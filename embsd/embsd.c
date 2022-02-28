/**********************************************************************
	FILE NAME: EMBSD.C
	DEVELOPER: YANG, WONHO
	USAGE: embsd original distorted flag
	where
	embsd : command for running the program
	original : filename of original speech
	distorted : filename of distorted speech
	flag : flag for speech data format
	(0 for MSB-LSB; 1 for LSB-MSB)
	EXAMPLE: mbsd f1.d f1_coder.d 0
**********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "embsd.h"
#include "embsd_s.c"

main( argc, argv )
int argc;
char *argv[];
{
	FILE *fp1, *fp2;    //, *fp3; 									/* FILE POINTERS */
	char *FLAG; 											/* FLAG FOR DATA FORMAT */
	int i;
	int j;
	int p;
	int q;
	int flag;
	int tframe; 											/* TOTAL NUMBER OF FRAMES TO BE PROCESSED */
	double distortion;
	double MBSD; 											/* MBSD VALUE */
	double alpha;
	double pre_MBSD;
	double temp1;
	double pcount;

	if ( argc >= 3 ) { 										/* THERE MUST BE AT LEAST TWO PARAMETERS */
		if (( fp1 = fopen( *++argv, "rb" )) == NULL ) {
			/* OPEN THE ORIGINAL SPEECH FILE */
			printf("ERROR : can't open %s!!\n", *argv );
			return 1;
		}
		if (( fp2 = fopen( *++argv, "rb" )) == NULL ) {
			/* OPEN THE DISTORTED SPEECH FILE */
			printf("ERROR : can't open %s!!\n", *argv );
			return 1;
		}
		else {
//			fp3 = fopen("result.res","a");
			if ( argc == 4 )
				FLAG = *++argv; 							/* FLAG FOR DATA FORMAT */
			else
				FLAG = "1"; 								/* DEFAULT */
			prepare_for_normalization( fp1, fp2, FLAG );
			tframe = 1 + floor( (Nz - FRAME)/(FRAME/2) );	/* TOTAL NUMBER OF FRAMES */
			initialization( fp1, fp2, tframe, FLAG );
			read_original_speech( fp1, FLAG, 2 );			/* READ ONE FRAME OF ORIGINAL SPEECH */
			read_distorted_speech( fp2, FLAG, 2 );			/* READ ONE FRAME OF DISTORTED SPEECH */
			pcount = 0.0;
			distortion = 0.0;
			pre_MBSD = 0.0;
			p = 0;
			q = 0;
			temp1 = 0.0;
			for ( i = 0; i < tframe; i++ ) {
				normalize();
				flag = check_frame();						/* CHECK IF THE FRAME IS TO BE PROCESSED */
				if ( flag == 1 ) { 							/* FOR NON-SILENCE FRAME */
															/* POWER SPECTRUM */
					fft_n01( 0 );
					fft_n01( 1 );							/* BARK SPECTRUM */
					bk_frq( 0 ); 							/* FOR ORIGINAL */
					bk_frq( 1 ); 							/* FOR DISTORTED */
					for ( j = 0; j < 18; j++ ) {
						CX[j] = BX[j];
						CY[j] = BY[j];
					}
															/* BARK SPECTRUM IN PHON SCALE */
					dbtophon( 0 ); 							/* FOR ORIGINAL */
					dbtophon( 1 ); 							/* FOR DISTORTED */
					alpha = sfm();							/* SPECTRAL FLANESS MEASURE*/
					thresh2( alpha );						/* NOISE MASKING THRESHOLD IN dB */
					dbtophon( 2 ); 							/* NOISE MASKING THRESHOLD */
															/* CONVERSION OF PHON LEVEL INTO SONE LEVEL */
					phontoson( 0 ); 						/* FOR ORIGINAL */
					phontoson( 1 ); 						/* FOR DISTORTED */
					phontoson( 2 );							/* FOR NOISE MASKING THRESHOLD */
					MBSD = measure(); 						/* MBSD FOR A FRAME */
															/* COGNITION MODEL */
					p++;
					if ( temp1 < MBSD )
					temp1 = MBSD;

					if ( p == 10 || q > 0 ) {
						pre_MBSD *= T_FACTOR;	
						if ( pre_MBSD < temp1 )
							pre_MBSD = temp1;
						distortion += pre_MBSD;
						p = 0;
						q = 0;
						temp1 = 0.0;
						pcount++;
					}
				}
				else { 										/* FOR SILENCE FRAME */
					q++;
					if ( p > 0 || q == 10 ) {
						pre_MBSD *= T_FACTOR;
						distortion += pre_MBSD;
						p = 0;
						q = 0;
						temp1 = 0.0;
						pcount++;
					}
				}
				read_original_speech( fp1, FLAG, 1 );		/* READ A HALF FRAME OF ORIGINAL SPEECH */
				read_distorted_speech( fp2, FLAG, 1 );		/* READ A HALF FRAME OF DISTORTED SPEECH */
			} /* END OF FOR */
			fclose( fp1 );
			fclose( fp2 );
//			fprintf(fp3, "%5.1f\n", distortion/pcount);
			fprintf(stderr, "EMBSD = %5.1f\n", distortion/pcount);
//			fclose( fp3 );
			return 0;
		} 													/* END OF ELSE */
	} 														/* END OF IF */
	return 1;
} 															/* END OF PROGRAM */
