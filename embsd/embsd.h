/**********************************************************************
 FILE NAME: EMBSD.H
 DEVELOPER: YANG, WONHO
**********************************************************************/
#include    <stdio.h>
#include    <stdlib.h>
#include    <math.h>
#define     FRAME       320                 /* FRAME SIZE IN SAMPLES */
#define     PI          3.14159265358979323846
#define     NORM        1000.0              /* NORM AMPLITUDE */
#define     BSIZE       18                  /* NUMBER OF BARK FREQUENCIES */
#define     FSIZE       512                 /* HALF OF FFT SIZE */
#define     N           1024                /* FFT SIZE */
#define     TWOPI       (2*3.14159265358979323846)
#define     SQRTHALF    0.70710678118654752440
#define     OFFSET      0                   /* HEADER LENGTH IN BYTES */
#define     T_FACTOR    0.8
double      XMEAN;                          /* DC OFFSET OF ORIGINAL SPEECH */
double      YMEAN;                          /* DC OFFSET OF DISTORTED SPEECH */
double      XRMS;                           /* RMS VALUE OF ORIGINAL SPEECH */
double      YRMS;                           /* RMS VALUE OF DISTORTED SPEECH */
double      XTHRESHOLD;                     /* SILENCE THRESHOLD FOR PROCESSING */
double      YTHRESHOLD;                     /* SILENCE THRESHOLD FOR PROCESSING */
double      W[FRAME];                       /* HANNING WINDOW */
double      FREQ[FSIZE];                    /* FREQUENCY SCALE */
double      Abs_thresh[BSIZE];              /* ABSOLUTE HEARING THRESHOLD IN BARK */
int         X[FRAME];                       /* ORIGINAL SPEECH */
int         Y[FRAME];                       /* DISTORTED SPEECH */
double      XX[FRAME];                      /* NORMALIZED ORIGINAL SPEECH */
double      YY[FRAME];                      /* NORMALIZED DISTORTED SPEECH */
double      PSX[FSIZE];                     /* POWER SPECTRUM OF ORIGINAL */
double      PSY[FSIZE];                     /* POWER SPECTRUM OF DISTORTED */
double      BX[BSIZE];                      /* BARK SPECTRUM OF ORIGINAL */
double      BY[BSIZE];                      /* BARK SPECTRUM OF DISTORTED */
double      CX[BSIZE];                      /* SPREAD BARK SPECTRUM OF ORIGINAL */
double      CX1[BSIZE];                     /* SPREAD BARK SPECTRUM FOR NMT */
double      CY[BSIZE];                      /* SPREAD BARK SPECTRUM OF DISTORTED */
double      PX[BSIZE-3];                    /* SPREAD BARK SPECTRUM OF ORIGINAL IN PHON SCALE */
double      PY[BSIZE-3];                    /* SPREAD BARK SPECTRUM OF DISTORTED IN PHON SCALE */
double      PN[BSIZE-3];                    /* SPREAD BARK SPECTRUM OF NOISE IN PHON SCALE */
double      SX[BSIZE-3];                    /* SPECIFIC LOUDNESS OF ORIGINAL */
double      SY[BSIZE-3];                    /* SPECIFIC LOUDNESS OF DISTORTED */
double      SN[BSIZE-3];                    /* SPECIFIC LOUDNESS OF NOISE */
double      CNMT[BSIZE];                    /* NOISE MASKING THRESHOLD IN SPREAD BARK SPECTRUM */
double      ABS_TH[BSIZE];                  /* ABSOLUTE HEARING THRESHOLD */
double      Nx;                             /* NUMBER OF SAMPLES IN ORIGINAL SPEECH */
double      Ny;                             /* NUMBER OF SAMPLES IN DISTORTED SPEECH */
double      Nz;                             /* NUMBER OF SAMPLES TO BE COMPARED */
int         cur_run = 0;
double      *sncos = NULL;
double      WEIGHT[15] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
int         BARK[BSIZE+1] = {0, 100, 200, 300, 400, 510, 630, 770, 920,
            1080, 1270,1480, 1720, 2000, 2320, 2700,
            3150, 3700, 4400 };             /* BARK FREQUENCY */
double      eqlcon[13][15] =                /* EQUAL-LOUDNESS CONTOURS */
            {   {12,7,4,1,0,0,0,-0.5,-2,-3,-7,-8,-8.5,-8.5,-8.5},
                {20,17,14,12,10,9.5,9,8.5,7.5,6.5,4,3,2.5,2,2.5},
                {29,26,23,21,20,19.5,19.5,19,18,17,15,14,13.5,13,13.5},
                {36,34,32,30,29,28.5,28.5,28.5,28,27.5,26,25,24.5,24,24.5},
                {45,43,41,40,40,40,40,40,40,39.5,38,37,36.5,36,36.5},
                {53,51,50,49,48.5,48.5,49,49,49,49,48,47,46.5,45.5,46},
                {62,60,59,58,58,58.5,59,59,59,59,58,57.5,57,56,56},
                {70,69,68,67.5,67.5,68,68,68,68,68,67,66,65.5,64.5,64.5},
                {79,79,79,79,79,79,79,79,78,77.5,76,75,74.5,73,73},
                {89,89,89,89.5,90,90,90,89.5,89,88.5,87,86,85.5,84,83.5},
                {100,100,100,100,100,99.5,99,99,98.5,98,96,95,94.5,93.5,93},
                {112,112,112,112,111,110.5,109.5,109,108.5,108,106,105,104.5,103,102.5},
                {122,122,121,121,120.5,120,119,118,117,116.5,114.5,113.5,113,111,110.5} 
            };
double      phons[13]=                      /* LOUDNESS LEVELS (PHON SCALES) */
            {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0};

/* FUNCTIONS */

void hanning_window( void );
void check_original_speech1( FILE * );
void check_distorted_speech1( FILE * );
int read_speech_sample( FILE *, char * );
void check_original_speech2( FILE *, char * );
void check_distorted_speech2( FILE *, char * );
void find_original_rms( FILE *, char * );
void find_distorted_rms( FILE *, char * );
void read_header( FILE *, FILE * );
void read_original_speech( FILE *, char *, int );
void read_distorted_speech( FILE *, char *, int );
void normalize( void );
void silence_threshold( FILE *, FILE *, int, char * );
double sfm( void );
void thresh2( double );
void init_sincos( void );
double s_sin( int );
double s_cos( int );
void scramble_real( double * );
void fft_real_to_hermitian( double * );
void fft_n01( int );
void bk_frq( int );
void thrshld( void );
void dbtophon( int );
void phontoson( int );
double measure( void );
void prepare_for_normalization( FILE *, FILE *, char * );
void initialization( FILE *, FILE *, int, char * );
int check_frame( void );
