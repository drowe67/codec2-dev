/**********************************************************************
 FILE NAME: EMBSD_S.C
 DEVELOPER: YANG, WONHO
**********************************************************************/

void hanning_window()
/* THIS FUNCTION CALCULATES HANNING WINDOW */
{
    extern double W[FRAME];
    int i;

    for ( i = 0; i < FRAME; i++ )
        W[i] = 0.5*(1.0-cos(2.0*PI*(i+1.0)/(FRAME+1.0)));
}

void check_original_speech1( fp )
FILE *fp;
/* THIS FUNCTION READS AN ORIGINAL BINARY SPEECH FILE
 AND FIND OUT THE NUMBER OF SAMPLES IN THAT FILE */
{
    extern double Nx; /* NUMBER OF SAMPLES IN ORIGINAL SPEECH */
    int t;
    double k;

    k = 0.0;
    while( !feof( fp ) ) {
        t = getc( fp ); /* GET 2 BYTES */
        t = getc( fp );
        k++;
    }
    Nx = k - (double)OFFSET;
    rewind( fp );
}

void check_distorted_speech1( fp )
FILE *fp;
/* THIS FUNCTION READS AN ORIGINAL BINARY SPEECH FILE
 AND FIND OUT THE NUMBER OF SAMPLES IN THAT FILE */
{
    extern double Ny;
    int t;
    double k;

    k = 0.0;
    while( !feof( fp ) ) {
        t = getc( fp ); /* GET 2 BYTES */
        t = getc( fp );
        k++;
    }
    Ny = k - (double)OFFSET;
    rewind( fp );
}

int read_speech_sample( fp, FLAG )
FILE *fp;
char *FLAG;
/* THIS FUNCTION READS A SPEECH SAMPLE FROM A FILE */
{
    int MSB, LSB, sign, n, n1, t;
    int check = 0x00ff;

    if ( *FLAG == '0' ) { /* MSB-LSB FORMAT */
        MSB = getc( fp ); /* GET ONE BYTE */
        LSB = getc( fp ); /* GET ONE BYTE */
        sign = MSB;
        sign = sign >> 7;
        if ( sign == 0 ) /* POSITIVE */
            n = MSB;
        else { /* NEGATIVE */
            t = ~MSB;
            n = t & check;
            n = -1 * n;
        }
        if ( sign == 1 ) { /* NEGATIVE */
            t = ~LSB;
            n1 = t & check;
            n1 = -1 * n1 - 1;
            n = n * 256 + n1;
        }
        else /* POSITIVE */
            n = n * 256 + LSB;
    } /* END OF IF */
    else { /* LSB-MSB FORMAT */
        LSB = getc( fp ); /* GET ONE BYTE */
        MSB = getc( fp ); /* GET ONE BYTE */
        sign = MSB;
        sign = sign >> 7;
        if ( sign == 0 ) /* POSITIVE */
            n = MSB;
        else { /* NEGATIVE */
            t = ~MSB;
            n = t & check;
            n = -1 * n;
        }
        if ( sign == 1 ) { /* NEGATIVE */
            t = ~LSB;
            n1 = t & check;
            n1 = -1 * n1 - 1;
            n = n * 256 + n1;
        }
        else /* POSITIVE */
            n = n * 256 + LSB;
    } /* END OF ELSE */
    return n;
}

void check_original_speech2( fp, FLAG )
FILE *fp;
char *FLAG;
/* THIS FUNCTION READS A BINARY SPEECH FILE AND FIND OUT
 DC OFFSET OF THE SPEECH SIGNAL */
{
    extern double XMEAN;
    extern double Nz;
    int n;
    double k;
    double temp1 = 0.0;

    k = 0.0;
    while( k < Nz + (double)OFFSET ) {
        n = read_speech_sample( fp, FLAG );
        if ( k >= OFFSET )
        temp1 += (double)n; /* SUM */
        k++;
    } /* END OF WHILE */
    XMEAN = temp1 / ( k - (double)OFFSET ); /* MEAN */
    rewind( fp );
}

void check_distorted_speech2( fp, FLAG )
FILE *fp;
char *FLAG;
/* THIS FUNCTION READS A BINARY SPEECH FILE AND FIND OUT
 DC OFFSET OF THE SPEECH SIGNAL */
{
 extern double YMEAN;
 extern double Nz;
 int n;
 double k;
 double temp1 = 0.0;

    k = 0.0;
    while( k < Nz + (double)OFFSET ) {
        n = read_speech_sample( fp, FLAG );
        if ( k >= OFFSET )
            temp1 += (double)n; /* SUM */
        k++;
    } /* END OF WHILE */
    YMEAN = temp1 / ( k - (double)OFFSET ); /* MEAN */
    rewind( fp );
}

void find_original_rms( fp, FLAG )
FILE *fp;
char *FLAG;
/* THIS FUNCTION READS A BINARY SPEECH FILE AND FIND OUT
 RMS VALUE OF THE SPEECH SIGNAL */
{
 extern double XMEAN; /* DC OFFSET OF ORIGINAL SPEECH */
 extern double XRMS; /* RMS VALUE OF ORIGINAL SPEECH */
 extern double Nz;
 int n;
 double k;
 double temp1;
 double temp2 = 0.0;

    k = 0.0;
    while( k < Nz + (double)OFFSET ) {
        n = read_speech_sample( fp, FLAG );
        if ( k >= OFFSET ) {
            temp1 = (double)n - XMEAN;
            temp2 += temp1 * temp1;
        }
        k++;
    } /* END OF WHILE */
    XRMS = sqrt(temp2 /( k - (double)OFFSET));
    rewind( fp );
}

void find_distorted_rms( fp, FLAG )
FILE *fp;
char *FLAG;
/* THIS FUNCTION READS A BINARY SPEECH FILE AND FIND OUT
 RMS VALUE OF THE SPEECH SIGNAL */
{
 extern double YMEAN; /* DC OFFSET OF DISTORTED SPEECH */
 extern double YRMS; /* RMS VALUE OF DISTORTED SPEECH */
 extern double Nz;
 int n;
 double k;
 double temp1;
 double temp2 = 0.0;

    k = 0;
    while( k < Nz + (double)OFFSET ) {
        n = read_speech_sample( fp, FLAG );
        if ( k >= OFFSET ) {
            temp1 = (double)n - YMEAN;
            temp2 += temp1 * temp1;
        }
        k++;
    } /* END OF WHILE */
    YRMS = sqrt(temp2 / ( k - (double)OFFSET ));
    rewind( fp );
}

void read_header( fp1, fp2 )
FILE *fp1;
FILE *fp2;
/* THIS FUNCTION READS HEADER OF BINARY SPEECH FILES */
{
 int t;
 int k;

    k = 0;
    while( k < OFFSET ) {
        t = getc( fp1 ); /* GET ONE BYTE */
        t = getc( fp1 ); /* GET ONE BYTE */
        t = getc( fp2 ); /* GET ONE BYTE */
        t = getc( fp2 ); /* GET ONE BYTE */
        k++;
    } /* END OF WHILE */
}

void read_original_speech( fp, FLAG, p )
FILE *fp;
char *FLAG;
int p; /* p = 1 FOR READING REAR HALF FRAME
 p = 2 FOR READING A FRAME */
/* THIS PROGRAM READS A BINARY SPEECH FILE IN WHICH
A SAMPLE IS A 2 BYTE INTEGER AS AN INPUT AND WRITES
THOSE INTEGERS. THESE 2 BYTES ARE STORED IN MSB-LSB
OR LSB-MSB. IF SAMPLES ARE STORED IN MSB-LSB, FLAG
SHOULD BE "0". OTHERWISE, FLAG SHOULD BE "1".
IF flag IS 0, THIS PROGRAM READS THE ORIGINAL SPEECH.
IF flag IS 1, THIS PROGRAM READS THE DISTORTED SPEECH. */
{
 extern int X[FRAME];
 int n;
 int k;
 int i;

    k = 0;
    if ( p == 1 ) /* READING HALF FRAME */
    for ( i = 0; i < FRAME/2; i++ ) /* OVERLAPPED HALF FRAME */
        X[i] = X[i+FRAME/2];
    while( k < p * (FRAME/2) ) {
        n = read_speech_sample( fp, FLAG );
        if ( p == 1 )
            X[(FRAME/2)+k] = n;
/* STORE A SPEECH SAMPLES IN AN ARRAY */
        else
            X[k] = n;
        k++;
    } /* END OF WHILE */
}

void read_distorted_speech( fp, FLAG, p )
FILE *fp;
char *FLAG;
int p; /* p = 1 FOR READING REAR HALF FRAME
 p = 2 FOR READING A FRAME */
/* THIS PROGRAM READS A BINARY SPEECH FILE IN WHICH
A SAMPLE IS A 2 BYTE INTEGER AS AN INPUT AND WRITES
THOSE INTEGERS. THESE 2 BYTES ARE STORED IN MSB-LSB
OR LSB-MSB. IF SAMPLES ARE STORED IN MSB-LSB, FLAG
SHOULD BE "0". OTHERWISE, FLAG SHOULD BE "1".
IF flag IS 0, THIS PROGRAM READS THE ORIGINAL SPEECH.
IF flag IS 1, THIS PROGRAM READS THE DISTORTED SPEECH. */
{
 extern int Y[FRAME];
 int n;
 int k;
 int i;

    k = 0;
    if ( p == 1 )
        for ( i = 0; i < FRAME/2; i++ )
            Y[i] = Y[i+FRAME/2];
    while( k < p * (FRAME/2) ) {
        n = read_speech_sample( fp, FLAG );
        if ( p == 1 )
            Y[(FRAME/2)+k] = n;
/* STORE A SPEECH SAMPLES IN AN ARRAY */
        else
            Y[k] = n;
        k++;
    } /* END OF WHILE */
}

void normalize()
/* THIS FUNCTION NORMALIZE TWO INPUT SIGNALS */
{
extern int X[FRAME]; /* ORIGINAL SPEECH */
extern int Y[FRAME]; /* DISTORTED SPEECH */
extern double XX[FRAME]; /* NORMALIZED ORIGINAL SPEECH */
extern double YY[FRAME]; /* NORMALIZED DISTORTED SPEECH */
extern double XMEAN;
extern double YMEAN;
extern double XRMS;
extern double YRMS;
int i;

    for ( i = 0; i < FRAME; i++ ) {
        XX[i] = (double)X[i] - XMEAN;
        YY[i] = (double)Y[i] - YMEAN;
    }
    for ( i = 0; i < FRAME; i++ ) {
        XX[i] = NORM * XX[i] / XRMS;
        YY[i] = NORM * YY[i] / YRMS;
    }
}

void silence_threshold( fp1, fp2, tframe, FLAG )
FILE *fp1;
FILE *fp2;
int tframe;
char *FLAG;
/* THIS FUNCTION DETERMINES THE THRESHOLD FOR A SILENCE FRAME */
{
 extern double W[FRAME];
 extern double XX[FRAME];
 extern double YY[FRAME];
 extern double XTHRESHOLD; /* SILENCE THRESHOLD FOR PROCESSING */
 extern double YTHRESHOLD; /* SILENCE THRESHOLD FOR PROCESSING */
 int i, j;
 double xenergy, max_xenergy;
 double yenergy, max_yenergy;

    read_header( fp1, fp2 );
    max_xenergy = 0.0;
    max_yenergy = 0.0;
    read_original_speech( fp1, FLAG, 2 );
    read_distorted_speech( fp2, FLAG, 2 );
    for ( j = 0; j < tframe; j++ ) {
        normalize();
        xenergy = 0.0;
        yenergy = 0.0;
        for ( i = 0; i < FRAME; i++ ) {
            xenergy += (XX[i] * W[i])*(XX[i] * W[i]);
            yenergy += (YY[i] * W[i])*(YY[i] * W[i]);
        }
        if ( xenergy > max_xenergy )
            max_xenergy = xenergy;
        if ( yenergy > max_yenergy )
            max_yenergy = yenergy;
        read_original_speech( fp1, FLAG, 1 );
        read_distorted_speech( fp2, FLAG, 1 );
    }
    XTHRESHOLD = pow(10.0, -1.5) * max_xenergy; /* 15dB BELOW */
    YTHRESHOLD = pow(10.0, -3.5) * max_yenergy; /* 35dB BELOW */
    rewind( fp1 );
    rewind( fp2 );
}

double sfm()
/* USING POWER SPECTRUM OF ORIGINAL SPEECH,
 THIS FUNCTION COMPUTES THE SPECTRAL FLATNESS MEASURE.
for alpha = 1, SFM <= -60 dB : entirely tone-like signal
for alpha = 0, SFM >= 0 dB : entirely noise-like signal
*/
{
 extern double PSX[FSIZE]; /* POWER SPECTRUM OF ORIGINAL */
 double alpha;
 double a_mean; /* ALGEBRAIC MEAN */
 double g_mean; /* GEOMETRIC MEAN */
 int i;
 double sum1, sum2;
 double sfm_db, sfm_db_ratio;
 double t;

    sum1 = 0.0;
    sum2 = 0.0;
    for ( i = 0; i < FSIZE; i++ ) {
        sum1 += PSX[i];
        sum2 += log10( PSX[i] );
    }
    a_mean = sum1 / (double)FSIZE;
    t = sum2 / (double)FSIZE;
    g_mean = pow( 10.0, t );
    sfm_db = 10.0 * log10( g_mean / a_mean );
    sfm_db_ratio = sfm_db / -60.0;
    if ( sfm_db_ratio < 1 )
        alpha = sfm_db_ratio;
    else
        alpha = 1;
    return alpha;
}

void thresh2( alpha )
double alpha;
/* USING SPREAD BARK SPECTRUM IN DB, THIS FUNCTION
CALCULATES NOISE MASKING THRESHOLD */
{
 extern double CX[BSIZE]; /* SPREAD BARK SPECTRUM OF ORIGINAL */
 extern double CNMT[BSIZE];
/* NOISE MASKING THRESHOLD IN SPREAD SPECTRUM */
 extern double ABS_TH[BSIZE];
 int i;
 double t, tt;

    for ( i = 0; i < BSIZE; i++ ) {
        t = alpha * ( 14.5 + (double)i + 1.0 ) + ( 1.0 - alpha ) * 5.5;
        tt = 10.0 * log10(CX[i]) - t;
        if ( tt < ABS_TH[i] )
            CNMT[i] = pow(10.0, ABS_TH[i]/10.0);
        else
            CNMT[i] = pow(10.0, tt/10.0);
    }
}

/* fft_real.c** Routines for split-radix, real-only transforms.
These routines are adapted from [Sorenson 1987] * * When all x[j] are
real the standard DFT of (x[0],x[1],...,x[N-1]),* call it x^, has the
property of Hermitian symmetry: x^[j] =x^[N-j].
Thus we only need to find the set (x^[0].re, x^[1].re,..., x^[N/2].re,
x^[N/2-1].im, ..., x^[1].im) * which, like the original signal x, has N
elements.* The two key routines perform forward (real-to-Hermitian)
FFT, and * backward (Hermitian-to-real) FFT, respectively. For example,
the* sequence: fft_real_to_hermitian(x, N);
fftinv_hermitian_to_real(x, N); is an identity operation on the
signal x. To convolve twopure-real signals x and y, one does:
fft_real_to_hermitian(x, N);fft_real_to_hermitian(y, N);
mul_hermitian(y, x, N);fftinv_hermitian_to_real(x, N); and x is the
pure-real cyclic convolution of x and y. */
void init_sincos()
{
 extern int cur_run;
 extern double *sncos;
 int j;
 double e = TWOPI / N;

    if ( N <= cur_run )
        return;
    cur_run = N;
    if ( sncos )
        free( sncos );
    sncos = (double *)malloc(sizeof(double) * ( 1 + ( N >> 2 )));
    for ( j = 0; j <= ( N >> 2 ); j++ )
        sncos[j] = sin( e * j );
}

double s_sin( n )
int n;
{
 extern int cur_run;
 extern double *sncos;
 int seg = n / (cur_run >> 2);

    switch (seg) {
    case 0:
        return (sncos[n]);
    case 1:
        return (sncos[(cur_run >> 1) - n]);
    case 2:
        return (-sncos[n - (cur_run >> 1)]);
    case 3:
        return (-sncos[cur_run - n]);
    }
}

double s_cos( n )
int n;
{
 extern int cur_run;
 int quart = (cur_run >> 2);
 
    if (n < quart)
        return (s_sin(n + quart));
    return (-s_sin(n - quart));
}

void scramble_real( x )
double *x;
{
 register int i, j, k;
 double tmp;

    for ( i = 0, j = 0; i < N - 1; i++ ) {
        if ( i < j ) {
            tmp = x[j];
            x[j] = x[i];
            x[i] = tmp;
        }   
        k = N / 2;
        while ( k <= j ) {
            j -= k;
            k >>= 1;
        }
        j += k;
    }
}

void fft_real_to_hermitian( z )
double *z;
/*
 * Output is {Re(z^[0]),...,Re(z^[n/2),Im(z^[n/2-1]),...,Im(z^[1]).
 * This is a decimation-in-time, split-radix algorithm.
 */
{
 extern int cur_run;
 register double cc1, ss1, cc3, ss3;
 register int is, id, i0, i1, i2, i3, i4, i5, i6, i7, i8, a, a3, b, b3, nminus = N - 1, dil, expand;
 register double *x, e;
 int nn = N >> 1;
 double t1, t2, t3, t4, t5, t6;
 register int n2, n4, n8, i, j;

    init_sincos();
    expand = cur_run / N;
    scramble_real( z );
    x = z - 1; /* FORTRAN compatibility. */
    is = 1;
    id = 4;
    do {
        for ( i0 = is; i0 <= N; i0 += id) {
            i1 = i0 + 1;
            e = x[i0];
            x[i0] = e + x[i1];
            x[i1] = e - x[i1];
        }
        is = ( id << 1 ) - 1;
        id <<= 2;
    } while ( is < N );
    n2 = 2;
    while ( nn >>= 1 ) {
        n2 <<= 1;
        n4 = n2 >> 2;
        n8 = n2 >> 3;
        is = 0;
        id = n2 << 1;
        do {
            for ( i = is; i < N; i += id ) {
                i1 = i + 1;
                i2 = i1 + n4;
                i3 = i2 + n4;
                i4 = i3 + n4;
                t1 = x[i4] + x[i3];
                x[i4] -= x[i3];
                x[i3] = x[i1] - t1;
                x[i1] += t1;
                if ( n4 == 1 )
                    continue;
                i1 += n8;
                i2 += n8;
                i3 += n8;
                i4 += n8;
                t1 = ( x[i3] + x[i4] ) * SQRTHALF;
                t2 = ( x[i3] - x[i4] ) * SQRTHALF;
                x[i4] = x[i2] - t1;
                x[i3] = -x[i2] - t1;
                x[i2] = x[i1] - t2;
                x[i1] += t2;
            }
            is = (id << 1) - n2;
            id <<= 2;
        } while ( is < N);
        dil = N / n2;
        a = dil;
        for ( j = 2; j <= n8; j++ ) {
            a3 = ( a + ( a << 1 )) & nminus;
            b = a * expand;
            b3 = a3 * expand;
            cc1 = s_cos(b);
            ss1 = s_sin(b);
            cc3 = s_cos(b3);
            ss3 = s_sin(b3);
            a = (a + dil) & nminus;
            is = 0;
            id = n2 << 1;
            do {
                for ( i = is; i < N; i += id ) {
                    i1 = i + j;
                    i2 = i1 + n4;
                    i3 = i2 + n4;
                    i4 = i3 + n4;
                    i5 = i + n4 - j + 2;
                    i6 = i5 + n4;
                    i7 = i6 + n4;
                    i8 = i7 + n4;
                    t1 = x[i3] * cc1 + x[i7] * ss1;
                    t2 = x[i7] * cc1 - x[i3] * ss1;
                    t3 = x[i4] * cc3 + x[i8] * ss3;
                    t4 = x[i8] * cc3 - x[i4] * ss3;
                    t5 = t1 + t3;
                    t6 = t2 + t4;
                    t3 = t1 - t3;
                    t4 = t2 - t4;
                    t2 = x[i6] + t6;
                    x[i3] = t6 - x[i6];
                    x[i8] = t2;
                    t2 = x[i2] - t3;
                    x[i7] = -x[i2] - t3;
                    x[i4] = t2;
                    t1 = x[i1] + t5;
                    x[i6] = x[i1] - t5;
                    x[i1] = t1;
                    t1 = x[i5] + t4;
                    x[i5] -= t4;
                    x[i2] = t1;
                }
                is = (id << 1) - n2;
                id <<= 2;
            } while ( is < N );
        } /* END OF for */
    } /* END OF while */
} /* END OF function */

void fft_n01( flag )
/* CALCULATE POWER SPECTRUM
 IF flag IS 0, CALCULATE POWER SPECTRUM OF ORIGINAL SPEECH
 IF flag IS 1, CALCULATE POWER SPECTRUM OF DISTORTED SPEECH */
int flag;
{
 extern double W[FRAME]; /* HANNING WINDOW */
 extern double FREQ[FSIZE]; /* FREQUENCY SCALE */
 extern double XX[FRAME]; /* NORMALIZED ORIGINAL SPEECH */
 extern double YY[FRAME]; /* NORMALIZED DISTORTED SPEECH */
 extern double PSX[FSIZE]; /* POWER SPECTRUM OF ORIGINAL */
 extern double PSY[FSIZE]; /* POWER SPECTRUM OF DISTORTED */
 int i;
 double xxa[N];
 double x[N];
 double t;

    if ( flag == 0 )
        for ( i = 0; i < FRAME; i++ )
            x[i] = XX[i] * W[i];
    else
        for ( i = 0; i < FRAME; i++ )
            x[i] = YY[i] * W[i];
    for ( i = FRAME; i < N; i++ )
        x[i] = 0.0;
    fft_real_to_hermitian( x );
    for ( i = 0; i < N; i++ ){
        if ( i == 0 ) /* || i == FSIZE/2 ) */
            xxa[i] = x[i] * x[i] / (double)N;
        else
            xxa[i] = ( x[i]*x[i] + x[N-i]*x[N-i] ) / (double)N;
        if ( i > 0 )
            xxa[i] *= 2.0;
    }
    for ( i = 0; i < FSIZE; i++ ) {
        t = 8000.0/ (double)N;
        FREQ[i] = i * t;
        if ( flag == 0 )
            PSX[i] = xxa[i];
        else
            PSY[i] = xxa[i];
    }
}

void bk_frq( flag )
int flag;
/* Computes Critcal Bands in the Bark Spectrum */
{
 extern int BARK[BSIZE+1];
 extern double FREQ[FSIZE];
 extern double PSX[FSIZE]; /* POWER SPECTRUM OF ORIGINAL */
 extern double PSY[FSIZE]; /* POWER SPECTRUM OF DISTORTED */
 extern double BX[BSIZE]; /* BARK SPECTRUM OF ORIGINAL */
 extern double BY[BSIZE]; /* BARK SPECTRUM OF DISTORTED */
 int i,j;

    if ( flag == 0 ) {
        for ( i = 0; i < BSIZE; i++ )
            BX[i] = 0.0;

        for ( i = 0; i < BSIZE; i++ )
            for( j = 0; j < FSIZE; j++ )
                if( BARK[i] <= FREQ[j] && FREQ[j] < BARK[i+1] )
/* redo this freq j */
        BX[i] += PSX[j];
    }
    else {
        for ( i = 0; i < BSIZE; i++ )
            BY[i] = 0.0;
        for ( i = 0; i < BSIZE; i++ )
            for( j = 0; j < FSIZE; j++ )
                if( BARK[i] <= FREQ[j] && FREQ[j] < BARK[i+1] )
/* redo this freq j */
        BY[i] += PSY[j];
    }
}

void thrshld()
/* Estimate the threshold of hearing in dB by the formula of Terhardt
thrshld(f) = { 3.64(f/1000)^(-0.8) - 6.5exp[-0.6(f/1000 - 3.3)^2]
 + 0.001(f/1000)^4 }
 This Formula produces threshold of hearing in dB
Reference : Terhardt, E., Stoll, G. and Seewann, M, "Algorithm for
extraction of pitch and pitch salience from complex
tonal signals", J. Acoust. Soc. Am., vol. 71(3), Mar., 1982
*/
{
 extern double Abs_thresh[BSIZE];
/* ABSOLUTE HEARING THRESHOLD IN BARK */
 extern int BARK[BSIZE+1]; /* BARK FREQUENCY */
 extern double FREQ[FSIZE]; /* FREQUENCY SCALE */
 int k = 0;
 int i;
 int j;
 double f;
 double L[FSIZE];
 double xox, xox1, xox2, SUM;

    SUM = 0.0;
    for( i = 0; i < FSIZE-1; i++ ) {
        f = FREQ[i+1]/1000.0;
        xox = f * f;
        xox *= xox;
        xox = 0.001 * xox;
        xox1 = pow( f, 0.8);
        xox1 = 3.64 / xox1;
        xox2 = f - 3.3;
        xox2 *= xox2;
        xox2 = .6 * xox2;
        xox2 = -1.0 * xox2;
        xox2 = 6.5 * exp( xox2 );
        L[i+1] = xox1 - xox2 + xox;
    }
    L[0] = 0.0;
    for( i = 1; i <= 18; i++ ){
        for( j = 1; j < FSIZE; j++ ){
            if ( BARK[i-1] <= FREQ[j] && FREQ[j] < BARK[i] ){
                SUM += L[j];
                k++;
            }
            else {
                SUM = 0.0;
                k = 1;
            }
        }
        Abs_thresh[i-1] = SUM / k;
    }
}

void dbtophon( flag )
/* CONVERT SPREAD BARK SPECTRUM INTO PHON SCALE */
int flag;
{
 extern double CX[BSIZE]; /* SPREAD BARK SPECTRUM OF ORIGINAL */
 extern double CY[BSIZE]; /* SPREAD BARK SPECTRUM OF DISTORTED */
 extern double CNMT[BSIZE];
 /* NOISE MASKING THRESHOLD IN SPREAD BARK SPECTRUM */
 extern double PX[BSIZE-3];
 /* SPREAD BARK SPECTRUM OF ORIGINAL IN PHON SCALE */
 extern double PY[BSIZE-3];
 /* SPREAD BARK SPECTRUM OF DISTORTED IN PHON SCALE */
 extern double PN[BSIZE-3];
 /* SPREAD BARK SPECTRUM OF NOISE IN PHON SCALE */
 int i;
 int j;
 double t1;
 double T[BSIZE-3]={0};

    if ( flag == 0 ) { /* FOR ORIGINAL SPEECH */
        for( i = 0; i < BSIZE-3; i++ )
            T[i] = 10.0 * log10( CX[i] );
        for( i = 0; i < BSIZE-3; i++ ){
            j = 0;
            while( T[i] >= eqlcon[j][i] )
                j++;
            if( j == BSIZE-3 )
                break;
            if( j == 0 )
                PX[i] = phons[0];
            else {
                t1 = ( T[i] - eqlcon[j-1][i] ) / ( eqlcon[j][i] - eqlcon[j-1][i] );
                PX[i] = phons[j-1] + t1 * (phons[j] - phons[j-1]);
            }
        }
    }
    else if ( flag == 1 ) { /* FOR DISTORTED SPEECH */
        for( i = 0; i < BSIZE-3; i++ )
            T[i] = 10.0 * log10( CY[i] );
        for( i = 0; i < BSIZE-3; i++ ){
            j = 0;
            while( T[i] >= eqlcon[j][i] )
                j++;
            if( j == BSIZE-3 )
                break;
            if( j == 0 )
                PY[i] = phons[0];
            else {
                t1 = ( T[i] - eqlcon[j-1][i] ) / ( eqlcon[j][i] - eqlcon[j-1][i] );
                PY[i] = phons[j-1] + t1 * (phons[j] - phons[j-1]);
            }
        }
    }
    else { /* FOR NOISE MASKING THRESHOLD */
        for( i = 0; i < BSIZE-3; i++ )
            T[i] = 10.0 * log10( CNMT[i] );
        for( i = 0; i < BSIZE-3; i++ ){
            j = 0;
            while( T[i] >= eqlcon[j][i] )
                j++;
            if( j == BSIZE-3 )
                break;
            if( j == 0 )
                PN[i] = phons[0];
            else {
                t1 = ( T[i] - eqlcon[j-1][i] ) / ( eqlcon[j][i] - eqlcon[j-1][i] );
                PN[i] = phons[j-1] + t1 * (phons[j] - phons[j-1]);
            }
        }
    }
}

void phontoson( flag )
/* CONVERT LOUDNESS LEVEL (PHON SCALE) INTO LOUDNESS (SONE SCALE) */
int flag;
{
extern double PX[BSIZE-3];
/* SPREAD BARK SPECTRUM OF ORIGINAL IN PHON SCALE */
extern double PY[BSIZE-3];
/* SPREAD BARK SPECTRUM OF DISTORTED IN PHON SCALE */
extern double PN[BSIZE-3];
/* SPREAD BARK SPECTRUM OF NOISE IN PHON SCALE */
extern double SX[BSIZE-3];
/* SPECIFIC LOUDNESS OF ORIGINAL */
extern double SY[BSIZE-3];
/* SPECIFIC LOUDNESS OF DISTORTED */
extern double SN[BSIZE-3];
/* SPECIFIC LOUDNESS OF NOISE */
int i;
double xox;

    if ( flag == 0 ) { /* FOR ORIGINAL SPEECH */
        for( i = 0; i < BSIZE-3; i++ )
            if( PX[i] >= 40.0 ){
                xox = PX[i] - 40.0;
                xox *= 0.1;
                SX[i] = pow( 2.0, xox );
            }
            else{
                xox = PX[i] / 40.0;
                SX[i] = pow( xox, 2.642 );
            }
    }      
    else if ( flag == 1 ) { /* FOR DISTORTED SPEECH */
        for( i = 0; i < BSIZE-3; i++ )
            if( PY[i] >= 40.0 ){
                xox = PY[i] - 40.0;
                xox *= 0.1;
                SY[i] = pow( 2.0, xox );
            }
            else{
                xox = PY[i] / 40.0;
                SY[i] = pow( xox, 2.642 );
            }
        }
    else { /* FOR NOISE MASKING THRESHOLD */
        for( i = 0; i < BSIZE-3; i++ )
            if( PN[i] >= 40.0 ){
                xox = PN[i] - 40.0;
                xox *= 0.1;
                SN[i] = pow( 2.0, xox );
            }
            else{
                xox = PN[i] / 40.0;
                SN[i] = pow( xox, 2.642 );
            }
    }
}

double measure()
/* METRIC ESTIMATING DISTORTION */
{
extern double SX[BSIZE-3];
/* SPECIFIC LOUDNESS OF ORIGINAL */
extern double SY[BSIZE-3];
/* SPECIFIC LOUDNESS OF DISTORTED */
extern double SN[BSIZE-3];
/* SPECIFIC LOUDNESS OF NOISE MASKING THRESHOLD */
 int i;
 double dist = 0.0;
 double temp;
 double x;
 double ww;
 double w1;
 double temp1, temp2;
 extern double WEIGHT[15];

    temp1 = 1.0;
    temp2 = 1.0;
    for ( i = 0; i < 15; i++ ) {
        temp1 += SX[i];
        temp2 += SY[i];
    }
    w1 = temp1 / temp2;
    dist = 0.0;
    for ( i = 0; i < BSIZE-3; i++ ) {
        ww = w1;
        temp = fabs( SX[i] - ww*SY[i] );
        x = temp - SN[i];
        if ( x > 0.0 )
            dist += WEIGHT[i]*x;
    }
    return dist;
}

void prepare_for_normalization( fp1, fp2, FLAG )
FILE *fp1;
FILE *fp2;
char *FLAG;
{
extern double Nx; /* NUMBER OF SAMPLES OF ORIGINAL */
extern double Ny; /* NUMBER OF SAMPLES OF DISTORTED */
extern double Nz; /* NUMBER OF SAMPLES TO BE COMPARED */

    check_original_speech1( fp1 );
    check_distorted_speech1( fp2 );
    if ( Nx < Ny )
        Nz = Nx;
    else
        Nz = Ny;
    check_original_speech2( fp1, FLAG );
    check_distorted_speech2( fp2, FLAG );
    find_original_rms( fp1, FLAG );
    find_distorted_rms( fp2, FLAG );
}

void initialization( fp1, fp2, tframe, FLAG )
FILE *fp1;
FILE *fp2;
int tframe;
char *FLAG;
{
int i;
double t;
extern double FREQ[FSIZE];

    for ( i = 0; i < FSIZE; i++ ) {
        t = 8000.0/ (double)N;
        FREQ[i] = i * t;
    }
    hanning_window(); /* HANNING WINDOW */
    thrshld(); /* ABSOLUTE HEARING THRESHOLD */
    silence_threshold( fp1, fp2, tframe, FLAG );
}

int check_frame()
{
extern double W[FRAME];
extern double XX[FRAME]; /* NORMALIZED ORIGINAL SPEECH */
extern double YY[FRAME]; /* NORMALIZED DISTORTED SPEECH */
extern double XTHRESHOLD; /* SILENCE THRESHOLD FOR PROCESSING */
extern double YTHRESHOLD; /* SILENCE THRESHOLD FOR PROCESSING */
double xenergy;
double yenergy;
int i;
int flag;

    xenergy = 0.0;
    yenergy = 0.0;
    for ( i = 0; i < FRAME; i++ ) {
        xenergy += (XX[i] * W[i])*(XX[i] * W[i]);
        yenergy += (YY[i] * W[i])*(YY[i] * W[i]);
    }
    if ( xenergy > XTHRESHOLD && yenergy > YTHRESHOLD )
        flag = 1;
    else
        flag = 0;
    return flag;
}
