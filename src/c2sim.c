/*---------------------------------------------------------------------------*\

  FILE........: c2sim.c
  AUTHOR......: David Rowe
  DATE CREATED: 20/8/2010

  Codec2 simulation.  Combines encoder and decoder and allows
  switching in and out various algorithms and quantisation steps. Used
  for algorithm development.  Various option to write/read vectors to
  external files for plotting and external processing/VQ training.
  Quite a lot of cruft and dead code over the years.  It isn't pretty
  - but that's OK!

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2009 David Rowe

  All rights reserved.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License version 2.1, as
  published by the Free Software Foundation.  This program is
  distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this program; if not, see <http://www.gnu.org/licenses/>.
*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>

#include "defines.h"
#include "sine.h"
#include "nlp.h"
#include "dump.h"
#include "lpc.h"
#include "lsp.h"
#include "quantise.h"
#include "phase.h"
#include "postfilter.h"
#include "interp.h"
#include "bpf.h"
#include "bpfb.h"
#include "newamp1.h"
#include "lpcnet_freq.h"
#include "sd.h"

void synth_one_frame(int n_samp, codec2_fft_cfg fft_inv_cfg, short buf[],
                     MODEL *model, float Sn_[], float Pn[], int prede, float *de_mem, float gain);
void print_help(const struct option *long_options, int num_opts, char* argv[]);

#define N_SAMP n_samp  /* quick fix for run time sample rate selection */

extern int phase_dispersion;

struct RATEKMEANCOMP {
    float noise_gate;
    float thresh;
    float slope;
    float makeup_gain;
};

/*---------------------------------------------------------------------------* \

				MAIN

\*---------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{

    int Fs = 8000;
    int set_fs;

    int lpc_model = 0, order = LPC_ORD;
    int lsp = 0, lspd = 0, lspvq = 0;
    int lspjmv = 0;
    int prede = 0;
    int   postfilt;
    int   hand_voicing = 0, hi = 0, simlpcpf = 0, modelin=0, modelout=0;
    int   lpcpf = 0;
    FILE *fvoicing = 0;
    int dec;
    int decimate = 1;
    int   amread, Woread, pahw;
    int   awread;
    int   hmread;
    int   phase0 = 0;
    int   scalar_quant_Wo_e = 0;
    int   scalar_quant_Wo_e_low = 0;
    int   vector_quant_Wo_e = 0;
    int   dump_pitch_e = 0;
    float gainoutlin = 1.0;
    float comp_limit_samples = 32767;
    float comp_gain_dB = 0.0;
    float comp_dr_lin = 0.0;
    int   comp_en = 0;
    int   bpf_en = 0;
    int   bpfb_en = 0;
    int   hpf200_en = 0;
    FILE *fam = NULL, *fWo = NULL;
    FILE *faw = NULL;
    FILE *fhm = NULL;
    FILE *fjmv = NULL;
    FILE *flspEWov = NULL;
    FILE *ften_ms_centre = NULL;
    FILE *fmodelout = NULL;
    FILE *fmodelin = NULL;
    #ifdef DUMP
    int   dump;
    #endif
    char  out_file[MAX_STR];
    FILE *fout = NULL;	/* output speech file */
    int   rateK = 0, newamp1vq = 0, rate_K_dec = 0, perframe=0, postfilter_newamp1_en=0, rateKlf=0, setK=0;
    float rateK_mean_min = -1E32, rateK_mean_max = 1E32;
    int   rateK_mean_min_en, rateK_mean_max_en, rateK_mean_comp_en;
    struct RATEKMEANCOMP rateKmeancomp = {0.0,0.0,0.0,0.0};
    int   K = 20, K_st=0;
    int   K_en=K-1;
    float framelength_s = N_S;
    int   lspEWov = 0, rateKWov = 0, first = 0;
    FILE  *frateKWov = NULL;
    int   ten_ms_centre = 0;
    FILE  *fphasenn = NULL;
    FILE  *frateK = NULL;
    FILE  *frateKnomean = NULL;
    FILE  *frateKin = NULL;
    FILE  *frateKnomeanin = NULL;
    FILE  *frateKmean = NULL;
    int   rateKout, rateKnomeanout, rateKmeanout_en, rateKin, rateKnomeanin;
    
    char* opt_string = "ho:";
    struct option long_options[] = {
        { "Fs", required_argument, &set_fs, 1 },
        { "rateK", no_argument, &rateK, 1 },
        { "setK", required_argument, &setK, 1 },
        { "K_st", required_argument, NULL, 0 },
        { "K_en", required_argument, NULL, 0 },
        { "perframe", no_argument, &perframe, 1 },
        { "newamp1vq", no_argument, &newamp1vq, 1 },
        { "rateKdec", required_argument, &rate_K_dec, 1 },
        { "rateKout", required_argument, &rateKout, 1 },
        { "rateKnomeanout", required_argument, &rateKnomeanout, 1 },
        { "rateKmeanout", required_argument, &rateKmeanout_en, 1 },
        { "rateKin", required_argument, &rateKin, 1 },
        { "rateKnomeanin", required_argument, &rateKnomeanin, 1 },
        { "rateK_mean_min", required_argument, &rateK_mean_min_en, 1 },
        { "rateK_mean_max", required_argument, &rateK_mean_max_en, 1 },
        { "rateK_mean_comp", required_argument, &rateK_mean_comp_en, 1 },
        { "postfilter_newamp1", no_argument, &postfilter_newamp1_en, 1 },
        { "rateKLF", no_argument, &rateKlf, 1 },
        { "lpc", required_argument, &lpc_model, 1 },
        { "lsp", no_argument, &lsp, 1 },
        { "lspd", no_argument, &lspd, 1 },
        { "lspvq", no_argument, &lspvq, 1 },
        { "lspjmv", no_argument, &lspjmv, 1 },
        { "phase0", no_argument, &phase0, 1 },
        { "dispersion", required_argument, &phase_dispersion, 1 },
        { "postfilter", no_argument, &postfilt, 1 },
        { "hand_voicing", required_argument, &hand_voicing, 1 },
        { "dec", required_argument, &dec, 1 },
        { "hi", no_argument, &hi, 1 },
        { "simlpcpf", no_argument, &simlpcpf, 1 },
        { "lpcpf", no_argument, &lpcpf, 1 },
        { "prede", no_argument, &prede, 1 },
        { "dump_pitch_e", required_argument, &dump_pitch_e, 1 },
        { "sq_pitch_e", no_argument, &scalar_quant_Wo_e, 1 },
        { "sq_pitch_e_low", no_argument, &scalar_quant_Wo_e_low, 1 },
        { "vq_pitch_e", no_argument, &vector_quant_Wo_e, 1 },
        { "rate", required_argument, NULL, 0 },
        { "gainoutlin", required_argument, NULL, 0 },
        { "comp", required_argument, NULL, 0 },
        { "comp_gain", required_argument, NULL, 0 },
        { "comp_dr", required_argument, NULL, 0 },
        { "bpf", no_argument, &bpf_en, 1 },
        { "bpfb", no_argument, &bpfb_en, 1 },
        { "hpf200", no_argument, &hpf200_en, 1 },
        { "amread", required_argument, &amread, 1 },
        { "hmread", required_argument, &hmread, 1 },
        { "awread", required_argument, &awread, 1 },
        { "Woread", required_argument, &Woread, 1 },
        { "pahw", required_argument, &pahw, 1 },
        { "lspEWov", required_argument, &lspEWov, 1 },
        { "rateKWov", required_argument, &rateKWov, 1 },
        { "first", no_argument, &first, 1 },
        { "ten_ms_centre", required_argument, &ten_ms_centre, 1 },
        { "framelength_s", required_argument, NULL, 0 },
        { "modelout",  required_argument, &modelout, 1 },
        { "modelin",   required_argument, &modelin, 1 },
        #ifdef DUMP
        { "dump", required_argument, &dump, 1 },
        #endif
        { "help", no_argument, NULL, 'h' },
        { NULL, no_argument, NULL, 0 }
    };
    int num_opts=sizeof(long_options)/sizeof(struct option);

    /*----------------------------------------------------------------*\

                     Interpret Command Line Arguments

    \*----------------------------------------------------------------*/

    if (argc < 2) {
        print_help(long_options, num_opts, argv);
    }

    while(1) {
        int option_index = 0;
        int opt = getopt_long(argc, argv, opt_string,
                    long_options, &option_index);
        if (opt == -1)
            break;
        switch (opt) {
         case 0:
            if(strcmp(long_options[option_index].name, "Fs") == 0) {
                Fs= atoi(optarg);
                if((Fs != 8000) && (Fs != 16000)) {
                    fprintf(stderr, "Error Fs must be 8000 or 16000\n");
                    exit(1);
                }
            } else if(strcmp(long_options[option_index].name, "lpc") == 0) {
                order = atoi(optarg);
            } else if(strcmp(long_options[option_index].name, "dispersion") == 0) {
                phase_dispersion = atoi(optarg);
            #ifdef DUMP
            } else if(strcmp(long_options[option_index].name, "dump") == 0) {
                if (dump)
	            dump_on(optarg);
            #endif
            } else if(strcmp(long_options[option_index].name, "lsp") == 0
                  || strcmp(long_options[option_index].name, "lspd") == 0
                  || strcmp(long_options[option_index].name, "lspvq") == 0) {
	        assert(order == LPC_ORD);
            } else if(strcmp(long_options[option_index].name, "setK") == 0) {
                K = atoi(optarg);
                fprintf(stderr, "K: %d\n", K);
                K_en = K-1;
            } else if(strcmp(long_options[option_index].name, "K_st") == 0) {
                K_st = atoi(optarg);
                fprintf(stderr, "K_st: %d\n", K_st);
            } else if(strcmp(long_options[option_index].name, "K_en") == 0) {
                K_en = atoi(optarg);
                fprintf(stderr, "K_en: %d\n", K_en);
            } else if(strcmp(long_options[option_index].name, "rateKdec") == 0) {
                rate_K_dec = atoi(optarg);
                fprintf(stderr, "rate_K_dec: %d\n", rate_K_dec);
	    } else if(strcmp(long_options[option_index].name, "rateKout") == 0) {
                /* write rateK vectors to file */
                if ((frateK = fopen(optarg,"wb")) == NULL) {
	            fprintf(stderr, "Error opening output rateK file: %s: %s\n",
		        optarg, strerror(errno));
                    exit(1);
                }
                fprintf(stderr, "rateKout: each record is %d bytes\n", (int)(K*sizeof(float)));
	    } else if(strcmp(long_options[option_index].name, "rateKnomeanout") == 0) {
                /* write rateK mean removed vectors to file */
                if ((frateKnomean = fopen(optarg,"wb")) == NULL) {
	            fprintf(stderr, "Error opening output rateKnomean file: %s: %s\n",
		        optarg, strerror(errno));
                    exit(1);
                }
                fprintf(stderr, "rateKnomeanout: each record is %d bytes\n", (int)(K*sizeof(float)));
	    } else if(strcmp(long_options[option_index].name, "rateKmeanout") == 0) {
                /* write rateK mean scalar to file */
                if ((frateKmean = fopen(optarg,"wb")) == NULL) {
	            fprintf(stderr, "Error opening output rateKmean file: %s: %s\n",
		        optarg, strerror(errno));
                    exit(1);
                }
                fprintf(stderr, "rateKmeanout: each record is %d bytes\n", (int)sizeof(float));
	    } else if(strcmp(long_options[option_index].name, "rateKin") == 0) {
                /* read rateK vectos from file */
                if ((frateKin = fopen(optarg,"rb")) == NULL) {
	            fprintf(stderr, "Error opening input rateK file: %s: %s\n",
		        optarg, strerror(errno));
                    exit(1);
                }
                fprintf(stderr, "rateKin: each record is %d bytes\n", (int)(K*sizeof(float)));
	    } else if(strcmp(long_options[option_index].name, "rateKnomeanin") == 0) {
                /* read mean removed rateK vectors from file */
                if ((frateKnomeanin = fopen(optarg,"rb")) == NULL) {
	            fprintf(stderr, "Error opening input rateK file: %s: %s\n",
		        optarg, strerror(errno));
                    exit(1);
                }
                fprintf(stderr, "rateKnomeanin: each record is %d bytes\n", (int)(K*sizeof(float)));
            } else if(strcmp(long_options[option_index].name, "rateK_mean_min") == 0) {
                rateK_mean_min = atof(optarg);
                fprintf(stderr, "rateK_mean_min: %f\n", rateK_mean_min);
            } else if(strcmp(long_options[option_index].name, "rateK_mean_max") == 0) {
                rateK_mean_max = atof(optarg);
                fprintf(stderr, "rateK_mean_max: %f\n", rateK_mean_max);
            } else if(strcmp(long_options[option_index].name, "rateK_mean_comp") == 0) {
                /* parse 4 values --rateK_mean_comp noisegate,thresh,slope,makeupgain */
                fprintf(stderr,"rateK mean compressor %s\n", optarg);
                char *comma, *p, val_str[80];
                p = optarg;
                float values[4];
                for(int n=0; n<4; n++) {
                    strcpy(val_str, p);
                    comma = strchr(val_str, ',');
                    if (comma) {
                        *comma = 0;
                        p = comma+1;
                    }
                    values[n] = atof(val_str);
                }
                rateKmeancomp.noise_gate = values[0];
                rateKmeancomp.thresh = values[1];
                rateKmeancomp.slope = values[2];
                rateKmeancomp.makeup_gain = values[3];
            } else if(strcmp(long_options[option_index].name, "dec") == 0) {

                decimate = atoi(optarg);
	        if ((decimate != 2) && (decimate != 3) && (decimate != 4)) {
	            fprintf(stderr, "Error in --dec, must be 2, 3, or 4\n");
	            exit(1);
	        }

                if (!phase0) {
                    fprintf(stderr, "needs --phase0 to resample phase when using --dec\n");
                    exit(1);
                }
                if (!lpc_model) {
                    fprintf(stderr, "needs --lpc [order] to resample amplitudes when using --dec\n");
                    exit(1);
                }

            } else if(strcmp(long_options[option_index].name, "hand_voicing") == 0) {
	        if ((fvoicing = fopen(optarg,"rt")) == NULL) {
	            fprintf(stderr, "Error opening voicing file: %s: %s.\n",
		        optarg, strerror(errno));
                    exit(1);
                }
            } else if(strcmp(long_options[option_index].name, "Woread") == 0) {
	        if ((fWo = fopen(optarg,"rb")) == NULL) {
	            fprintf(stderr, "Error opening float Wo file: %s: %s.\n",
		        optarg, strerror(errno));
                    exit(1);
                }
            } else if(strcmp(long_options[option_index].name, "amread") == 0) {
	        if ((fam = fopen(optarg,"rb")) == NULL) {
	            fprintf(stderr, "Error opening float Am file: %s: %s.\n",
		        optarg, strerror(errno));
                    exit(1);
                }
            } else if(strcmp(long_options[option_index].name, "hmread") == 0) {
	        if ((fhm = fopen(optarg,"rb")) == NULL) {
	            fprintf(stderr, "Error opening float Hm file: %s: %s.\n",
		        optarg, strerror(errno));
                    exit(1);
                }
            } else if(strcmp(long_options[option_index].name, "awread") == 0) {
	        if ((faw = fopen(optarg,"rb")) == NULL) {
	            fprintf(stderr, "Error opening float Aw file: %s: %s.\n",
                            optarg, strerror(errno));
                    exit(1);
                }
	    } else if(strcmp(long_options[option_index].name, "dump_pitch_e") == 0) {
	        if ((fjmv = fopen(optarg,"wt")) == NULL) {
	            fprintf(stderr, "Error opening pitch & energy dump file: %s: %s.\n",
		        optarg, strerror(errno));
                    exit(1);
                }
	    } else if(strcmp(long_options[option_index].name, "gainoutlin") == 0) {
		gainoutlin = atof(optarg);
	    } else if(strcmp(long_options[option_index].name, "comp") == 0) {
		comp_limit_samples = atof(optarg);
                comp_en = 1;
	    } else if(strcmp(long_options[option_index].name, "comp_gain") == 0) {
		comp_gain_dB = atof(optarg);
                comp_en = 1;
        } else if(strcmp(long_options[option_index].name, "comp_dr") == 0) {
		float comp_dr_dB = atof(optarg);
                comp_dr_lin = pow(10,-comp_dr_dB/20);
	    } else if(strcmp(long_options[option_index].name, "framelength_s") == 0) {
		framelength_s = atof(optarg);
	    } else if(strcmp(long_options[option_index].name, "pahw") == 0) {

                /* set up a bunch of arguments instead of having to enter them on cmd line every time */

                phase0 = postfilt = amread = hmread = Woread = 1;
                char file_name[MAX_STR];
                sprintf(file_name, "%s_am.out", optarg);
                fprintf(stderr, "reading %s", file_name);
	        if ((fam = fopen(file_name,"rb")) == NULL) {
	            fprintf(stderr, "Error opening float Am file: %s: %s.\n",
		        file_name, strerror(errno));
                    exit(1);
                }
                sprintf(file_name, "%s_hm.out", optarg);
                fprintf(stderr, " %s", file_name);
	        if ((fhm = fopen(file_name,"rb")) == NULL) {
	            fprintf(stderr, "Error opening float Hm file: %s: %s.\n",
		        file_name, strerror(errno));
                    exit(1);
                }
                sprintf(file_name, "%s_Wo.out", optarg);
                fprintf(stderr, " %s\n", file_name);
 	        if ((fWo = fopen(file_name,"rb")) == NULL) {
	            fprintf(stderr, "Error opening float Wo file: %s: %s.\n",
		        file_name, strerror(errno));
                    exit(1);
                }
	    } else if(strcmp(long_options[option_index].name, "lspEWov") == 0) {
                /* feature file for deep learning experiments */
                lpc_model = 1; phase0 = 1;
	        if ((flspEWov = fopen(optarg,"wb")) == NULL) {
	            fprintf(stderr, "Error opening lspEWov float file: %s: %s\n",
		        optarg, strerror(errno));
                    exit(1);
                }
	    } else if(strcmp(long_options[option_index].name, "rateKWov") == 0) {
                /* feature file for deep learning experiments */
                rateK = 1; newamp1vq = 1;
	        if ((frateKWov = fopen(optarg,"wb")) == NULL) {
	            fprintf(stderr, "Error opening rateKWov float file: %s: %s\n",
		        optarg, strerror(errno));
                    exit(1);
                }
	    } else if(strcmp(long_options[option_index].name, "ten_ms_centre") == 0) {
                /* dump 10ms of audio centred on analysis frame to check time alignment with
                   16 kHz source audio */
                ten_ms_centre = 1;
	        if ((ften_ms_centre = fopen(optarg,"wb")) == NULL) {
	            fprintf(stderr, "Error opening ten_ms_centre short file: %s: %s\n",
		        optarg, strerror(errno));
                    exit(1);
                }
	    } else if(strcmp(long_options[option_index].name, "modelout") == 0) {
                /* write model records to file or stdout */
                modelout = 1;
                if (strcmp(optarg, "-") == 0) fmodelout = stdout;
	        else if ((fmodelout = fopen(optarg,"wb")) == NULL) {
	            fprintf(stderr, "Error opening modelout file: %s: %s\n",
		        optarg, strerror(errno));
                    exit(1);
                }
                fprintf(stderr, "modelout: each model record is %d bytes\n", (int)sizeof(MODEL));
	    } else if(strcmp(long_options[option_index].name, "modelin") == 0) {
                /* read model records from file or stdin */
                modelin = 1;
                if (strcmp(optarg, "-") == 0) fmodelin = stdin;
	        else if ((fmodelin = fopen(optarg,"rb")) == NULL) {
	            fprintf(stderr, "Error opening modelin file: %s: %s\n",
		        optarg, strerror(errno));
                    exit(1);
                }
                fprintf(stderr, "modelin: each model record is %d bytes\n", (int)sizeof(MODEL));
            } else if(strcmp(long_options[option_index].name, "rate") == 0) {
                if(strcmp(optarg,"3200") == 0) {
	            lpc_model = 1;
		    scalar_quant_Wo_e = 1;
	            lspd = 1;
	            phase0 = 1;
	            postfilt = 1;
	            decimate = 1;
		    lpcpf = 1;
               } else if(strcmp(optarg,"2400") == 0) {
	            lpc_model = 1;
		    vector_quant_Wo_e = 1;
	            lsp = 1;
	            phase0 = 1;
	            postfilt = 1;
	            decimate = 2;
		    lpcpf = 1;
               } else if(strcmp(optarg,"1400") == 0) {
	            lpc_model = 1;
		    vector_quant_Wo_e = 1;
	            lsp = 1;
	            phase0 = 1;
	            postfilt = 1;
	            decimate = 4;
 		    lpcpf = 1;
               } else if(strcmp(optarg,"1300") == 0) {
	            lpc_model = 1;
		    scalar_quant_Wo_e = 1;
	            lsp = 1;
	            phase0 = 1;
	            postfilt = 1;
	            decimate = 4;
 		    lpcpf = 1;
               } else if(strcmp(optarg,"1200") == 0) {
	            lpc_model = 1;
		    scalar_quant_Wo_e = 1;
	            lspjmv = 1;
	            phase0 = 1;
	            postfilt = 1;
	            decimate = 4;
 		    lpcpf = 1;
                } else {
                    fprintf(stderr, "Error: invalid output rate (3200|2400|1400|1200) %s\n", optarg);
                    exit(1);
                }
            }
            break;

         case 'h':
             print_help(long_options, num_opts, argv);
             break;

         case 'o':
	     if (strcmp(optarg, "-") == 0) fout = stdout;
	     else if ((fout = fopen(optarg,"wb")) == NULL) {
	        fprintf(stderr, "Error opening output speech file: %s: %s.\n",
		    optarg, strerror(errno));
	        exit(1);
	     }
	     strcpy(out_file,optarg);
	     break;

         default:
            /* This will never be reached */
            break;
        }
    }

    /* Input file */

    FILE *fin;		/* input speech file                     */
    if (strcmp(argv[optind], "-")  == 0) fin = stdin;
    else if ((fin = fopen(argv[optind],"rb")) == NULL) {
	fprintf(stderr, "Error opening input speech file: %s: %s.\n",
		argv[optind], strerror(errno));
	exit(1);
    }

    C2CONST c2const = c2const_create(Fs, framelength_s);
    int   n_samp = c2const.n_samp;
    int   m_pitch = c2const.m_pitch;

    short buf[N_SAMP];	/* input/output buffer                   */
    float buf_float[N_SAMP];
    float Sn[m_pitch];	/* float input speech samples            */
    float Sn_pre[m_pitch];	/* pre-emphasised input speech samples   */
    COMP  Sw[FFT_ENC];	/* DFT of Sn[]                           */
    codec2_fft_cfg  fft_fwd_cfg;
    codec2_fftr_cfg  fftr_fwd_cfg;
    codec2_fft_cfg  fft_inv_cfg;
    float w[m_pitch];	        /* time domain hamming window            */
    float W[FFT_ENC];	/* DFT of w[]                            */
    MODEL model;
    float Pn[2*N_SAMP];	/* trapezoidal synthesis window          */
    float Sn_[2*N_SAMP];	/* synthesised speech */
    int   i,m;		/* loop variable                         */
    int   frames;
    float prev_f0;
    float pitch;
    float snr;
    float sum_snr;

    float pre_mem = 0.0, de_mem = 0.0;
    float ak[1+order];
    // COMP  Sw_[FFT_ENC];
    // COMP  Ew[FFT_ENC];

    float ex_phase[MAX_AMP+1];

    float bg_est = 0.0;

    MODEL prev_model;
    float lsps[order];
    float e, prev_e;
    int   lsp_indexes[order];
    float lsps_[order];
    float Woe_[2];

    float lsps_dec[4][order], e_dec[4], weight, weight_inc, ak_dec[4][order];
    MODEL model_dec[4], prev_model_dec;
    float prev_lsps_dec[order], prev_e_dec;

    void *nlp_states;
    float bpf_buf[BPF_N+N_SAMP];

    COMP Aw[FFT_ENC];
    COMP H[MAX_AMP];

    float sd_sum = 0.0; int sd_frames = 0;
    
    for(i=0; i<m_pitch; i++) {
	Sn[i] = 1.0;
	Sn_pre[i] = 1.0;
    }
    for(i=0; i<2*N_SAMP; i++)
	Sn_[i] = 0;

    prev_f0 = 1/P_MAX_S;

    prev_model.Wo = c2const.Wo_max;
    prev_model.L = floor(PI/prev_model.Wo);
    for(i=1; i<=prev_model.L; i++) {
	prev_model.A[i] = 0.0;
	prev_model.phi[i] = 0.0;
    }
    for(i=1; i<=MAX_AMP; i++) {
	//ex_phase[i] = (PI/3)*(float)rand()/RAND_MAX;
	ex_phase[i] = 0.0;
    }
    e = prev_e = 1;

    nlp_states = nlp_create(&c2const);

    ex_phase[0] = 0;
    Woe_[0] = Woe_[1] = 1.0;

    /* Initialise ------------------------------------------------------------*/

    fft_fwd_cfg = codec2_fft_alloc(FFT_ENC, 0, NULL, NULL);   /* fwd FFT,used in several places   */
    fftr_fwd_cfg = codec2_fftr_alloc(FFT_ENC, 0, NULL, NULL); /* fwd FFT,used in several places   */
    fft_inv_cfg = codec2_fft_alloc(FFT_DEC, 1, NULL, NULL);   /* inverse FFT, used just for synth */
    codec2_fft_cfg phase_fft_fwd_cfg = codec2_fft_alloc(NEWAMP1_PHASE_NFFT, 0, NULL, NULL);
    codec2_fft_cfg phase_fft_inv_cfg = codec2_fft_alloc(NEWAMP1_PHASE_NFFT, 1, NULL, NULL);

    make_analysis_window(&c2const, fft_fwd_cfg, w, W);
    make_synthesis_window(&c2const, Pn);

    assert(BPF_N == BPFB_N); assert(BPF_N == HPF200_N);
    for(i=0; i<BPF_N; i++)
        bpf_buf[i] = 0.0;

    for(i=0; i<LPC_ORD; i++) {
        prev_lsps_dec[i] = i*PI/(LPC_ORD+1);
    }
    prev_e_dec = 1;
    for(m=1; m<=MAX_AMP; m++)
	prev_model_dec.A[m] = 0.0;
    prev_model_dec.Wo = c2const.Wo_min;
    prev_model_dec.L = PI/prev_model_dec.Wo;
    prev_model_dec.voiced = 0;

    /* mel resampling experiments */

    float rate_K_sample_freqs_kHz[K]; float se = 0.0; int nse = 0;
    if (rateK) {
	mel_sample_freqs_kHz(rate_K_sample_freqs_kHz, K, ftomel(200.0), ftomel(3700.0) );
    }
    float rate_K_vec_delay[rate_K_dec+1][K];
    float rate_K_vec_delay_[rate_K_dec+1][K];
    MODEL rate_K_model_delay[rate_K_dec+1];
    for (int d=0; d<=rate_K_dec; d++) {
        for(int k=0; k<K; k++) {
            rate_K_vec_delay[d][k] = 0;
            rate_K_vec_delay_[d][k] = 0;
        }
        for(m=1; m<=MAX_AMP; m++)
            rate_K_model_delay[d].A[m] = 0.0;
        rate_K_model_delay[d].Wo = c2const.Wo_min;
        rate_K_model_delay[d].L = M_PI/prev_model_dec.Wo;
        rate_K_model_delay[d].voiced = 0;
    }
    float eq[K];
    for(int k=0; k<K; k++) eq[k] = 0;

    /*----------------------------------------------------------------* \

                            Main Loop

    \*----------------------------------------------------------------*/

    frames = 0;
    sum_snr = 0;
    while(fread(buf,sizeof(short),N_SAMP,fin)) {
	frames++;

	for(i=0; i<N_SAMP; i++)
	    buf_float[i] = buf[i];

	/* optionally filter input speech */

        if (prede) {
           pre_emp(Sn_pre, buf_float, &pre_mem, N_SAMP);
           for(i=0; i<N_SAMP; i++)
                buf_float[i] = Sn_pre[i];
        }

        if (bpf_en || bpfb_en || hpf200_en) {
            /* filter input speech to create buf_float_bpf[], this is fed to the
               LPC modelling.  Unfiltered speech in in buf_float[], which is
               delayed to match that of the BPF */

            /* BPF speech */

            for(i=0; i<BPF_N; i++)
                bpf_buf[i] =  bpf_buf[N_SAMP+i];
            for(i=0; i<N_SAMP; i++)
                bpf_buf[BPF_N+i] = buf_float[i];
            if (bpf_en)
                inverse_filter(&bpf_buf[BPF_N], bpf, N_SAMP, buf_float, BPF_N);
            if (bpfb_en)
                inverse_filter(&bpf_buf[BPF_N], bpfb, N_SAMP, buf_float, BPFB_N);
            if (hpf200_en)
                inverse_filter(&bpf_buf[BPF_N], hpf200, N_SAMP, buf_float, HPF200_N);
        }

        /* shift buffer of input samples, and insert new samples */

	for(i=0; i<m_pitch-N_SAMP; i++) {
	    Sn[i] = Sn[i+N_SAMP];
	}
	for(i=0; i<N_SAMP; i++) {
	    Sn[i+m_pitch-N_SAMP] = buf_float[i];
        }

	/*------------------------------------------------------------*\

                      Estimate Sinusoidal Model Parameters

	\*------------------------------------------------------------*/

        nlp(nlp_states, Sn, N_SAMP, &pitch, Sw, W, &prev_f0);
	model.Wo = TWO_PI/pitch;

        dft_speech(&c2const, fft_fwd_cfg, Sw, Sn, w);
	two_stage_pitch_refinement(&c2const, &model, Sw);
	estimate_amplitudes(&model, Sw, W, 1);

        /* estimate voicing - do this all the time so model.voicing
	 * is set, useful for machine learning work */
	snr = est_voicing_mbe(&c2const, &model, Sw, W);
        
        /*------------------------------------------------------------*\

                            Hilbert Compressor

	\*------------------------------------------------------------*/

        if (comp_en) {
            float Acomp[MAX_AMP+1];
            float s_max = 0.0;
            float comp_gain_lin = pow(10.0, comp_gain_dB/20.0);
            
            // apply compressor gain
            for(int m=1; m<=model.L; m++) {
                Acomp[m] = 2.0*comp_gain_lin*model.A[m];
                s_max += Acomp[m];
            }
            float g = 1.0;
            if (s_max > comp_limit_samples)
                g = comp_limit_samples/s_max;
            /*fprintf(stderr, "f: %d s_max_low: %6.0f s_max_high: %6.0f gains: %3.2f %3.2f\n", 
                            frames, s_max_low, s_max_high, g_low, g_high); */
            float mx = 0.0;
            for(int m=1; m<=model.L; m++) {
                model.A[m] = g*Acomp[m];
                if (model.A[m] > mx) mx = model.A[m];
            }
            
            // limit dynamic range
            for(int m=1; m<=model.L; m++)
              if (model.A[m] < mx*comp_dr_lin) model.A[m] = mx*comp_dr_lin;
                 
        }
        
        #ifdef DUMP
        dump_Sn(m_pitch, Sn); dump_Sw(Sw); dump_model(&model);
        #endif

	/* speech centred on analysis frame for Deep Learning work */

	if (ten_ms_centre) {
	    int n_10_ms = Fs*0.01;
	    int n_5_ms = Fs*0.005;
	    short buf[n_10_ms];
	    for(i=0; i<n_10_ms; i++) {
		buf[i] = Sn[m_pitch/2-n_5_ms+i];
	    }
	    fwrite(buf, n_10_ms, sizeof(short), ften_ms_centre);
	}

	if (hi) {
	    int m;
	    for(m=1; m<model.L/2; m++)
		model.A[m] = 0.0;
	    for(m=3*model.L/4; m<=model.L; m++)
		model.A[m] = 0.0;
	}
        
	/*------------------------------------------------------------*\

                            Zero-phase modelling

	\*------------------------------------------------------------*/

	if (phase0) {
            #ifdef DUMP
	    dump_phase(&model.phi[0], model.L);
            #endif

	    if (dump_pitch_e)
		fprintf(fjmv, "%f %f %d ", model.Wo, snr, model.voiced);

            #ifdef DUMP
	    dump_snr(snr);
            #endif

	    /* just to make sure we are not cheating - kill all phases */

	    for(i=0; i<=MAX_AMP; i++)
	    	model.phi[i] = 0;

	    if (hand_voicing) {
		int ret = fscanf(fvoicing,"%d\n",&model.voiced);
                assert(ret == 1);
	    }
	}

	/*------------------------------------------------------------*\

	        LPC model amplitudes and LSP quantisation

	\*------------------------------------------------------------*/

	if (lpc_model) {
            float ak_[LPC_ORD+1];

            e = speech_to_uq_lsps(lsps, ak, Sn, w, m_pitch, order);
            for(i=0; i<order; i++)
                lsps_[i] = lsps[i];

            #ifdef DUMP
	    dump_ak(ak, order);
            dump_E(e);
            #endif

	    if (dump_pitch_e)
		fprintf(fjmv, "%f\n", e);

            #ifdef DUMP
            dump_lsp(lsps);
            #endif

	    /* various LSP quantisation schemes */

	    if (lsp) {
		encode_lsps_scalar(lsp_indexes, lsps, LPC_ORD);
		decode_lsps_scalar(lsps_, lsp_indexes, LPC_ORD);
		bw_expand_lsps(lsps_, LPC_ORD, 50.0, 100.0);
		lsp_to_lpc(lsps_, ak_, LPC_ORD);
	    }

	    if (lspd) {
		encode_lspds_scalar(lsp_indexes, lsps, LPC_ORD);
		decode_lspds_scalar(lsps_, lsp_indexes, LPC_ORD);
		lsp_to_lpc(lsps_, ak_, LPC_ORD);
	    }

	    if (lspjmv) {
		/* Jean-Marc's multi-stage, split VQ */
		lspjmv_quantise(lsps, lsps_, LPC_ORD);
		{
		    float lsps_bw[LPC_ORD];
		    memcpy(lsps_bw, lsps_, sizeof(float)*order);
		    bw_expand_lsps(lsps_bw, LPC_ORD, 50.0, 100.0);
		    lsp_to_lpc(lsps_bw, ak_, LPC_ORD);
		}
	    }

            if (lsp || lspd || lspjmv) {
                sd_sum += spectral_dist(ak, ak_, LPC_ORD, fft_fwd_cfg, FFT_ENC);
                sd_frames ++;
            }

            memcpy(ak, ak_, (LPC_ORD+1)*sizeof(float));

	    if (scalar_quant_Wo_e) {
		e = decode_energy(encode_energy(e, E_BITS), E_BITS);
                model.Wo = decode_Wo(&c2const, encode_Wo(&c2const, model.Wo, WO_BITS), WO_BITS);
		model.L  = PI/model.Wo; /* if we quantise Wo re-compute L */
	    }

	    if (scalar_quant_Wo_e_low) {
                int ind;
		e = decode_energy(ind = encode_energy(e, 3), 3);
                model.Wo = decode_log_Wo(&c2const, encode_log_Wo(&c2const, model.Wo, 5), 5);
		model.L  = PI/model.Wo; /* if we quantise Wo re-compute L */
	    }

	    if (vector_quant_Wo_e) {
		/* JVM's experimental joint Wo & LPC energy quantiser */
		quantise_WoE(&c2const, &model, &e, Woe_);
	    }

	}

        if (amread) {
            int ret = fread(model.A, sizeof(float), MAX_AMP, fam);
            assert(ret == MAX_AMP);
        }

        if (Woread) {
            int ret = fread(&model.Wo, sizeof(float), 1, fWo);
            model.L = floor(PI/model.Wo);
            assert(ret == 1);
        }

        /* dump features for Deep learning, placed here so we can get quantised features */

        if (lspEWov) {
            /* order LSPs - energy - Wo - voicing flag - order LPCs */
            if (lsp)
                fwrite(lsps_, order, sizeof(float), flspEWov);
            else
                fwrite(lsps, order, sizeof(float), flspEWov);

            fwrite(&e, 1, sizeof(float), flspEWov);
            fwrite(&model.Wo, 1, sizeof(float), flspEWov);
            float voiced_float = model.voiced;
            fwrite(&voiced_float, 1, sizeof(float), flspEWov);
            fwrite(&ak[1], order, sizeof(float), flspEWov);
        }

	/*------------------------------------------------------------*\

	   newamp1/rateK, as used in 700C.  We resample the rate L {Am}
           samples to rate K.  L is dependant on the pitch Wo, so is
           time varying. K is constant, but sampled at non-linear 
           points on the frequency axis.

	\*------------------------------------------------------------*/

        if (rateK) {
            float rate_K_vec[K];
            resample_const_rate_f(&c2const, &model, rate_K_vec, rate_K_sample_freqs_kHz, K);

	    if (frateKin != NULL) {
		assert(fread(rate_K_vec, sizeof(float), K, frateKin) == K);
	    }
        
            /* remove mean (note EQ and post filter also work on mean removed vector) */            
	    float sum = 0.0;
	    for(int k=K_st; k<=K_en; k++)
	      sum += rate_K_vec[k];
	    float mean = sum/(K_en-K_st+1);
	    float rate_K_vec_no_mean[K]; float rate_K_vec_no_mean_[K];
	    for(int k=0; k<K; k++)
	      rate_K_vec_no_mean[k] = rate_K_vec[k] - mean;

            /* mean compressor */
            if (rateK_mean_comp_en) {
                if (mean < rateKmeancomp.noise_gate)
                    mean = 0;
                if (mean > rateKmeancomp.thresh)
                    mean = rateKmeancomp.thresh + (mean-rateKmeancomp.thresh)*rateKmeancomp.slope;
                mean += rateKmeancomp.makeup_gain;
            }
            
            /* clip (limit) mean */
            if (mean < rateK_mean_min) mean = rateK_mean_min;
	    if (mean > rateK_mean_max) mean = rateK_mean_max;

            if (frateKmean != NULL) assert(fwrite(&mean, sizeof(float), 1, frateKmean) == 1);
            
            // dump output for VQ training
	    if (frateK != NULL) {
                float rate_K_vec[K];
                for(int k=0; k<K; k++)
                    rate_K_vec[k] = rate_K_vec_no_mean[k] + mean;
		assert(fwrite(rate_K_vec, sizeof(float), K, frateK) == K);
            }
	    if (frateKnomean != NULL) {
                assert(fwrite(rate_K_vec_no_mean, sizeof(float), K, frateKnomean) == K);
            }
	    if (frateKnomeanin != NULL) {
		assert(fread(rate_K_vec_no_mean, sizeof(float), K, frateKnomeanin) == K);
	    }

            /* newamp1 VQ */
            float rate_K_vec_[K];
            if (newamp1vq) {
		newamp1_eq(rate_K_vec_no_mean, eq, K, 1);

                /* two stage VQ */
                int indexes[2];
                rate_K_mbest_encode(indexes, rate_K_vec_no_mean, rate_K_vec_no_mean_, K, NEWAMP1_VQ_MBEST_DEPTH);

                /* running sum of squared error for variance calculation */
                for(int k=0; k<K; k++)
                    se += pow(rate_K_vec_no_mean[k]-rate_K_vec_no_mean_[k],2.0);
                nse += K;
            }
            else {
                for(int k=0; k<K; k++)
                    rate_K_vec_no_mean_[k] = rate_K_vec_no_mean[k];
            }

	    if (postfilter_newamp1_en)
	        post_filter_newamp1(rate_K_vec_no_mean_, rate_K_sample_freqs_kHz, K, 1.5);
	    for(int k=0; k<K; k++)
                rate_K_vec_[k] = rate_K_vec_no_mean_[k] +  mean;

            /* decimation/interpolation of rateK vectors */
            if (rate_K_dec) {
                // update delay lines
                for(int d=0; d<rate_K_dec; d++) {
                    rate_K_model_delay[d] = rate_K_model_delay[d+1];
                    memcpy(&rate_K_vec_delay[d][0], &rate_K_vec_delay[d+1][0], sizeof(float)*K);
                }
                rate_K_model_delay[rate_K_dec] = model;
                memcpy(&rate_K_vec_delay[rate_K_dec][0], rate_K_vec_, sizeof(float)*K);

                if ((frames % rate_K_dec) == 0) {
                    // every rate_K_dec frames, calculate interpolated output values
                    if (perframe) {
                        // calculate interpolation coeff c for each frame
                        float *A = &rate_K_vec_delay[0][0];
                        float *B = &rate_K_vec_delay[rate_K_dec][0];
                        for(int d=0; d<=rate_K_dec; d++) {
                            float *T = &rate_K_vec_delay[d][0];
                            float num = 0.0, den = 0.0;
                            for(int k=0; k<K; k++) {
                                num += (B[k]-T[k])*(A[k]-B[k]);
                                den += (A[k]-B[k])*(A[k]-B[k]);
                            }
                            float c = -num/den;
                            for(int k=0; k<K; k++)
                                rate_K_vec_delay_[d][k] = c*A[k] + (1.0-c)*B[k];
                        }
                    }
                    else {
                        // use linear interpolation
                        float c=0.0, inc = 1.0/rate_K_dec;
                        for(int d=0; d<=rate_K_dec; d++) {
                            for(int k=0; k<K; k++)
                                rate_K_vec_delay_[d][k] = (1.0-c)*rate_K_vec_delay[0][k] + c*rate_K_vec_delay[rate_K_dec][k];
                            c += inc;
                        }
                    }
                } else {
                    // otherwise just shift out frames we have already interpolated
                    for(int d=0; d<rate_K_dec; d++) {
                        memcpy(&rate_K_vec_delay_[d][0], &rate_K_vec_delay_[d+1][0], sizeof(float)*K);
                    }
                }

                // output from delay line
                model = rate_K_model_delay[0];
                for(int k=0; k<K; k++)
                    rate_K_vec_[k] = rate_K_vec_delay_[0][k];
            }
	    
            float A_backup[MAX_AMP];
            for(int m=1; m<=model.L; m++)
                A_backup[m] = model.A[m];
            
            resample_rate_L(&c2const, &model, rate_K_vec_, rate_K_sample_freqs_kHz, K);
            
            if (rateKlf) {
                for(int m=model.L/4; m<=model.L; m++)
                    model.A[m] = A_backup[m];
            }
            
        } // rateK

	/*------------------------------------------------------------*\

          Synthesise and optional decimation to 20 or 40ms frame rate

	\*------------------------------------------------------------*/

        /*
           if decimate == 2, we interpolate frame n from frame n-1 and n+1
           if decimate == 4, we interpolate frames n, n+1, n+2, from frames n-1 and n+3

           This is meant to give identical results to the implementations of various modes
           in codec2.c
        */

        /* delay line to keep frame by frame voicing decisions */

        for(i=0; i<decimate-1; i++)
            model_dec[i] = model_dec[i+1];
        model_dec[decimate-1] = model;

        if ((frames % decimate) == 0) {
            for(i=0; i<order; i++)
                lsps_dec[decimate-1][i] = lsps_[i];
            e_dec[decimate-1] = e;
            model_dec[decimate-1] = model;

            /* interpolate the model parameters */

            weight_inc = 1.0/decimate;
            for(i=0, weight=weight_inc; i<decimate-1; i++, weight += weight_inc) {
                //model_dec[i].voiced = model_dec[decimate-1].voiced;
                interpolate_lsp_ver2(&lsps_dec[i][0], prev_lsps_dec, &lsps_dec[decimate-1][0], weight, order);
                interp_Wo2(&model_dec[i], &prev_model_dec, &model_dec[decimate-1], weight, c2const.Wo_min);
                e_dec[i] = interp_energy2(prev_e_dec, e_dec[decimate-1],weight);
            }

            /* then recover spectral amplitudes and synthesise */

            for(i=0; i<decimate; i++) {
                if (lpc_model) {
                    lsp_to_lpc(&lsps_dec[i][0], &ak_dec[i][0], order);
                    aks_to_M2(fftr_fwd_cfg, &ak_dec[i][0], order, &model_dec[i], e_dec[i],
                              &snr, 0, simlpcpf, lpcpf, 1, LPCPF_BETA, LPCPF_GAMMA, Aw);
                    apply_lpc_correction(&model_dec[i]);
                    sum_snr += snr;
                    #ifdef DUMP
                    dump_lsp_(&lsps_dec[i][0]);
                    dump_ak_(&ak_dec[i][0], order);
                    #endif
                }

                if (modelin) {
                    int nrec;
                    nrec = fread(&model_dec[i],sizeof(MODEL),1,fmodelin);
                    if (nrec != 1) {
			fprintf(stderr, "Warning - error reading model in record in frame %d - do you have enough records in file?\n", frames);
		    }
                }

                if (phase0) {
                    /* optionally read in Aw, replacing values generated using LPC */
                    if (awread) {
                        int ret = fread(Aw, sizeof(COMP), FFT_ENC, faw);
                        assert(ret == FFT_ENC);
                    }

                    /* optionally read in Hm directly, bypassing sampling of Aw[] */
                    if (hmread) {
                        int ret = fread(H, sizeof(COMP), MAX_AMP, fhm);
                        assert(ret == MAX_AMP);
                    } else {
                        determine_phase(&c2const, H, &model_dec[i], NEWAMP1_PHASE_NFFT, phase_fft_fwd_cfg, phase_fft_inv_cfg);
                    }
                    phase_synth_zero_order(n_samp, &model_dec[i], ex_phase, H);
                    #ifdef DUMP
                    dump_H(H, model.L);
                    dump_phase_(&model_dec[i].phi[0], model.L);
                    #endif
                }

                if (postfilt)
                    postfilter(&model_dec[i], &bg_est);
                synth_one_frame(n_samp, fft_inv_cfg, buf, &model_dec[i], Sn_, Pn, prede, &de_mem, gainoutlin);
                #ifdef DUMP
                dump_quantised_model(&model_dec[i]);
                #endif
                if (fout != NULL)
                    fwrite(buf,sizeof(short),N_SAMP,fout);
                if (modelout)
                    fwrite(&model_dec[i],sizeof(MODEL),1,fmodelout);
            }

            /* update memories for next frame ----------------------------*/

            prev_model_dec = model_dec[decimate-1];
            prev_e_dec = e_dec[decimate-1];
            for(i=0; i<LPC_ORD; i++)
                prev_lsps_dec[i] = lsps_dec[decimate-1][i];
       }

    }

    /*----------------------------------------------------------------*\

                            End Main Loop

    \*----------------------------------------------------------------*/

    fclose(fin);

    if (fout != NULL)
	fclose(fout);

    if (lpc_model) {
    	fprintf(stderr, "LPC->{Am} SNR av: %5.2f dB over %d frames\n", sum_snr/frames, frames);
        if (lsp || lspd || lspjmv)
            fprintf(stderr, "LSP quantiser SD: %5.2f dB*dB over %d frames\n", sd_sum/sd_frames, sd_frames);     
    }
    if (newamp1vq) {
    	fprintf(stderr, "var: %3.2f dB*dB\n", se/nse);
    }
    #ifdef DUMP
    if (dump)
	dump_off();
    #endif

    if (hand_voicing)
	fclose(fvoicing);

    nlp_destroy(nlp_states);

    if (fam     != NULL) fclose(fam);
    if (fWo     != NULL) fclose(fWo);
    if (faw     != NULL) fclose(faw);
    if (fhm     != NULL) fclose(fhm);
    if (fjmv    != NULL) fclose(fjmv);
    if (flspEWov != NULL) fclose(flspEWov);
    if (fphasenn != NULL) fclose(fphasenn);
    if (frateK != NULL) fclose(frateK);
    if (frateKnomean != NULL) fclose(frateKnomean);
    if (frateKmean != NULL) fclose(frateKmean);
    if (frateKin != NULL) fclose(frateKin);
    if (frateKnomeanin != NULL) fclose(frateKnomeanin);
    if (ften_ms_centre != NULL) fclose(ften_ms_centre);
    if (fmodelout != NULL) fclose(fmodelout);
    if (frateKWov != NULL) fclose(frateKWov);

    return 0;
}

void synth_one_frame(int n_samp, codec2_fft_cfg fft_inv_cfg, short buf[], MODEL *model, float Sn_[],
                     float Pn[], int prede, float *de_mem, float gain)
{
    int     i;

    synthesise(n_samp, fft_inv_cfg, Sn_, model, Pn, 1);
    if (prede)
        de_emp(Sn_, Sn_, de_mem, n_samp);

    for(i=0; i<n_samp; i++) {
	Sn_[i] *= gain;
	if (Sn_[i] > 32767.0)
	    buf[i] = 32767;
	else if (Sn_[i] < -32767.0)
	    buf[i] = -32767;
	else
	    buf[i] = Sn_[i];
    }

}

void print_help(const struct option* long_options, int num_opts, char* argv[])
{
	int i;
	char *option_parameters;

	fprintf(stderr, "\nCodec2 - low bit rate speech codec - Simulation Program\n"
		"\thttp://rowetel.com/codec2.html\n\n"
		"usage: %s [OPTIONS] <InputFile>\n\n"
                "Options:\n"
                "\t-o <OutputFile>\n", argv[0]);
        for(i=0; i<num_opts-1; i++) {
		if(long_options[i].has_arg == no_argument) {
			option_parameters="";
		} else if (strcmp("lpc", long_options[i].name) == 0) {
			option_parameters = " <Order>";
		} else if (strcmp("dec", long_options[i].name) == 0) {
			option_parameters = " <2|4>";
		} else if (strcmp("hand_voicing", long_options[i].name) == 0) {
			option_parameters = " <VoicingFile>";
		} else if (strcmp("dump_pitch_e", long_options[i].name) == 0) {
			option_parameters = " <Dump File>";
		} else if (strcmp("rate", long_options[i].name) == 0) {
			option_parameters = " <3200|2400|1400|1300|1200>";
		} else if (strcmp("dump", long_options[i].name) == 0) {
			option_parameters = " <DumpFilePrefix>";
		} else {
			option_parameters = " <UNDOCUMENTED parameter>";
		}
		fprintf(stderr, "\t--%s%s\n", long_options[i].name, option_parameters);
	}

	exit(1);
}
