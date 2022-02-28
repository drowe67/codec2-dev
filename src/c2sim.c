/*---------------------------------------------------------------------------*\

  FILE........: c2sim.c
  AUTHOR......: David Rowe
  DATE CREATED: 20/8/2010

  Codec2 simulation.  Combines encoder and decoder and allows
  switching in and out various algorithms and quantisation steps. Used
  for algorithm development.

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

void synth_one_frame(int n_samp, codec2_fftr_cfg fftr_inv_cfg, short buf[], MODEL *model, float Sn_[], float Pn[], int prede, float *de_mem, float gain);
void print_help(const struct option *long_options, int num_opts, char* argv[]);

#define N_SAMP n_samp  /* quick fix for run time sample rate selection */

/*---------------------------------------------------------------------------*\

				MAIN

\*---------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{

    int Fs = 8000;
    int set_fs;

    int   lpc_model = 0, order = LPC_ORD;
    int   lsp = 0, lspd = 0, lspvq = 0;
    int   lspjmv = 0;
    int   prede = 0;
    int   postfilt;
    int   hand_voicing = 0, hi = 0, simlpcpf = 0, modelin=0, modelout=0;
    int   lpcpf = 0;
    FILE *fvoicing = 0;
    int   dec;
    int   decimate = 1;
    int   amread, Woread, pahw;
    int   awread;
    int   hmread;
    int   phase0 = 0;
    int   scalar_quant_Wo_e = 0;
    int   scalar_quant_Wo_e_low = 0;
    int   vector_quant_Wo_e = 0;
    int   dump_pitch_e = 0;
    float gain = 1.0;
    int   bpf_en = 0;
    int   bpfb_en = 0;
    int   cnw = 0, custom_nw = 0;
    int   dmp = 0;
    int   framecounter = 0;
    int   minmax = 0;
    int   statt = 1;
    int	  sqd = 0;
    int   sqdd = 255;
    int   sqddopt = 0;
    int   lspvt = 0;
    int   opti = 0;
    int   clr = 0;
    long  vtl = 0L;
    long  vtl1 = 0L;
    long  vtl2 = 0L;
    long  vtlb = 0L;
    long  sqerr = 0L, sqermn = 0L, bstix = 0L;
    long  sqerr1 = 0L, sqermn1 = 0L, bstix1 = 0L;
    long  sqerr2 = 0L, sqermn2 = 0L, bstix2 = 0L;
    int   maxoc = 0;
    int   maxoc1 = 0;
    int   maxoc2 = 0;
    long  frmcntl = 0L;
    int   vq10 = 0;
    int   cdb = 0, cdbsz = 0;
    int   cdb1 = 0, cdbsz1 = 0;
    int   cdb2 = 0, cdbsz2 = 0;
    int   wgt = 0;
    int   dfr = 0;
    int   lbg = 0;
    float lbgd0 = 0.0;
    int   n;
    int   dbl = 0;
    int   prev_vo = 0;
    int   trn = 0;
    float v_thresh = 6.0;
    int   thrs = 0;
    FILE *fam = NULL, *fWo = NULL;
    FILE *faw = NULL;
    FILE *fhm = NULL;
    FILE *fjmv = NULL;
    FILE *flspEWov = NULL;
    FILE *ften_ms_centre = NULL;
    FILE *fmodelout = NULL;
    FILE *fmodelin = NULL;
    FILE *fdmp = NULL;
    FILE *fdmp1 = NULL;
    FILE *fdmp2 = NULL;
    FILE *faux = NULL;
    FILE *fwght = NULL;
    #ifdef DUMP
    int   dump;
    #endif
    char  out_file[MAX_STR];
    FILE *fout = NULL;	/* output speech file */
    int   rateK = 0, newamp1vq = 0, rate_K_dec = 0, perframe=0;
    int   bands = 0, bands_lower_en;
    float bands_lower = -1E32;
    int   K = 20;
    float framelength_s = N_S;
    int   lspEWov = 0, rateKWov = 0, first = 0;
    FILE  *frateKWov = NULL;
    int   ten_ms_centre = 0;
    FILE  *fphasenn = NULL;
    FILE  *frateK = NULL;
    FILE  *frateKin = NULL;
    int   rateKout, rateKin;
    FILE *fbands = NULL;
    int   bands_resample = 0;
    int   verbose = 0;
    
    char* opt_string = "ho:";
    struct option long_options[] = {
        { "Fs", required_argument, &set_fs, 1 },
        { "rateK", no_argument, &rateK, 1 },
        { "perframe", no_argument, &perframe, 1 },
        { "newamp1vq", no_argument, &newamp1vq, 1 },
        { "rateKdec", required_argument, &rate_K_dec, 1 },
        { "rateKout", required_argument, &rateKout, 1 },
        { "rateKin", required_argument, &rateKin, 1 },
        { "bands",required_argument, &bands, 1 },
        { "bands_lower",required_argument, &bands_lower_en, 1 },
        { "bands_resample", no_argument, &bands_resample, 1 },
        { "lpc", required_argument, &lpc_model, 1 },
        { "nw", required_argument, &cnw, 1 },                       // flag if we wnat to customize lsp measurement window
        { "minmax", required_argument, &minmax, 1},                 // reads input file and counts minimum and maximum value for each lsp in order
        { "sqd", required_argument, &sqd, 1},			    // counts scalar quantizer distribution for given logic and dumps it into lspsqd.csv
        { "sqdd", required_argument, &sqddopt, 1},                  // sets scalar quantization dimension, defaults to 255 in --minmax, it is written in minmax.csv to read in other commands
        { "lspvt", required_argument, &lspvt, 1},                   // reads minmax.csv and reads input file and saves all subsequent vectors into lspvt.csv
        { "opti", required_argument, &opti, 1},                     // reads minmax.csv and --cdb count vectors from lspvt.csv, than counts distribution of closest matches of vectors built on input file, outputs lspvtopt.csv
        { "dbl", no_argument, &dbl, 1},                             // used with opti switches double closest vector to count
        { "clr", no_argument, &clr, 1},                             // used with opti clears vectors count in the vector table read, usefull for optimisation on joined previously optimised tables
        { "vq10", required_argument, &vq10, 1},                     // vq10 quantizer, number tells the "sysytem", together with --cdb N reads N vectors from lspvtopt.csv codebook, together with --lbg counts next LBG step to codebook.csv
        { "cdb", required_argument, &cdb, 1},                       // --cdb X together with opti reads X vectors from lspvt.csv file, together with vq10 reads X vectors from lspvtopt.csv.
        { "cdb1", required_argument, &cdb1, 1},                     // tha same for first split vector on system 2
        { "cdb2", required_argument, &cdb2, 1},                     // the same for the second split vector on system 2
        { "weight", no_argument, &wgt, 1},                          // reads weights and can be used for weighted distance calcuation
        { "dfr", no_argument, &dfr, 1},                             // doubles frequency sampling on the output file
        { "lbg", no_argument, &lbg, 1},                             // works together with vq10 and performs next iteration of LBG to the codebook.csv
        { "verbose", no_argument, &verbose, 1},                     // verobose output, mostly counts frames
        { "trn", no_argument, &trn, 1},                             // switch on the third transcient state - mixed excitation
        { "thrs", required_argument, &thrs, 1},                     // V/UV treshold for cwsim (defaults to 6.0)
        { "lsp", no_argument, &lsp, 1 },
        { "lspd", no_argument, &lspd, 1 },
        { "lspvq", no_argument, &lspvq, 1 },
        { "lspjmv", no_argument, &lspjmv, 1 },
        { "phase0", no_argument, &phase0, 1 },
        { "postfilter", no_argument, &postfilt, 1 },
        { "dmp", required_argument, &dmp, 1 },                      // my dump of the model vector for each frame into CSV format
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
        { "gain", required_argument, NULL, 0 },
        { "bpf", no_argument, &bpf_en, 1 },
        { "bpfb", no_argument, &bpfb_en, 1 },
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
            } else if(strcmp(long_options[option_index].name, "nw") == 0) {
                custom_nw = atoi(optarg);
            } else if(strcmp(long_options[option_index].name, "sqd") == 0) {
                statt = atoi(optarg);
            } else if(strcmp(long_options[option_index].name, "sqdd") == 0) {
                sqdd = atoi(optarg);
            } else if(strcmp(long_options[option_index].name, "cdb") == 0) {
                cdbsz = atoi(optarg);
            } else if(strcmp(long_options[option_index].name, "cdb1") == 0) {
                cdbsz1 = atoi(optarg);
            } else if(strcmp(long_options[option_index].name, "cdb2") == 0) {
                cdbsz2 = atoi(optarg);
            } else if(strcmp(long_options[option_index].name, "thrs") == 0) {
                v_thresh = atof(optarg);
            } else if(strcmp(long_options[option_index].name, "minmax") == 0) {
                statt = atoi(optarg);
            #ifdef DUMP
            } else if(strcmp(long_options[option_index].name, "dump") == 0) {
                if (dump)
	            dump_on(optarg);
            #endif
            } else if(strcmp(long_options[option_index].name, "lsp") == 0
                  || strcmp(long_options[option_index].name, "lspd") == 0
                  || strcmp(long_options[option_index].name, "lspvq") == 0) {
	        assert(order == LPC_ORD);
            } else if(strcmp(long_options[option_index].name, "rateKdec") == 0) {
                rate_K_dec = atoi(optarg);
                fprintf(stderr, "rate_K_dec: %d\n", rate_K_dec);
	    } else if(strcmp(long_options[option_index].name, "rateKout") == 0) {
                /* read model records from file or stdin */
                if ((frateK = fopen(optarg,"wb")) == NULL) {
	            fprintf(stderr, "Error opening output rateK file: %s: %s\n",
		        optarg, strerror(errno));
                    exit(1);
                }
                fprintf(stderr, "each record is %d bytes\n", (int)(K*sizeof(float)));
	    } else if(strcmp(long_options[option_index].name, "rateKin") == 0) {
                /* read model records from file or stdin */
                if ((frateKin = fopen(optarg,"rb")) == NULL) {
	            fprintf(stderr, "Error opening input rateK file: %s: %s\n",
		        optarg, strerror(errno));
                    exit(1);
                }
                fprintf(stderr, "each record is %d bytes\n", (int)(K*sizeof(float)));
	    } else if(strcmp(long_options[option_index].name, "bands") == 0) {
                /* write mel spaced band energies to file or stdout */
                if ((fbands = fopen(optarg,"wb")) == NULL) {
	            fprintf(stderr, "Error opening bands file: %s: %s\n",
		        optarg, strerror(errno));
                    exit(1);
                }
	    } else if(strcmp(long_options[option_index].name, "bands_lower") == 0) {
		bands_lower = atof(optarg);
		fprintf(stderr, "bands_lower: %f\n", bands_lower);
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
            } else if(strcmp(long_options[option_index].name, "dmp") == 0) {
                if ((fdmp = fopen(optarg,"at")) == NULL) {
                    fprintf(stderr, "Error opening my dump file: %s: %s. \n", optarg, strerror(errno));
                    exit(1);
                }
            } else if(strcmp(long_options[option_index].name, "lspvt") == 0) {
                statt = atoi(optarg);
                if (statt == 1 || statt == 4) {
                    if ((fdmp = fopen("lspvt.csv","wt")) == NULL) {
                        fprintf(stderr, "Error opening lspvt.csv %s. \n", strerror(errno));
                        exit(1);
                    }
                } else if (statt == 2) {
                    if ((fdmp1 = fopen("frmpos.csv","wt")) == NULL) {
                        fprintf(stderr, "Error opening frmpos.csv %s. \n", strerror(errno));
                        exit(1);
                    }
                    if ((fdmp2 = fopen("frmwdt.csv","wt")) == NULL) {
                        fprintf(stderr, "Error opening frmwdt %s. \n", strerror(errno));
                        exit(1);
                    }
                } else {
                    exit(1);
                }
            } else if(strcmp(long_options[option_index].name, "opti") == 0) {
                statt = atoi(optarg);
                if (statt == 1 || statt == 4) {
                    if ((fdmp = fopen("lspvt.csv","rt")) == NULL) {
                        fprintf(stderr, "Error opening lspvt.csv %s. \n", strerror(errno));
                        exit(1);
                    }
                    vtl1 = vtl2 = 1L;
                    vtl = 32768L;
                } else if (statt == 2) {
                    if ((fdmp1 = fopen("frmpos.csv","rt")) == NULL) {
                        fprintf(stderr, "Error opening frmpos.csv %s. \n", strerror(errno));
                        exit(1);
                    }
                    vtl1 = 4098L;
                    if ((fdmp2 = fopen("frmwdt.csv","rt")) == NULL) {
                        fprintf(stderr, "Error opening frmwdt.csv %s. \n", strerror(errno));
                        exit(1);
                    }
                    vtl2 = 4096L;
                }
            } else if(strcmp(long_options[option_index].name, "vq10") == 0) {
                statt = atoi(optarg);
                if (statt == 1 || statt == 4) {
                    if ((fdmp = fopen("lspvtopt.csv","rt")) == NULL) {
                        fprintf(stderr, "Error opening lspvtopt.csv %s. \n", strerror(errno));
                        exit(1);
                    }
                    vtl = 4096L;
                } else if (statt == 2) {
                    if ((fdmp1 = fopen("frmposopt.csv","rt")) == NULL) {
                        fprintf(stderr, "Error opening frmpos.csv %s. \n", strerror(errno));
                        exit(1);
                    }
                    if ((fdmp2 = fopen("frmwdtopt.csv","rt")) == NULL) {
                        fprintf(stderr, "Error opening frmwdt.csv %s. \n", strerror(errno));
                        exit(1);
                    }
                    vtl1 = vtl2 = 4098;
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
	    } else if(strcmp(long_options[option_index].name, "gain") == 0) {
		gain = atof(optarg);
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
                fprintf(stderr, "each model record is %d bytes\n", (int)sizeof(MODEL));
	    } else if(strcmp(long_options[option_index].name, "modelin") == 0) {
                /* read model records from file or stdin */
                modelin = 1;
                if (strcmp(optarg, "-") == 0) fmodelin = stdin;
	        else if ((fmodelin = fopen(optarg,"rb")) == NULL) {
	            fprintf(stderr, "Error opening modelin file: %s: %s\n",
		        optarg, strerror(errno));
                    exit(1);
                }
                fprintf(stderr, "each model record is %d bytes\n", (int)sizeof(MODEL));
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

    int i,m;

    float lspmin[order];
    float lspmax[order];
    float lspspn[order];

    if (sqd || lspvt || opti || vq10) {
        if ((faux = fopen("minmax.csv","rt")) == NULL) {
	    fprintf(stderr, "Error opening minmax.csv file: %s.\n",strerror(errno));
            exit(1);
        }
        n = fscanf(faux, "%d\n", &i);
        if (!sqddopt) sqdd = i;
        for(i=0; i<order; i++) {
            n = fscanf(faux, "%f;%f\n", &lspmin[i], &lspmax[i]);
            lspspn[i] = lspmax[i] - lspmin[i];
        }
        fclose(faux);
        faux = NULL;
    } else {
        for(i=0; i<order; i++) {
            lspmin[i] = 100.0;
            lspmax[i] = 0.0;
        }
    }


    if (statt == 2) {
        if (cdb) {
            vtl1 = cdbsz;
            vtl2 = cdbsz;
        } else {
            if (cdb1) vtl1 = cdbsz1;
            if (cdb2) vtl2 = cdbsz2;
        }
    } else {
        if (cdb) vtl = cdbsz;
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

    if (cnw) c2const.nw = custom_nw;

    short buf[N_SAMP];	/* input/output buffer                   */
    short buf16[N_SAMP * 2];	/* output buffer 16ksps                  */
    float buf_float[N_SAMP];
    float Sn[m_pitch];	/* float input speech samples            */
    float Sn_pre[m_pitch];	/* pre-emphasised input speech samples   */
    COMP  Sw[FFT_ENC];	/* DFT of Sn[]                           */
    codec2_fft_cfg  fft_fwd_cfg;
    codec2_fftr_cfg  fftr_fwd_cfg;
    codec2_fftr_cfg  fftr_inv_cfg;
    float w[m_pitch];	        /* time domain hamming window            */
    float W[FFT_ENC];	/* DFT of w[]                            */
    MODEL model;
    float Pn[2*N_SAMP];	/* trapezoidal synthesis window          */
    float Sn_[2*N_SAMP];	/* synthesised speech */
    long  l;
    float f, g;
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
    float hpf_states[2];
    #if 0
    struct PEXP *pexp = NULL;
    struct AEXP *aexp = NULL;
    #endif
    float bpf_buf[BPF_N+N_SAMP];

    COMP Aw[FFT_ENC];
    COMP H[MAX_AMP];

    float sd_sum = 0.0; int sd_frames = 0;

    if (lbg) vtlb = vtl; else vtlb = 1L;

    int vtraw[vtl][order+1];
    int vtraw1[vtl1][order/2+1];
    int vtraw2[vtl2][order/2+1];
    int lbgt[vtlb][order+1];
    int maxlbg = 0;

    int lspsqd[order][sqdd];

    float wght[order];

    if (wgt) {
        if ((fwght = fopen("weights.csv","rt")) == NULL) {
           fprintf(stderr, "Error opening input weight.csv file: %s\n", strerror(errno));
           exit(1);
        }
        for(i=0;i<order;i++) n = fscanf(fwght, "%f;", &wght[i]);
        fclose(fwght);
    } else {
        for(i=0;i<order;i++) wght[i] = 1.0;
    }
       
    for(i=0; i<order; i++) for(m=0; m<sqdd; m++) lspsqd[i][m] = 0;
    for(l=0L; l<vtl; l++) for(i=0; i<=order;i++) lbgt[l][i] = 0;

    if (vq10) {
        if (statt == 1 || statt == 4) {
            for(l=0L; l<vtl; l++) {
                n = fscanf(fdmp, "%d;", &vtraw[l][order]);
                for(i=0; i<order; i++) n = fscanf(fdmp, "%d;", &vtraw[l][i]);
                n = fscanf(fdmp, "\n");
            }
            fclose(fdmp);
            fdmp = NULL;
            if (verbose) fprintf(stderr, "\n");
        } else if (statt == 2) {
            for(l=0L; l<vtl1; l++) {
                n = fscanf(fdmp1, "%d;", &vtraw1[l][order/2]);
                for(i=0; i<order/2; i++) n = fscanf(fdmp1, "%d;", &vtraw1[l][i]);
                n = fscanf(fdmp1, "\n");
            }
            for(l=0L; l<vtl2; l++) {
                n = fscanf(fdmp2, "%d;", &vtraw2[l][order/2]);
                for(i=0; i<order/2; i++) n = fscanf(fdmp2, "%d;", &vtraw2[l][i]);
                n = fscanf(fdmp2, "\n");
            }
            fclose(fdmp1);
            fdmp1 = NULL;
            fclose(fdmp2);
            fdmp2 = NULL;
            if (verbose) fprintf(stderr, "\n");
        }
    }

    if (opti) {
        if (statt == 1 || statt == 4) {
            for(l=0L; l<vtl; l++) {
                n = fscanf(fdmp, "%d;", &vtraw[l][order]);
                for(i=0; i<order; i++) n = fscanf(fdmp, "%d;", &vtraw[l][i]);
                if (clr) vtraw[l][order] = 0;
                n = fscanf(fdmp, "\n");
            }
            fclose(fdmp);
            fdmp = NULL;
            if (verbose) fprintf(stderr, "\n");
        } else if (statt == 2) {
            for(l=0L; l<vtl1; l++) {
                n = fscanf(fdmp1, "%d;", &vtraw1[l][order/2]);
                for(i=0; i<order/2; i++) n = fscanf(fdmp1, "%d;", &vtraw1[l][i]);
                n = fscanf(fdmp1, "\n");
                if (clr) vtraw1[l][order/2] = 0;
            }
            for(l=0L; l<vtl2; l++) {
                n = fscanf(fdmp2, "%d;", &vtraw2[l][order/2]);
                for(i=0; i<order/2; i++) n = fscanf(fdmp2, "%d;", &vtraw2[l][i]);
                n = fscanf(fdmp2, "\n");
                if (clr) vtraw2[l][order/2] = 0;
            }
            fclose(fdmp1);
            fdmp1 = NULL;
            fclose(fdmp2);
            fdmp2 = NULL;
            if (verbose) fprintf(stderr, "\n");
        }
    }

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
    hpf_states[0] = hpf_states[1] = 0.0;

    nlp_states = nlp_create(&c2const);

    ex_phase[0] = 0;
    Woe_[0] = Woe_[1] = 1.0;

    /* Initialise ------------------------------------------------------------*/

    fft_fwd_cfg = codec2_fft_alloc(FFT_ENC, 0, NULL, NULL);   /* fwd FFT,used in several places   */
    fftr_fwd_cfg = codec2_fftr_alloc(FFT_ENC, 0, NULL, NULL); /* fwd FFT,used in several places   */
    fftr_inv_cfg = codec2_fftr_alloc(FFT_DEC, 1, NULL, NULL); /* inverse FFT, used just for synth */
    codec2_fft_cfg phase_fft_fwd_cfg = codec2_fft_alloc(NEWAMP1_PHASE_NFFT, 0, NULL, NULL);
    codec2_fft_cfg phase_fft_inv_cfg = codec2_fft_alloc(NEWAMP1_PHASE_NFFT, 1, NULL, NULL);

    make_analysis_window(&c2const, fft_fwd_cfg, w, W);
    make_synthesis_window(&c2const, Pn);

    if (bpfb_en)
        bpf_en = 1;
    if (bpf_en) {
        for(i=0; i<BPF_N; i++)
            bpf_buf[i] = 0.0;
    }

    for(i=0; i<order; i++) {
        prev_lsps_dec[i] = i*PI/(order+1);
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
	mel_sample_freqs_kHz(rate_K_sample_freqs_kHz, NEWAMP1_K, ftomel(200.0), ftomel(3700.0) );
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

        if (bpf_en) {
            /* filter input speech to create buf_float_bpf[], this is fed to the
               LPC modelling.  Unfiltered speech in in buf_float[], which is
               delayed to match that of the BPF */

            /* BPF speech */

            for(i=0; i<BPF_N; i++)
                bpf_buf[i] =  bpf_buf[N_SAMP+i];
            for(i=0; i<N_SAMP; i++)
                bpf_buf[BPF_N+i] = buf_float[i];
            if (bpfb_en)
                inverse_filter(&bpf_buf[BPF_N], bpfb, N_SAMP, buf_float, BPF_N);
            else
                inverse_filter(&bpf_buf[BPF_N], bpf, N_SAMP, buf_float, BPF_N);
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

	/* estimate voicing - do this all the time so model.voicing
	 * is set, useful for machine learning work */
	snr = est_voicing_mbe(&c2const, &model, Sw, W, v_thresh);

        if (trn && (model.voiced != prev_vo)) model.transcient = 1;
        else model.transcient = 0;
        prev_vo = model.voiced;

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
            float ak_[order+1];

            e = speech_to_uq_lsps(lsps, ak, Sn, w, m_pitch, order);
            for(i=0; i<order; i++)
                lsps_[i] = lsps[i];


/*******************************************************************************************************
            THIS IS MY PART WHERE WHOLE PARAMETER SET FOR A FRAME IS CALCULATED
            SO HERE MOST OF THE ADDITIONS WORK
********************************************************************************************************/

            if (dmp) {
                framecounter++;
                fprintf(fdmp, "%d;%d;%f;%f;%f", framecounter, model.voiced, snr, model.Wo, e);
                for(i=0; i<order; i++) fprintf(fdmp, ";%f", lsps[i]);
                fprintf(fdmp, "\n");
            }

            if (minmax) {
                if (statt == 1) {
                    for (i=0; i<order; i++) {
                        if (lsps[i]<lspmin[i]) lspmin[i]=lsps[i];
                        if (lsps[i]>lspmax[i]) lspmax[i]=lsps[i];
                    }
                } else if (statt == 2) {
                    for(i=0; i<order; i+=2) {
                        f = (lsps[i+1]-lsps[i])/2.0+lsps[i];
                        if (f < lspmin[i]) lspmin[i]=f;
                        if (f > lspmax[i]) lspmax[i]=f;
                        f = (lsps[i+1]-lsps[i])/2.0;
                        if (f < lspmin[i+1]) lspmin[i+1]=f;
                        if (f > lspmax[i+1]) lspmax[i+1]=f;
                    }
                } else if (statt == 3) {
                    for (i=0; i<order; i+=2) {
                        if (i==0) f = (lsps[1]-lsps[0])/2.0+lsps[0];
                        else f = ((lsps[i+1]-lsps[i])/2.0+lsps[i])-((lsps[i-1]-lsps[i-2])/2.0+lsps[i-2]);
                        if (f < lspmin[i/2]) lspmin[i/2]=f;
                        if (f > lspmax[i/2]) lspmax[i/2]=f;
                    }
                    for (i=0; i<order; i+=2) {
                        if (i==0) f = (lsps[1]-lsps[0])/2.0;
                        else f = ((lsps[i+1]-lsps[i])/2.0-((lsps[i-1]-lsps[i-2])/2.0));
                        if (f < lspmin[i/2+order/2]) lspmin[i/2+order/2] = f;
                        if (f > lspmax[i/2+order/2]) lspmax[i/2+order/2] = f;
                    }
                }
                framecounter++;
                if (verbose) fprintf(stderr, "\r%d", framecounter);
            }

            if (sqd) {
                if (statt == 1) {
                    for(i=0; i<order; i++) lspsqd[i][(int)floor((lsps[i]-lspmin[i])/(lspmax[i]-lspmin[i])*sqdd)]++;
                } else if (statt == 2) {
                    for(i=0; i<order; i+=2) {
                        f = (lsps[i+1]-lsps[i])/2.0+lsps[i];
                        m = (int)floor((f-lspmin[i])/(lspmax[i]-lspmin[i])*sqdd);
                        lspsqd[i][m]++;
                        f = (lsps[i+1]-lsps[i])/2.0;
                        m = (int)floor((f-lspmin[i+1])/(lspmax[i+1]-lspmin[i+1])*sqdd);
                        lspsqd[i+1][m]++;
                    }
                } else if (statt == 3) {
                    for(i=0; i<order; i+=2) {
                        if (i==0) f = (lsps[1]-lsps[0])/2.0+lsps[0];
                        else f = ((lsps[i+1]-lsps[i])/2.0+lsps[i])-((lsps[i-1]-lsps[i-2])/2.0+lsps[i-2]);
                        m = (int)floor((f-lspmin[i/2])/(lspmax[i/2]-lspmin[i/2])*sqdd);
                        lspsqd[i/2][m]++;
                        if (i==0) f = (lsps[1]-lsps[0]);
                        else f = ((lsps[i+1]-lsps[i])/2.0-((lsps[i-1]-lsps[i-2])/2.0));
                        m = (int)floor((f-lspmin[i/2+order/2])/(lspmax[i/2+order/2]-lspmin[i/2+order/2])*sqdd);
                        lspsqd[i/2+order/2][m]++;
                    }
                }
                framecounter++;
                if (verbose) fprintf(stderr, "\r%d", framecounter);
            }

            if (lspvt) {
                framecounter++;
                if (statt == 1) {
                    fprintf(fdmp, "0;");
                    for(i=0; i<order; i++) fprintf(fdmp, "%d;", (int)floor((lsps[i]-lspmin[i])/(lspmax[i]-lspmin[i])*sqdd));
                    fprintf(fdmp, "\n");
                } else if (statt == 2) {
                    fprintf(fdmp1, "0;");
                    fprintf(fdmp2, "0;");
                    for(i=0; i<order; i+=2) {
                        f = (lsps[i+1]-lsps[i])/2.0+lsps[i];
                        m = (int)floor((f-lspmin[i])/(lspmax[i]-lspmin[i])*sqdd);
                        fprintf(fdmp1, "%d;", m);
                        f = (lsps[i+1]-lsps[i])/2.0;
                        m = (int)floor((f-lspmin[i+1])/(lspmax[i+1]-lspmin[i+1])*sqdd);
                        fprintf(fdmp2, "%d;", m);
                    }
                    fprintf(fdmp1, "\n");
                    fprintf(fdmp2, "\n");
                } else if (statt == 3) {
                    for(i=0; i<order/2; i++) {
                        if (i==0) f = (lsps[1]-lsps[0])/2.0+lsps[0];
                        else f = ((lsps[i*2+1]-lsps[i*2])/2.0+lsps[i*2])-((lsps[i*2-1]-lsps[i*2-2])/2.0+lsps[i*2-2]);
                        m = (int)floor((f-lspmin[i])/(lspmax[i]-lspmin[i])*sqdd);
                        fprintf(fdmp, "%d;", m);
                    }
                    for(i=0; i<order/2; i++) {
                        if (i==0) f = (lsps[1]-lsps[0])/2.0;
                        else f = ((lsps[i*2+1]-lsps[i*2])/2.0-((lsps[i*2-1]-lsps[i*2-2])/2.0));
                        m = (int)floor((f-lspmin[i+order/2])/(lspmax[i+order/2]-lspmin[i+order/2])*sqdd);
                        fprintf(fdmp, "%d;", m);
                    }
                    fprintf(fdmp, "\n");
                }
                if (verbose) fprintf(stderr, "\r%d", framecounter);
            }

            if (opti) {
                if (statt == 1) {
                    framecounter++;
                    sqermn = 2147483647L;
                    bstix = 0L;
                    bstix2 = 0L;
                    sqerr = 0L;
                    for(l=0L; l<vtl; l++) {
                        sqerr = 0L;
                        for(i=0; i<order; i++) {                        
                            if (lsps[i] < lspmin[i])
                                m = 0;
                            else if (lsps[i] > lspmax[i])
                                m = sqdd;
                            else 
                                m = (int)floor((lsps[i]-lspmin[i])/lspspn[i]*sqdd);
                            m = m - vtraw[l][i];
                            sqerr += (m * m);
                            if (sqerr >= sqermn) break;
                        }
                        if(sqerr < sqermn) {
                            sqermn = sqerr;
                            if (dbl) bstix2 = bstix;
                            bstix = l;
                        }
                    }
                    vtraw[bstix][order]++;
                    if (dbl) vtraw[bstix2][order]++;
                    if (vtraw[bstix][order] > maxoc) maxoc = vtraw[bstix][order];
                    if (verbose) fprintf(stderr, "\r%d", framecounter);
                } else if (statt == 2) {
                    framecounter++;
                    sqermn1 = 2147483647L;
                    bstix1 = 0L;
                    sqermn2 = 2147483647L;
                    bstix2 = 0L;
                    for(l=0L; l<vtl1; l++) {
                        sqerr1 = 0L;
                        for(i=0; i<order; i+=2) {
                            f = (lsps[i+1]-lsps[i])/2.0+lsps[i];
                            if (f < lspmin[i])
                                m = 0;
                            else if (f > lspmax[i])
                                m = sqdd;
                            else 
                                m = (int)floor((f-lspmin[i])/(lspmax[i]-lspmin[i])*sqdd);
                            m = m - vtraw1[l][i/2];
                            sqerr1 += m * m;
                            if (sqerr1 >= sqermn1) break;
                        }
                        if(sqerr1 < sqermn1) {
                            sqermn1 = sqerr1;
                            bstix1 = l;
                        }
                    }
                    for(l=0L; l<vtl2; l++) {
                        sqerr2 = 0L;
                        for(i=0; i<order; i+=2) {
                            f = (lsps[i+1]-lsps[i])/2.0;
                            if (f < lspmin[i+1])
                                m = 0;
                            else if (f > lspmax[i+1])
                                m = sqdd;
                            else
                                m = (int)floor((f-lspmin[i+1])/(lspmax[i+1]-lspmin[i+1])*sqdd);
                            m = m - vtraw2[l][i/2];
                            sqerr2 += m * m;
                            if (sqerr2 >= sqermn2) break;
                        }
                        if(sqerr2 < sqermn2) {
                            sqermn2 = sqerr2;
                            bstix2 = l;
                        }
                    }
                    vtraw1[bstix1][order/2]++;
                    vtraw2[bstix2][order/2]++;
                    if (vtraw1[bstix1][order/2] > maxoc1) maxoc1 = vtraw1[bstix1][order/2];
                    if (vtraw2[bstix2][order/2] > maxoc2) maxoc2 = vtraw2[bstix2][order/2];
                    if (verbose) fprintf(stderr, "\r%d", framecounter);
                } else if (statt == 4) {
                    framecounter++;
                    sqermn = 2147483647L;
                    bstix = 0L;
                    sqerr = 0L;
                    for(l=0L; l<vtl; l++) {
                        sqerr = 0L;
                        for(i=0; i<order/2; i++) {                        
                            m = (int)floor((lsps[i*2]-lspmin[i*2])/(lspmax[i*2]-lspmin[i*2])*sqdd);
                            n = (int)floor((lsps[i*2+1]-lspmin[i*2+1])/(lspmax[i*2+1]-lspmin[i*2+1])*sqdd);
                            m = ((n-m)/2+m) - ((vtraw[l][i*2+1]-vtraw[l][i*2])/2+vtraw[l][i*2]);
                            sqerr += (m * m);
                            m = (int)floor((lsps[i*2]-lspmin[i*2])/(lspmax[i*2]-lspmin[i*2])*sqdd);
                            n = (int)floor((lsps[i*2+1]-lspmin[i*2+1])/(lspmax[i*2+1]-lspmin[i*2+1])*sqdd);
                            m = ((n-m)/2) - ((vtraw[l][i*2+1]-vtraw[l][i*2])/2);
                            sqerr += (m * m);
                            if (sqerr >= sqermn) break;
                        }
                        if(sqerr < sqermn) {
                            sqermn = sqerr;
                            bstix = l;
                        }
                    }
                    vtraw[bstix][order]++;
                    if (vtraw[bstix][order] > maxoc) maxoc = vtraw[bstix][order];
                    if (verbose) fprintf(stderr, "\r%d", framecounter);
                }
            }

            if (vq10) {
                framecounter++;
                if (statt == 1) {
                    sqermn = 2147483647L;
                    bstix = 0L;
                    sqerr = 0L;
                    for(l=0L; l<vtl; l++) {
                        sqerr = 0L;
                        for(i=0; i<order; i++) {                        
                            if (lsps[i] < lspmin[i])
                                m = 0;
                            else if (lsps[i] > lspmax[i])
                                m = sqdd;
                            else
                                m = (int)floor((lsps[i]-lspmin[i])/(lspmax[i]-lspmin[i])*sqdd);
                            m = m - vtraw[l][i];
                            sqerr += (m * m);
                            if (sqerr >= sqermn) break;
                        }
                        if(sqerr < sqermn) {
                            sqermn = sqerr;
                            bstix = l;
                        }
                    }
                    for (i=0; i<order; i++) lsps_[i] = ((float)vtraw[bstix][i]/sqdd)*(lspmax[i]-lspmin[i])+lspmin[i];
	            lsp_to_lpc(lsps_, ak_, order);
                    if (lbg) {
                        lbgt[bstix][order]++;
                        if (lbgt[bstix][order] > maxlbg) maxlbg = lbgt[bstix][order];
                        for(i=0;i<order;i++) {
                            m = (int)floor((lsps[i]-lspmin[i])/(lspmax[i]-lspmin[i])*sqdd);
                            lbgt[bstix][i] += m;
                            m = m - vtraw[bstix][i];
                            lbgd0 += m * m;
                        }
                    }
                } else if (statt == 2) {
                    sqermn1 = 2147483647L;
                    bstix1 = 0L;
                    for(l=0L; l<vtl1; l++) {
                        sqerr1 = 0L;
                        for(i=0; i<order; i+=2) {
                            f = (lsps[i+1]-lsps[i])/2.0+lsps[i];
                            if (f < lspmin[i])
                                m = 0;
                            else if (f > lspmax[i])
                                m = sqdd;
                            else
                                m = (int)floor((f-lspmin[i])/(lspmax[i]-lspmin[i])*sqdd);
                            m = m - vtraw1[l][i/2];
                            sqerr1 += m * m;
                            if (sqerr1 >= sqermn1) break;
                        }
                        if(sqerr1 < sqermn1) {
                            sqermn1 = sqerr1;
                            bstix1 = l;
                        }
                    }
                    sqermn2 = 2147483647L;
                    bstix2 = 0L;
                    for(l=0L; l<vtl2; l++) {
                        sqerr2 = 0L;
                        for(i=0; i<order; i+=2) {
                            f = (lsps[i+1]-lsps[i])/2.0;
                            if (f < lspmin[i+1])
                                m = 0;
                            else if (f > lspmax[i+1])
                                m = sqdd;
                            else
                                m = (int)floor((f-lspmin[i+1])/(lspmax[i+1]-lspmin[i+1])*sqdd);
                            m = m - vtraw2[l][i/2];
                            sqerr2 += m * m;
                            if (sqerr2 >= sqermn2) break;
                        }
                        if(sqerr2 < sqermn2) {
                            sqermn2 = sqerr2;
                            bstix2 = l;
                        }
                    }
                    for(i=0; i<order; i+=2) {
                        f = ((float)vtraw1[bstix1][i/2]/sqdd)*(lspmax[i]-lspmin[i])+lspmin[i];
                        g = ((float)vtraw2[bstix2][i/2]/sqdd)*(lspmax[i+1]-lspmin[i+1])+lspmin[i+1];
                        lsps_[i] = f - g;
                        lsps_[i+1] = f + g;
                    }
		            lsp_to_lpc(lsps_, ak_, order);
                } else if (statt == 4) {
                    sqermn = 2147483647L;
                    bstix = 0L;
                    sqerr = 0L;
                    for(l=0L; l<vtl; l++) {
                        sqerr = 0L;
                        for(i=0; i<order/2; i++) {                        
                            m = (int)floor((lsps[i*2]-lspmin[i*2])/(lspmax[i*2]-lspmin[i*2])*sqdd);
                            n = (int)floor((lsps[i*2+1]-lspmin[i*2+1])/(lspmax[i*2+1]-lspmin[i*2+1])*sqdd);
                            m = ((n-m)/2+m) - ((vtraw[l][i*2+1]-vtraw[l][i*2])/2+vtraw[l][i*2]);
                            sqerr += (m * m);
                            m = (int)floor((lsps[i*2]-lspmin[i*2])/(lspmax[i*2]-lspmin[i*2])*sqdd);
                            n = (int)floor((lsps[i*2+1]-lspmin[i*2+1])/(lspmax[i*2+1]-lspmin[i*2+1])*sqdd);
                            m = ((n-m)/2) - ((vtraw[l][i*2+1]-vtraw[l][i*2])/2);
                            sqerr += (m * m);
                            if (sqerr >= sqermn) break;
                        }
                        if(sqerr < sqermn) {
                            sqermn = sqerr;
                            bstix = l;
                        }
                    }
                    for (i=0; i<order; i++) lsps_[i] = ((float)vtraw[bstix][i]/sqdd)*(lspmax[i]-lspmin[i])+lspmin[i];
		            lsp_to_lpc(lsps_, ak_, order);
                }
                if (verbose) fprintf(stderr, "\r%d", framecounter);
            }


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

            if (lsp || lspd || lspjmv || vq10) {
                sd_sum += spectral_dist(ak, ak_, order, fft_fwd_cfg, FFT_ENC);
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

	/* LPCNet type mel spaced band ML data */
	float bands_mean = 0.0;
	if (fbands) {
	    float bandE[LPCNET_FREQ_MAX_BANDS];
            float freqkHz[LPCNET_FREQ_MAX_BANDS];
	    int nbands = lpcnet_compute_band_energy(bandE, freqkHz, Sw, Fs, FFT_ENC);
	    for(int i=0; i<nbands; i++)
		bands_mean += bandE[i];
	    bands_mean /= nbands;
	    //fprintf(stderr, "bands_mean: %f bands_lower %f\n", bands_mean,  bands_lower);
	    if (bands_mean > bands_lower)
 		assert(fwrite(bandE, sizeof(float), nbands, fbands) == nbands);
            // optionally reconstruct [Am} by linear interpolation of band energies,
            // this doesn't sound very Good
            if (bands_resample)
                resample_rate_L(&c2const, &model, &bandE[1], &freqkHz[1], nbands-2);
	}

	/*------------------------------------------------------------*\

	            Optional newamp1 simulation, as used in 700C

	\*------------------------------------------------------------*/

        if (rateK) {
            float rate_K_vec[K];
            resample_const_rate_f(&c2const, &model, rate_K_vec, rate_K_sample_freqs_kHz, K);

	    if (frateK != NULL)
		assert(fwrite(rate_K_vec, sizeof(float), K, frateK) == K);
	    
	    if (frateKin != NULL) {
		assert(fread(rate_K_vec, sizeof(float), K, frateKin) == K);
		/* apply newamp1 postfilter - this helped male samples with VQVAE work */
                float sum = 0.0;
                for(int k=0; k<K; k++)
                    sum += rate_K_vec[k];
                float mean = sum/K;
                float rate_K_vec_no_mean[K];
                for(int k=0; k<K; k++)
                    rate_K_vec_no_mean[k] = rate_K_vec[k] - mean;
		post_filter_newamp1(rate_K_vec_no_mean,  rate_K_sample_freqs_kHz, K, 1.5);
                for(int k=0; k<K; k++)
                    rate_K_vec[k] = rate_K_vec_no_mean[k] +  mean;
	    }
	    
            float rate_K_vec_[K];
            if (newamp1vq) {
                /* remove mean */
                float sum = 0.0;
                for(int k=0; k<K; k++)
                    sum += rate_K_vec[k];
                float mean = sum/K;
                float rate_K_vec_no_mean[K];
                for(int k=0; k<K; k++)
                    rate_K_vec_no_mean[k] = rate_K_vec[k] - mean;

		newamp1_eq(rate_K_vec_no_mean, eq, K, 1);

                /* two stage VQ */
                float rate_K_vec_no_mean_[K]; int indexes[2];
                rate_K_mbest_encode(indexes, rate_K_vec_no_mean, rate_K_vec_no_mean_, K, NEWAMP1_VQ_MBEST_DEPTH);
                for(int k=0; k<K; k++)
                    rate_K_vec_[k] = rate_K_vec_no_mean_[k] + mean;

                /* running sum of squared error for variance calculation */
                for(int k=0; k<K; k++)
                    se += pow(rate_K_vec_no_mean[k]-rate_K_vec_no_mean_[k],2.0);
                nse += K;
            }
            else {
                for(int k=0; k<K; k++)
                    rate_K_vec_[k] = rate_K_vec[k];
            }

	    if (frateKWov != NULL) {
		/* We use standard nb_features=55 feature records for compatability with train_lpcnet.py */
		float features[55] = {0};
		/* just using 18/20 for compatability with LPCNet, coarse scaling for NN imput */
		for(int i=0; i<18; i++)
		    features[i] = (rate_K_vec_[i]-30)/40;
		// keep in range of 40 ... 255 for pitch embedding
		int pitch_index = 21 + 2.0*M_PI/model.Wo;
		features[36] = 0.02*(pitch_index-100);
		//features[36] = (model.Wo - c2const.Wo_min)/(c2const.Wo_max - c2const.Wo_min) - 0.5;
		features[37] = model.voiced;
		if (first)
		    features[18] = -0.9;
		if (lpc_model) {
		    MODEL model_;
		    model_.Wo = model.Wo;
		    model_.L  = model.L;
		    model_.voiced = model.voiced;
		    float Rk[order+1], ak[order+1];
		    resample_rate_L(&c2const, &model_, rate_K_vec_, rate_K_sample_freqs_kHz, K);
		    determine_autoc(&c2const, Rk, order, &model_, NEWAMP1_PHASE_NFFT, phase_fft_fwd_cfg, phase_fft_inv_cfg);
		    /* -40 dB noise floor and Lag windowing from LPCNet/freq.c - helps reduce large spikes in spectrum when LPC
                       analysis loses it. */
		    Rk[0] += Rk[0]*1e-4 + 320/12/38.;
		    for (i=1;i<order+1;i++) Rk[i] *= (1 - 6e-5*i*i);
		    levinson_durbin(Rk, ak, order);

		    for(int i=0; i<order; i++)
			features[18+i] = ak[i+1];
		}
		fwrite(features, 55, sizeof(float), frateKWov);
	    }

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
	    
            resample_rate_L(&c2const, &model, rate_K_vec_, rate_K_sample_freqs_kHz, K);
        }

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
                    dump_quantised_model(&model_dec[i]);
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
                }

                if (postfilt)
                    postfilter(&model_dec[i], &bg_est);

                synth_one_frame(n_samp, fftr_inv_cfg, buf, &model_dec[i], Sn_, Pn, prede, &de_mem, gain);

                if (fout != NULL) {
                    if (!dfr) {
                        fwrite(buf,sizeof(short),N_SAMP,fout);
                    } else {
                        for(m=0; m<N_SAMP-1; m++) {
                            buf16[m*2] = buf[m];
                            buf16[m*2+1] = (buf[m+1]+buf[m])/2;
                        }
                        buf16[N_SAMP*2-2] = buf[N_SAMP-1];
                        buf16[N_SAMP*2-1] = buf[N_SAMP-1];
                        fwrite(buf16,sizeof(short),N_SAMP*2,fout);
                    }
                }

                if (modelout) {
		    /* optionally don't write to filter out low energy frames */
	            if (bands) {
	                if (bands_mean > bands_lower)
	                    fwrite(&model_dec[i],sizeof(MODEL),1,fmodelout);
		    }
		    else
		        fwrite(&model_dec[i],sizeof(MODEL),1,fmodelout);
		}
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

    if (fout != NULL) fclose(fout);

    if (minmax) {
        if (fdmp == NULL) fdmp = fopen("minmax.csv", "wt");
        fprintf(fdmp, "%d\n", sqdd);
        for(i=0; i<order; i++) fprintf(fdmp, "%f;%f\n", lspmin[i], lspmax[i]);
    }

    if (fdmp != NULL) fclose(fdmp);
    if (fdmp1 != NULL) fclose(fdmp1);
    if (fdmp2 != NULL) fclose(fdmp2);


    if (sqd) {
        if ((faux = fopen("lspsqd.csv","wt")) == NULL) {
	    fprintf(stderr, "Error opening lspsqd.csv file: %s.\n",strerror(errno));
            exit(1);
        }
        for(i=0; i<sqdd; i++) {
            for(m=0; m<order; m++) fprintf(faux, "%d;", lspsqd[m][i]);
            fprintf(faux, "\n");
        }
        fclose(faux);
    }

    if (opti) {
        if (statt == 1) {
            frmcntl = 0L;
            if (faux != NULL) fclose(faux);
            if ((faux = fopen("lspvtopt.csv","wt")) == NULL) {
	        fprintf(stderr, "Error opening lspvtopt.csv file: %s.\n",strerror(errno));
                exit(1);
            }

            if (verbose) fprintf(stderr, "\nmaxoc=%d vtl=%ld  order=%d\n", maxoc, vtl, order);

            for(m=maxoc; m>0; m--) {
                for(l=0L; l<vtl; l++) {
                    if (vtraw[l][order] == m) {
                        frmcntl++;
                        fprintf(faux, "%d;", m);
                        for(i=0; i<order; i++) fprintf(faux, "%d;", vtraw[l][i]);
                        fprintf(faux, "\n");
                    }
                }
                if (verbose) fprintf(stderr, "\r%d", m);
            }
            fclose(faux);
            fprintf(stderr, "\r%ld meaningfull frames\n", frmcntl);
        } else if (statt == 2) {
            frmcntl = 0L;
            if (faux != NULL) fclose(faux);
            if ((faux = fopen("frmposopt.csv","wt")) == NULL) {
	        fprintf(stderr, "Error opening frmposopt.csv file: %s.\n",strerror(errno));
                exit(1);
            }
            if (verbose) fprintf(stderr, "\nmaxoc1=%d vtl1=%ld  order=%d\n", maxoc1, vtl1, order);
            for(m=maxoc1; m>0; m--) {
                for(l=0L; l<vtl1; l++) {
                    if (vtraw1[l][order/2] == m) {
                        frmcntl++;
                        fprintf(faux, "%d;", m);
                        for(i=0; i<order/2; i++) fprintf(faux, "%d;", vtraw1[l][i]);
                        fprintf(faux, "\n");
                    }
                }
                if (verbose) fprintf(stderr, "\r%d", m);
            }
            fclose(faux);
            faux = NULL;
            fprintf(stderr, "\r%ld meaningfull frames\n", frmcntl);
            frmcntl = 0L;
            if (faux != NULL) fclose(faux);
            if ((faux = fopen("frmwdtopt.csv","wt")) == NULL) {
	        fprintf(stderr, "Error opening frmposopt.csv file: %s.\n",strerror(errno));
                exit(1);
            }
            if (verbose) fprintf(stderr, "\nmaxoc2=%d vtl2=%ld  order=%d\n", maxoc2, vtl2, order);
            for(m=maxoc2; m>0; m--) {
                for(l=0L; l<vtl2; l++) {
                    if (vtraw2[l][order/2] == m) {
                        frmcntl++;
                        fprintf(faux, "%d;", m);
                        for(i=0; i<order/2; i++) fprintf(faux, "%d;", vtraw2[l][i]);
                        fprintf(faux, "\n");
                    }
                }
                if (verbose) fprintf(stderr, "\r%d", m);
            }
            fclose(faux);
            faux = NULL;
            fprintf(stderr, "\r%ld meaningfull frames\n", frmcntl);

        }
    }

    if (vq10 && lbg && (statt == 1)) {
        frmcntl = 0L;
        if (faux != NULL) fclose(faux);
        if ((faux = fopen("codebook.csv","wt")) == NULL) {
	        fprintf(stderr, "Error opening lspvtopt.csv file: %s.\n",strerror(errno));
            exit(1);
        }

        fprintf(stderr, "\nmax hits=%d vtl=%ld  order=%d\noverall codebook distortion = %f\n", maxlbg, vtl, order, lbgd0/framecounter);

        for(m=maxlbg; m>0; m--) {
            for(l=0L; l<vtl; l++) {
                if (lbgt[l][order] == m) {
                    frmcntl++;
                    fprintf(faux, "%d;", m);
                    for(i=0; i<order; i++) {
                        lbgt[l][i] /= lbgt[l][order];
                        fprintf(faux, "%d;", lbgt[l][i]);
                    }
                    fprintf(faux, "\n");
                }
            }
//            if (verbose) fprintf(stderr, "\r%d", m);
        }
        fclose(faux);
        fprintf(stderr, "\n%ld meaningfull frames\n", frmcntl);
    }


    if (lpc_model) {
    	fprintf(stderr, "LPC->{Am} SNR av: %5.2f dB over %d frames\n", sum_snr/frames, frames);
        if (lsp || lspd || lspjmv || vq10)
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
    if (frateKin != NULL) fclose(frateKin);
    if (ften_ms_centre != NULL) fclose(ften_ms_centre);
    if (fmodelout != NULL) fclose(fmodelout);
    if (fbands != NULL) fclose(fbands);
    if (frateKWov != NULL) fclose(frateKWov);

    return 0;
}

void synth_one_frame(int n_samp, codec2_fftr_cfg fftr_inv_cfg, short buf[], MODEL *model, float Sn_[],
                     float Pn[], int prede, float *de_mem, float gain)
{
    int     i;

    synthesise(n_samp, fftr_inv_cfg, Sn_, model, Pn, 1);
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
