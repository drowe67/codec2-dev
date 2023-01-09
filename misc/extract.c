/*
  extract.c
  david Rowe Jan 2019

  Vector Quantiser tool that pre-processes training data before vq_train.
*/

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <getopt.h>

#define NB_FEATURES 55 /* number of cols per row */

int main(int argc, char *argv[]) {
    FILE *fin, *fout;
    int st = 0;
    int en = 17;
    int stride = NB_FEATURES;
    float gain = 1.0;
    int frame_delay = 1;
    float pred = 0.0;
    int removemean = 0;
    float lower = -1E32;
    int writeall = 0;
    float dynamicrange = 100.0;
    
    static struct option long_options[] = {
        {"startcol",      required_argument, 0, 's'},
        {"endcol",        required_argument, 0, 'e'},
        {"stride",        required_argument, 0, 't'},
        {"gain",          required_argument, 0, 'g'},
        {"pred",          required_argument, 0, 'p'},
        {"delay",         required_argument, 0, 'd'},
        {"removemean",    no_argument,       0, 'm'},
        {"lower",         required_argument, 0, 'l'},
        {"dynamicrange", required_argument, 0, 'y'},
        {"writeall",      no_argument,       0, 'w'},
        {0, 0, 0, 0}
    };

    int opt_index = 0;
    int c;
    
    while ((c = getopt_long (argc, argv, "s:e:t:g:p:d:ml:y:", long_options, &opt_index)) != -1) {
        switch (c) {
        case 's':
            st = atoi(optarg);
            break;
        case 'e':
            en = atoi(optarg);
            break;
        case 't':
            stride = atoi(optarg);
            break;
        case 'g':
            gain = atof(optarg);
            break;
        case 'p':
            pred = atof(optarg);
            break;
        case 'd':
            frame_delay = atoi(optarg);
            break;
        case 'm':
            removemean = 1;
            break;
        case 'l':
            lower = atof(optarg);
            break;
        case 'w':
            writeall = 1;
            break;
        case 'y':
            dynamicrange = atof(optarg);
            break;
        default:
        helpmsg:
            fprintf(stderr, "usage: %s  -s startCol -e endCol [options] input.f32 output.f32\n"
                            "\n"
                            "-t strideCol           Vector dimenstion, e.g. 30\n"
                            "-g gain                Gain applied to vectors\n"
                            "-p predCoeff           Coefficient used for prediction, e.g. 0.9\n"
                            "-d framesDelay Delay   Delay (frames) between vectors used for prediction\n"
                            "--removemean           Remove mean from vectors\n"
                            "--lower minEnergy      Remove all vectors less than minEnergy\n"
                            "--writeall             Write all K outputs for each vectoir, even if EndCol-StartCol+1 < stride\n"
                            "--dynamicrange RangedB Clip min value of each vector to max - RangedB\n"
                            "input.f32 output.f32\n", argv[0]);
            exit(1);
        }
    }
    if ( (argc - optind) < 2) {
        fprintf(stderr, "Too few arguments\n");
        goto helpmsg;
    }
 
    if ((lower != -1E32) && (pred != 0.0)) {
        /* These two options cause unexpected jumps in predictive VQ training data */ 
        fprintf(stderr, "Warning - low energy vector removal and prediction enabled\n");
    }
    
    fin = fopen(argv[optind],"rb"); assert(fin != NULL);
    fout = fopen(argv[optind+1],"wb"); assert(fout != NULL);
    printf("extracting from %d to %d inclusive (stride %d) ... \n"
           "gain = %f pred = %f frame_delay = %d dynamic range = %f\n",
           st, en, stride, gain, pred, frame_delay, dynamicrange);
   
    float features[stride], features_prev[frame_delay][stride], delta[stride];
    int i,f,rd=0,wr=0;
    
    for(i=0; i<stride; i++)
        delta[i] = 0.0;
    for (f=0; f<frame_delay; f++)
        for(i=0; i<stride; i++)
            features_prev[f][i] = 0.0;

    while((fread(features, sizeof(float), stride, fin) == stride)) {
	float mean = 0.0;
	for(i=st; i<=en; i++)
	    mean += features[i];
	mean /= (en-st+1);
 	if (removemean) {
	    for(i=0; i<stride; i++)
		features[i] -= mean;
	}
        
	for(i=st; i<=en; i++) {
	    delta[i] = gain*(features[i] - pred*features_prev[frame_delay-1][i]);
	}
        
        if (dynamicrange != 100.0) {
            float max = -1E32;
	    for(i=st; i<=en; i++)
	        if (features[i] > max) max = features[i];
	    for(i=st; i<=en; i++) {
                float min = max-dynamicrange;
                if (features[i] < min ) features[i] = min;
            }        
        }
        
	if (mean > lower) {
	    if (writeall)
                fwrite(delta, sizeof(float), stride, fout);
            else
                fwrite(&delta[st], sizeof(float), en-st+1, fout);
	    wr++;
	}
	for (f=frame_delay-1; f>0; f--)
	    for(i=0; i<stride; i++)
		features_prev[f][i] = features_prev[f-1][i];
	for(i=0; i<stride; i++)
	    features_prev[0][i] = features[i];
        rd++;    
    }

    fclose(fin); fclose(fout);
    fprintf(stderr, "%d input %d extracted\n", rd, wr);
    return 0;
}

