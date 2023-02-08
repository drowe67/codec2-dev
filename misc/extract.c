/*
  extract.c
  David Rowe Jan 2019

  Vector Quantiser tool that pre-processes training data before vq_train.
*/

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <getopt.h>

#define NB_FEATURES 55 /* number of cols per row */

void restrict_dynamic_range(float features[], int st, int en, float dynamicrange) {
    int i;
    if (dynamicrange != 100.0) {
        float max = -1E32;
        for(i=st; i<=en; i++)
            if (features[i] > max) max = features[i];
        for(i=st; i<=en; i++) {
            float min = max-dynamicrange;
            if (features[i] < min ) features[i] = min;
        }        
    }
}

int main(int argc, char *argv[]) {
    FILE *fin, *fout;
    int stride = NB_FEATURES;
    int st = 0;
    int en = stride-1;
    float gain = 1.0;
    int frame_delay = 1;
    float pred = 0.0;
    int removemean = 0;
    float lower = -1E32;
    int writeall = 0;
    float dynamicrange = 100.0;
    int drlate = 0;
    int timestep = 1;
    int timeoffset = 0;
    int mean_l2 = 0;
    int pred2 = 0;
    
    static struct option long_options[] = {
        {"startcol",      required_argument, 0, 's'},
        {"endcol",        required_argument, 0, 'e'},
        {"stride",        required_argument, 0, 't'},
        {"timestep",      required_argument, 0, 'i'},
        {"timeoffset",    required_argument, 0, 'o'},
        {"gain",          required_argument, 0, 'g'},
        {"pred",          required_argument, 0, 'p'},
        {"delay",         required_argument, 0, 'd'},
        {"removemean",    no_argument,       0, 'm'},
        {"meanl2",        no_argument,       0, '2'},
        {"lower",         required_argument, 0, 'l'},
        {"dynamicrange",  required_argument, 0, 'y'},
        {"drlate",        no_argument,       0, 'z'},
        {"writeall",      no_argument,       0, 'w'},
        {"pred2",         no_argument,       0, 'r'},
        {0, 0, 0, 0}
    };

    int opt_index = 0;
    int c;
    
    while ((c = getopt_long (argc, argv, "s:e:t:g:p:d:ml:y:i:o:2zw", long_options, &opt_index)) != -1) {
        switch (c) {
        case 's':
            st = atoi(optarg);
            break;
        case 'e':
            en = atoi(optarg);
            break;
        case 't':
            stride = atoi(optarg);
            st = 0;
            en = stride-1;
            break;
        case 'i':
            timestep = atoi(optarg);
            break;
        case 'o':
            timeoffset = atoi(optarg);
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
        case '2':
            mean_l2 = 1;
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
        case 'z':
            drlate = 1;
            break;
        case 'r':
            pred2 = 1;
            break;
        default:
        helpmsg:
            fprintf(stderr, "VQ pre-processing tool\n\nusage: %s  -s startCol -e endCol [options] input.f32 output.f32\n"
                            "\n"
                            "-t K                   Vector length, e.g. 20\n"
                            "-g gain                Gain applied to vectors\n"
                            "-p predCoeff           Coefficient used for prediction, e.g. 0.9\n"
                            "-d predDelay           Delay (frames) between vectors used for prediction\n"
                            "--pred2                Output x[n-1]-(x[n]+x[n-3])/2, use with -d 4\n"
                            "--removemean           Remove mean from vectors\n"
                            "--lower minEnergy      Remove all vectors less than minEnergy\n"
                            "--writeall             Write all K outputs for each vector, even if EndCol-StartCol+1 < stride\n"
                            "--dynamicrange RangedB Clip min value of each vector to max - RangedB\n"
                            "--drlate               Dynamic range clip after mean removed\n"
                            "--timestep frames      time step (frames) between vectors\n"
                            "--timeoffset frames    Ignore this many initial frames\n"
                            "input.f32 output.f32\n", argv[0]);
            exit(1);
        }
    }
    if ( (argc - optind) < 2) {
        fprintf(stderr, "Too few arguments\n");
        goto helpmsg;
    }
 
    if (removemean && (pred != 0.0)) {
        /* These two options cause unexpected jumps in predictive VQ training data */ 
        fprintf(stderr, "Warning - low energy vector removal and prediction enabled\n");
    }
    
    fin = fopen(argv[optind],"rb"); assert(fin != NULL);
    fout = fopen(argv[optind+1],"wb"); assert(fout != NULL);
    printf("extracting from %d to %d inclusive (K=%d, timestep=%d) ... \n"
           "gain = %f pred = %f pred_frame_delay = %d dynamic range = %f drlate = %d\n",
           st, en, stride, timestep, gain, pred, frame_delay, dynamicrange, drlate);
   
    float features[stride], features_prev[frame_delay][stride], out[stride];
    int i,f,rd=0,wr=0;
    
    for(i=0; i<stride; i++)
        out[i] = 0.0;
    for (f=0; f<frame_delay; f++)
        for(i=0; i<stride; i++)
            features_prev[f][i] = 0.0;

    for(int i=0; i<timeoffset; i++)
        assert(fread(features, sizeof(float), stride, fin) == stride);

    while((fread(features, sizeof(float), stride, fin) == stride)) {
        rd++;

        if (drlate == 0)
          restrict_dynamic_range(features, st, en, dynamicrange);

	float mean = 0.0;
	if (mean_l2) {
            // if features are in dB, mean is L2 norm of linear vector
            for(i=st; i<=en; i++)
	        mean += pow(10,features[i]/10);
	    mean = 10*log10(mean/(en-st+1));
        } else {
            for(i=st; i<=en; i++)
	        mean += features[i];
	    mean /= (en-st+1);
        }
          
 	if (removemean) {
	    for(i=0; i<stride; i++)
		features[i] -= mean;
	}
        
	if (pred2) {
            /* prediction error from mean of past and future */
            assert(frame_delay == 4);
            for(i=st; i<=en; i++)
                out[i] = features_prev[1][i] - 0.5*(features[i] + features_prev[3][i]);
        }
        else {
            /* no prediction with pred == 0 */
            for(i=st; i<=en; i++)
                out[i] = gain*(features[i] - pred*features_prev[frame_delay-1][i]);
        }
        
        if (drlate)
          restrict_dynamic_range(features, st, en, dynamicrange);

	if (mean > lower) {
	    if (writeall)
                fwrite(out, sizeof(float), stride, fout);
            else
                fwrite(&out[st], sizeof(float), en-st+1, fout);
	    wr++;
	}
	for (f=frame_delay-1; f>0; f--)
	    for(i=0; i<stride; i++)
		features_prev[f][i] = features_prev[f-1][i];
	for(i=0; i<stride; i++)
	    features_prev[0][i] = features[i];

        for(int i=0; i<timestep-1; i++) {
            int ret = fread(features, sizeof(float), stride, fin);
            if (ret != stride) {
                fprintf(stderr, "warning: end of input reached ...\n");
            } else rd++;
        }
    }

    fclose(fin); fclose(fout);
    fprintf(stderr, "%d input %d extracted\n", rd, wr);
    return 0;
}

