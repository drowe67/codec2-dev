/*

  extract.c
  david Rowe Jan 2019

  Extracts sub sets of vectors from .f32 files, used for LPCNet VQ experiments.
*/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#define NB_FEATURES 55 /* number of cols per row */

int main(int argc, char *argv[]) {
    FILE *fin, *fout;

    if (argc < 7) {
        fprintf(stderr, "usage: %s input.f32 output.f32 startCol endCol scaleFactor predCoeff [framesDelay]\n", argv[0]);
        exit(1);
    }
    fin = fopen(argv[1],"rb"); assert(fin != NULL);
    fout = fopen(argv[2],"wb"); assert(fout != NULL);
    int st = atoi(argv[3]); int en = atoi(argv[4]); float scale = atof(argv[5]);
    float pred = atof(argv[6]);
    int frame_delay;
    if (argc == 8)
        frame_delay = atoi(argv[7]);
    else
        frame_delay = 1;
    printf("extracting from %d to %d inclusive... scale factor = %f pred = %f frame_delay = %d\n", st, en, scale, pred, frame_delay);
   
    float features[NB_FEATURES], features_prev[frame_delay][NB_FEATURES], delta[NB_FEATURES];
    int i,f;
    
    for (f=0; f<frame_delay; f++)
        for(i=0; i<NB_FEATURES; i++)
            features_prev[f][i] = 0.0;

    while((fread(features, sizeof(float), NB_FEATURES, fin) == NB_FEATURES)) {
        for(i=st; i<=en; i++) {
            delta[i] = scale*(features[i] - pred*features_prev[frame_delay-1][i]);
        }
        fwrite(&delta[st], sizeof(float), en-st+1, fout);
        for (f=frame_delay-1; f>0; f--)
            for(i=0; i<NB_FEATURES; i++)
                features_prev[f][i] = features_prev[f-1][i];
        for(i=0; i<NB_FEATURES; i++)
            features_prev[0][i] = features[i];        
    }

    fclose(fin); fclose(fout);
    
    return 0;
}

