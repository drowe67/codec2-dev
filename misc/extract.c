#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#define NB_FEATURES 55
#define NB_BANDS    18

int main(int argc, char *argv[]) {
    FILE *fin, *fout;

    fin = fopen(argv[1],"rb"); assert(fin != NULL);
    fout = fopen(argv[2],"wb"); assert(fout != NULL);
    int st = atoi(argv[3]); int en = atoi(argv[4]); float scale = atof(argv[5]);
    float pred = atof(argv[6]);
    printf("extracting from %d to %d inclusive... scale factor = %f pred = %f\n", st, en, scale, pred);
   
    float features[NB_FEATURES], features_prev[NB_FEATURES], delta[NB_FEATURES];
    int i;
    
    for(i=0; i<NB_FEATURES; i++)
        features_prev[i] = 0.0;

    while((fread(features, sizeof(float), NB_FEATURES, fin) == NB_FEATURES)) {
        for(i=st; i<=en; i++) {
            delta[i] = scale*(features[i] - pred*features_prev[i]);
            features_prev[i] = features[i];
        }
        fwrite(&delta[st], sizeof(float), en-st+1, fout);
    }

    fclose(fin); fclose(fout);
    
    return 0;
}

