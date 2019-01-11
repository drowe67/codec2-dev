#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#define NB_FEATURES 55
#define NB_BANDS    18

int main(int argc, char *argv[]) {
    FILE *fin, *fout;

    fin = fopen(argv[1],"rb"); assert(fin != NULL);
    fout = fopen(argv[2],"wb"); assert(fout != NULL);
    int st = atoi(argv[3]); int en = atoi(argv[4]);
    printf("extracting from %d to %d...\n", st, en);
    
    float features[NB_FEATURES];

    while((fread(features, sizeof(float), NB_FEATURES, fin) == NB_FEATURES))
        fwrite(&features[st], sizeof(float), en-st+1, fout);

    fclose(fin); fclose(fout);
    
    return 0;
}

