/* compare ints - a test utility */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <errno.h>

/*  Declarations */

/* Globals */

/* Functions */
int get_data(FILE *f, uint32_t *d, int bytes) {
    // TODO Loop on reads until, but catch EOF!!
    int read = fread(d, bytes, 1, f);
    if (read != 1) {
        return(0);
        }
    else return(1);
    }


/* Main */

int main(int argc, char *argv[]) {

    char usage[] = "Usage: %s [-s size_in_bytes] [-t tolerance] file1 file2\n";

    int bytes = 1;
    int tol = 1;

    int opt;
    while ((opt = getopt(argc, argv, "s:t:")) != -1) {
        switch (opt) {
            case 's':
                bytes = atoi(optarg);
                break;
            case 't':
                tol = atof(optarg);
                break;
            default:
                fprintf(stderr, usage, argv[0]);
                exit(1);
            }
        }

    if ((optind + 2) > argc) {
        fprintf(stderr, usage, argv[0]);
        exit(1);
        }
    char *fname1 = argv[optind++];
    char *fname2 = argv[optind++];

    FILE *f1 = fopen(fname1, "rb");
    if (f1 == NULL) {
        fprintf(stderr, "Error opening file1 \"%s\": ", fname1);
        perror(NULL);
        exit(1);
        }

    FILE *f2 = fopen(fname2, "rb");
    if (f2 == NULL) {
        fprintf(stderr, "Error opening file2 \"%s\": ", fname2);
        perror(NULL);
        exit(1);
        }

    uint32_t data1, data2;

    int count = 0;
    int errors = 0;
    int rms_sum = 0;

    while (get_data(f1, &data1, bytes)) {
        if (!get_data(f2, &data2, bytes)) {
            fprintf(stderr, "Error: file2 is shorter!");
            exit(1);
            }
        switch (bytes) {
            case 1: {
                uint8_t err = abs((uint8_t)data1 - (uint8_t)data2);
                if (err > tol) {
                    errors ++;
                    printf("%d %2d %2d\n", count, data1, data2);
                    }
                rms_sum += (err * err);
                count ++;
                } break;
            case 2: {
                uint16_t err = abs((uint16_t)data1 - (uint16_t)data2);
                if (err > tol) {
                    errors ++;
                    printf("%d %2d %2d\n", count, data1, data2);
                    }
                rms_sum += (err * err);
                count ++;
                } break;
            default: 
                fprintf(stderr, "Error: unsupported size %d bytes\n", bytes);
                exit(1);
            }
        }

    if (errors) {
        printf("Fail: %d errors\n", errors);
        printf("      rms error = %f\n", ((double)rms_sum/count));
        exit(1);
        }
    else printf("Pass\n");
    exit(0);

    } // main


/* vi:set ts=4 et sts=4: */
