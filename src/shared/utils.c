#include "shared.h"
#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/time.h>


void print_usage(char *app_name) {
    printf("Usage: %s [-h -v] -f filename\n", app_name);
}


void parse_arguments(int argc, char **argv, FILE **input_file){
    int option = 0;
    int file_specified = 0;

    while ((option = getopt(argc, argv,"hvf:")) != -1) {
        switch (option) {
            case 'f' : {
                *input_file = fopen(optarg, "r");  // I need to return this
                file_specified = 1;
            } break;
            case 'v' : {
                verbose = 1;
            } break;
            case 'h' : {
                print_usage(argv[0]);
                exit(0);
            } break;
            default:
                break;
        }
    }

    // Do some checks
    if (file_specified == 0) {
        printf("ERROR: Please specify the input file with -f\n");
        print_usage(argv[0]);
        exit(1);
    };

    if (input_file == NULL) {
        printf("ERROR: File doesn not exist\n");
        exit(1);
    };
};


// This function returns the wall clock time with micro seconds accuracy.
// The data type of the returned value is "double".
double cclock(void) {
    const double    micro = 1.0e-06;    // Conversion constant
    struct timeval  tp;                 // Structure used by gettimeofday
    double          wall_time;          // To hold the result
    static long     start = 0L, startu;

    // gettimeofday gives the number of seconds and microseconds since the
    // Epoch
    if (gettimeofday(&tp, NULL) == -1)
        wall_time = -1.0e0;
    else if (!start) {
        start = tp.tv_sec;          // seconds
        startu = tp.tv_usec;        // microseconds
        wall_time = 0.0e0;          // reset wall_time
    } else {
        // Let's calculate the wall_time with micro seconds accuracy
        wall_time = (double)(tp.tv_sec - start) + micro * (tp.tv_usec - startu);
    }
    return wall_time;
}
