#include <stdlib.h>
#include <sys/types.h>
#include <sys/time.h>


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
