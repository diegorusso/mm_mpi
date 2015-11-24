#include <math.h>
#include <float.h>
#include <stdio.h>

// Function 'check' checks the correctness for the matrix-matrix
// multiplication in function 'mxm'. Because of the particular
// data that are generatedin function 'gendat', this check can be
// used.
int check(int m, int l, int n, double *c) {
    double  tvalue = (double)l;
    int     i, j, ok = 0;
    int     ELEM = (int)(m * n);
    // DBL_EPSILON the min positive number such that 1.0 + DBL_EPSILON != 1.0
    // eps is a number infinitely small
    double  eps = DBL_EPSILON;

    // The data I need to check are in array, so I iterate over every element
    // in c to check its value
    for (i = 0; i < ELEM; i++) {
        // fabs returns the absolute value
        // if the absolute value of the difference between tvalue and the
        // value of the cell is bigger than eps, it means there is
        // something wrong in the check.
        // tvalue - c[i] shouldn't be bigger than the minimum positive
        // number in the system. Ideally the difference should be 0.
        if (fabs(tvalue - c[i]) > eps)
            ok++;
    }
    // To be successfull, ok should be 0
    return (ok);
}
