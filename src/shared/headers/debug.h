#include <stdio.h>

#ifdef DEBUG
    #define DEBUG 1     // Enable debug
#else
    #define DEBUG 0     // Disable debug
#endif

// debug_print function: the do { ... } while (0) idiom ensures that the code
// acts like a statement (function call)
#define debug_print(fmt, ...) \
            do { if (DEBUG) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, \
                                    __LINE__, __func__, __VA_ARGS__); } \
            while (0)

