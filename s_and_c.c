#include "spkmeans.h"
#ifndef S_AND_C_IS_DEFINED
#define S_AND_C_IS_DEFINED
typedef struct S_and_C
{
    double s;
    double c;
} S_and_C;
#endif

/* S_and_C API */
double s_and_c_get_s(S_and_C s_and_c) {
    return s_and_c.s;
}

double s_and_c_get_c(S_and_C s_and_c) {
    return s_and_c.c;
}