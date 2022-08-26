#include "spkmeans.h"
#include <stdlib.h>
#ifndef S_AND_C_IS_DEFINED
#define S_AND_C_IS_DEFINED
typedef struct S_and_C
{
    double s;
    double c;
} S_and_C;
#endif

/* S_and_C API */
S_and_C *create_S_and_C(double s, double c) {
    S_and_C *s_and_c = (S_and_C *)malloc(sizeof(S_and_C));
    s_and_c->s = s;
    s_and_c->c = c;
    return s_and_c;
}

double s_and_c_get_s(S_and_C *s_and_c) {
    return s_and_c->s;
}

double s_and_c_get_c(S_and_C *s_and_c) {
    return s_and_c->c;
}

void print_s_and_c(S_and_C *s_and_c) {
    space();
    printf("c is: %f s is: %f", s_and_c_get_c(s_and_c), s_and_c_get_s(s_and_c));
    space();
}