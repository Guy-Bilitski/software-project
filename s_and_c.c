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
S_and_C *create_empty_S_and_C() {
    S_and_C *s_and_c = (S_and_C *)malloc(sizeof(S_and_C));
    if (s_and_c == NULL){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    return s_and_c;
}

/* Getters */
double s_and_c_get_s(S_and_C *s_and_c) {
    return s_and_c->s;
}

double s_and_c_get_c(S_and_C *s_and_c) {
    return s_and_c->c;
}

/* Setters */
void S_and_C_set_values(S_and_C *s_and_c, double s, double c) {
    s_and_c->s = s;
    s_and_c->c = c;
}
