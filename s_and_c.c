typedef struct S_and_C
{
    double s;
    double c;
} S_and_C;


/* S_and_C API */
double s_and_c_get_s(S_and_C s_and_c);
double s_and_c_get_c(S_and_C s_and_c);

/* S_and_C API */
double s_and_c_get_s(S_and_C s_and_c) {
    return s_and_c.s;
}

double s_and_c_get_c(S_and_C s_and_c) {
    return s_and_c.c;
}