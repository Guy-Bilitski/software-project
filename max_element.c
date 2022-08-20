#include <stdio.h>
#include <stdlib.h>

typedef struct MaxElement
{
    int i;
    int j;
    double value;
} MaxElement;

/* MaxElemnt API */
MaxElement *create_max_element(double value, int i, int j);

/* getters */
double max_element_get_value(MaxElement *max_element);
int max_element_get_index1(MaxElement *max_element);
int max_element_get_index2(MaxElement *max_element);

/* setters */
void max_element_set_new_values(MaxElement *max_element, double value, int i, int j);
void max_element_set_value(MaxElement *max_element, double value);
void max_element_set_index1(MaxElement *max_element, int i);
void max_element_set_index2(MaxElement *max_element, int j);

/* debugging */
void print_max_element(MaxElement *max_element);


/* MaxElemnt API */
MaxElement *create_max_element(double value, int i, int j) {
    MaxElement *max_element = (MaxElement *)malloc(sizeof(MaxElement));
    max_element->value = value;
    max_element->i = i;
    max_element->j = j;
    return max_element;
}

/* getters */
double max_element_get_value(MaxElement *max_element) {
    return max_element->value;
}

int max_element_get_index1(MaxElement *max_element) {
    return max_element->i;
}

int max_element_get_index2(MaxElement *max_element) {
    return max_element->j;
}

/* setters */
void max_element_set_new_values(MaxElement *max_element, double value, int i, int j) {
    max_element_set_index1(max_element, i);
    max_element_set_index2(max_element, j);
    max_element_set_value(max_element, value);
}

void max_element_set_value(MaxElement *max_element, double value) {
    max_element->value = value;
}

void max_element_set_index1(MaxElement *max_element, int i) {
    max_element->i = i;
}

void max_element_set_index2(MaxElement *max_element, int j) {
    max_element->j = j;
}

/* debugging */
void print_max_element(MaxElement *max_element) { 
    printf("max element value: %f ", max_element_get_value(max_element));
    printf("max element i: %d ", max_element_get_index1(max_element));
    printf("max element j: %d", max_element_get_index2(max_element));
}