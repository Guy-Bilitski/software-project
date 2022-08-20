typedef struct MaxElement
{
    int i;
    int j;
    double value;
} MaxElement;

/* MaxElemnt API TODO: move to new file*/
int max_element_get_value(MaxElement max_element);
int max_element_get_index1(MaxElement max_element);
int max_element_get_index2(MaxElement max_element);

/* MaxElemnt API */
int max_element_get_value(MaxElement max_element) {
    return max_element.value;
}

int max_element_get_index1(MaxElement max_element) {
    return max_element.i;
}

int max_element_get_index2(MaxElement max_element) {
    return max_element.j;
}