#ifndef MATRIX
#define MATRIX

#include <algorithm>
#include <fstream>
#include <iostream>
#include <typeinfo>
#include <iomanip>

void matrix_input(char* filename, double* matrix, int n);

void matrix_output(double *matrix, int l, int n, int m);

double f(int k, int n, int i, int j);

#endif