#ifndef FUNCS
#define FUNCS

#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double norm(double* A, int n);

void three_diag(double* A, int n);

int n_(double* A, int n, double lambda, double eps);

int eigen_values(double* A, double* vals, int n, double eps);

void eigen_values_int(double * A, double * B,  int n, int &iters, int &count, double eps, double a, double b);

double get_r1(double* A, double* vals, int n);

double get_r2(double* A, double* vals, int n);

#endif