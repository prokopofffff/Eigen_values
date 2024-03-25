#include "funcs.h"
#define EPSILON 1e-16
#define MAX_ITERS 300000

double norm(double* A, int n){
	int i;
	int j;
	double sum;
	double res = 0;

	for (i = 0; i < n; ++i){
		sum = 0;
		for (j = 0; j < n; ++j)
			sum += fabs(A[i * n + j]);

		if (res < sum)
			res = sum;
	}

	return res;
}

void three_diag(double* A, int n){
	int i;
	int j;
	int k;
	double sum1;
	double sum2;

	for (i = 0; i < n - 2; ++i){
		sum1 = 0;
		for (j = i + 2; j < n; ++j) 
            sum1 += A[j * n + i] * A[j * n + i];

		sum2 = sqrt(A[(i + 1) * n + i] * A[(i + 1) * n + i] + sum1);

		if (sum2 < EPSILON){
			A[(i + 1) * n + i] = 0;
			A[(i + 2) * n + i] = 0;
			continue;
		}
		if (sum1 < EPSILON){
			A[(i + 2) * n + i] = 0;
			continue;
		}

		A[(i + 1) * n + i] -= sum2;

		sum1 = 1. / sqrt(A[(i + 1) * n + i] * A[(i + 1) * n + i] + sum1);
		for (j = i + 1; j < n; ++j) 
            A[j * n + i] *= sum1;

		for (j = i + 1; j < n; ++j){
			sum1 = 0;
			for (k = i + 1; k < n; k++) 
                sum1 += A[j * n + k] * A[k * n + i];
			A[i * n + j] = sum1;
		}

		sum1 = 0;
		for (j = i + 1; j < n; ++j) 
            sum1 += A[i * n + j] * A[j * n + i];
		sum1 *= 2;

		for (j = i + 1; j < n; ++j) 
            A[i * n + j] = 2 * A[i * n + j] - sum1 * A[j * n + i];

		for (j = i + 1; j < n; ++j)
			for (k = i + 1; k < n; ++k) 
                A[j * n + k] -= A[i * n + j] * A[k * n + i] + A[i * n + k] * A[j * n + i];

		A[(i + 1) * n + i] = sum2;
		A[i * n + i + 1] = sum2;
		for (j = i + 2; j < n; ++j){
			A[j * n + i] = 0;
			A[i * n + j] = 0;
		}
	}
}

int n_(double *A, int n, double lambda, double eps){
    if (n <= 0) return 0;
    int i;
    double up;
    double mid;
    double down;
    int s = 1;
    int count = 0;

    mid = A[0] - lambda;
    if (fabs(mid) < eps) mid = eps;
    if (s != s * (mid > 0 ? 1 : -1)) count++;
    up = A[1] / mid;
    s = s * (mid > 0 ? 1 : -1);

    for (i = 1; i < n - 1; i++){
        down = A[i * n + i - 1];
        mid = A[i * n + i] - lambda - down * up;
        if (fabs(mid) < eps) mid = eps;
        if (s != s * (mid > 0 ? 1 : -1)) count++;
        s = s * (mid > 0 ? 1 : -1);
        up = A[i * n + i + 1] / mid;
    }
    down = A[n * n - 2];
    mid = A[n * n - 1] - lambda - down * up;
    if (fabs(mid) < eps) mid = eps;

    if (s != s * (mid > 0 ? 1 : -1)) count++;
    
    return count;
}

int eigen_values(double * A, double * vals, int n, double eps){
    double a = -norm(A, n) * (1 + eps), b = -a;
    int iters = 0, count = 0;
    eps *= b;
    eigen_values_int(A, vals, n, iters, count,  eps, a, b);
    if (count != n) return -1;
    return iters;
}

void eigen_values_int(double *A, double *vals,  int n, int &iters, int &count, double eps, double a, double b){
    iters++;
    if (iters < MAX_ITERS){
		int k = n_(A, n, b, eps) - n_(A, n, a, eps);
		if (k != 0 && b - a > eps){
			eigen_values_int(A, vals, n, iters, count, eps, a, (a + b) / 2);
			eigen_values_int(A, vals, n, iters, count, eps, (a + b) / 2, b);
		}
		else if(k != 0){
			for (int i = count; i < count + k; i++){
				vals[i] =  (a + b) / 2;
			}
			count += k;
		}
	}
}

double get_r1(double* A, double* vals, int n){
    double sum1 = 0;
    double sum2 = 0;

    for(int i = 0; i < n; i++){
        sum1 += A[i * n + i];
        sum2 += vals[i];
    }

    return fabs(sum1 - sum2) / norm(A, n);
}

double get_r2(double* A, double* vals, int n){
    double sum1 = 0;
    double sum2 = 0;

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            sum1 += A[i * n + j] * A[j * n + i];
        }
        sum2 += vals[i] * vals[i];
    }

    return fabs(sqrt(sum1) - sqrt(sum2)) / norm(A, n);
}