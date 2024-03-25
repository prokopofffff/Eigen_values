#include <iostream>
#include "matrix.h"
#include "funcs.h"

int main(int argc, char* argv[]){
    try{
        if(argc != 5 && argc != 6){
            throw "Incorrect input";
            return -1;
        }

        int n = atoi(argv[1]), m = atoi(argv[2]), k = atoi(argv[4]);
        double eps = atof(argv[3]);
        char* filename;
        double *A = new double[n * n];

        if(k == 0){
            filename = argv[5];
            matrix_input(filename, A, n);
        }
        else{
            for(int i = 0; i < n; i++){
                for(int j = 0; j < n; j++){
                    A[i * n + j] = f(k, n, i + 1, j + 1);
                    if(k != 1 && k != 2 && k != 3 && k != 4){
                        delete[] A;
                        throw "Undefined behaivor";
                        return -1;
                    }
                }
            }
        }

        double *vals = new double[n];

        for(int i = 0; i < n; i++){
            vals[i] = 0;
        }

        matrix_output(A, n, n, m);

        double time1 = clock();
        three_diag(A, n);
        time1 = (clock() - time1) / CLOCKS_PER_SEC;

        double time2 = clock();
        int iters = eigen_values(A, vals, n , eps);
        time2 = (clock() - time2) / CLOCKS_PER_SEC;

        std::cout << std::endl;

        matrix_output(vals, 1, n, m);

        double r1 = get_r1(A, vals, n);
        double r2 = get_r2(A, vals, n);

        printf ("%s : Residual1 = %e Residual2 = %e Iterations = %d Iterations1 = %d Elapsed1 = %.2f Elapsed2 = %.2f\n",
                argv[0], r1, r2, iters, iters / n, time1, time2);

        delete[] A;
        delete[] vals;

    }
    catch(const char *err){
        std::cerr << err;
    }
}