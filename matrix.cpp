#include "matrix.h"

void matrix_input(char* filename, double* matrix, int n){
    std::ifstream file(filename);

    if(!(file.is_open())){
        delete[] matrix;
        throw "Can not open file";
    }

    int count = 0;

    while(file){
        if(!(file >> matrix[count])){
            delete[] matrix;
            throw "Incorrect type";
        }
        count++;
        if(count == n * n) break;
    }

    if(count < n * n){
        delete[] matrix;
        throw "Not enough data";
    }

    file.close();
}

double f(int k, int n, int i, int j){
    if (k == 1)
        return (n - (i > j ? i : j) + 1);
    if (k == 2)
        {
            if (i == j)
                return 2;
            if (i - j == 1 || j - i == 1)
                return -1;
            return 0;
        }
    if (k == 3)
        {
            if (i == j && j < n)
                return 1;
            if (j == n)
                return i;
            if (i == n)
                return j;
            return 0;
        }
    if (k == 4)
        return (double) (1. / (i + j - 1));
    return 0;
}

void matrix_output(double *matrix, int l, int n, int m){
    for(int i = 0; i < l && i < m; i++){
        for(int j = 0; j < n && j < m; j++){
            std::cout << std::scientific << matrix[i * n + j] << " ";
        }
        std::cout << std::endl;
    }
}