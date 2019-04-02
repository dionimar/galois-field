#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>
#include "int.h"
#include "polynom.h"


class Matrix{
public:
    int num_rows, num_cols;
    std::vector<Int> rows;
    std::vector<std::vector<Int>> mat;

    Matrix();
    
    explicit Matrix(int row, int col);

    explicit Matrix(int row, int col, int mod);

    ~Matrix() = default;
    
    Matrix(const Matrix & M) = default;
    Matrix & operator=(const Matrix & M) = default;

    inline Int & operator()(int row, int col){
	return mat[row][col];
    }

    inline std::vector<Int> & operator[](int j){
    	return mat[j];
    }

    void swap(int row1, int row2);
    
    Matrix& transpose();

    Matrix& operator+(Matrix& M);

    Matrix& operator-(Matrix& M);

    Matrix& gaussianElim();

    std::vector<std::vector<Int>> kernel_basis();

//    std::vector<std::vector<Int>> kernel_basis2();
    
};

Matrix& MatrixIn(int n);

Matrix& MatrixIn_mod(int n, int mod);

std::ostream & operator<<(std::ostream& os, Matrix & M);



#endif
