#include <iostream>
#include <iomanip>
#include <vector>

#include "matrix.h"


void print(Matrix & M, Matrix & N){
    for(int i = 0; i < M.num_rows; i++){
	std::cout<<"(";
	for(int j = 0; j < M.num_cols; j++){
	    if(j == M.num_cols -1){
		std::cout<<std::setw(1)<<std::setfill(' ')<<M[i][j];
	    }
	    else
		std::cout<<std::setw(1)<<std::setfill(' ')<<M[i][j]<<"  ";
	}
	
	std::cout<<")"<<std::setw(7)<<std::setfill(' ');
	
	std::cout<<"(";

	for(int j = 0; j < N.num_cols; j++){
	    if(j == N.num_cols -1){
		std::cout<<std::setw(1)<<std::setfill(' ')<<N[i][j];
	    }
	    else
		std::cout<<std::setw(1)<<std::setfill(' ')<<N[i][j]<<"  ";
	}
		
	std::cout<<")"<<std::endl<<std::endl;
    }
}

Matrix::Matrix(){
    num_rows = 0;
    num_cols = 0;
    rows = std::vector<Int>(0);
    mat = std::vector<std::vector<Int>>();
}

Matrix::Matrix(int row, int col){
    num_rows = row;
    num_cols = col;
    rows = std::vector<Int>(col, Int(0));
    mat = std::vector<std::vector<Int>>(row, rows);
}

Matrix::Matrix(int row, int col, int mod){
    num_rows = row;
    num_cols = col;
    rows = std::vector<Int>(col, Int(0, mod));
    mat = std::vector<std::vector<Int>>(row, rows);
}

void Matrix::swap(int row1, int row2){
    mat[row1].swap(mat[row2]);
}

Matrix& Matrix::transpose(){
    auto result = new Matrix(num_rows, num_cols);
    for(int i = 0; i < num_rows; i++){
	for(int j = 0; j < num_cols; j++){
	    (*result)[i][j] = mat[j][i];
	}
    }
    return *result;
}

Matrix& Matrix::operator+(Matrix& M){
    auto result = new Matrix(num_rows, num_cols);
    for(int i = 0; i < num_rows; i++){
	for(int j = 0; j < num_cols; j++){
	    (*result)[i][j] = (*this)[i][j] + M[i][j];
	}
    }
    return *result;
}

Matrix& Matrix::operator-(Matrix& M){
    auto result = new Matrix(num_rows, num_cols);
    for(int i = 0; i < num_rows; i++){
	for(int j = 0; j < num_cols; j++){
	    (*result)[i][j] = (*this)[i][j] - M[i][j];
	}
    }
    return *result;    
}

std::vector<Int> mult(Int a, std::vector<Int> v){
    std::vector<Int> result;
    for(int i = 0; i < v.size(); i++){
	Int aux = a * v[i];
	result.push_back(aux);
    }
    return result;
}

std::vector<Int> sub(std::vector<Int> u, std::vector<Int> v){
    std::vector<Int> result;
    for(int i = 0; i < v.size(); i++){
	Int aux = u[i] - v[i];
	result.push_back(aux);
    }
    return result;
}

Matrix& Matrix::gaussianElim(){
    auto result = new Matrix(num_rows, num_cols);
    *result = *this;
    int prim = ((*this)[0][0]).getMod();

    auto X = new Matrix;
    *X = MatrixIn_mod(num_rows, prim);

    int i = 0;
    int j = 0;
    
    
    while(i < num_rows){
    	bool T = false;
    	while((j < num_cols) && !T){
    	    if((*result)[i][j].getValue() != 0){
    		T = true;
    	    }
	    
    	    else{
    		int max_row = i;
    		int max_value = 0;
    		for(int k = i+1; k < num_rows; k++){
		    
    		    if((*result)[k][j].getValue() > max_value){
    			max_row = k;
    			max_value = (*result)[k][j].getValue();
    		    }
    		}
		
    		if(max_row != i){
		    (*X).swap(max_row,i);
		    (*result).swap(max_row,i);    
    		    T = true;
    		}

    		else{
    		    j++;
    		}
    	    }
    	}
	
    	if(T){
    	    Int d = Int(invMult((*result)[i][j].getValue(), prim), prim);

	    (*X)[i] = mult(d, (*X)[i]);
	    (*result)[i] = mult(d, (*result)[i]);

	    for(int s = 0; s < num_rows; s++){
		
		if(i != s){
		    (*X)[s] = sub((*X)[s], mult((*result)[s][j], (*X)[i]));
		    (*result)[s] = sub((*result)[s], mult((*result)[s][j], (*result)[i]));
		}
	    }
    	}
    	i++;
    	j++;
    }
    // std::cout<<"result"<<std::endl;
    // std::cout<<(*result)<<std::endl;
    // std::cout<<"X"<<std::endl;
    // std::cout<<(*X)<<std::endl;
    return *result;
}

std::vector<std::vector<Int>> Matrix::kernel_basis(){
    std::vector<std::vector<Int>> basis;
    Matrix result = *this;
//    std::cout<<(*this)<<std::endl;
    int prim = ((*this)[0][0]).getMod();

    Matrix X = MatrixIn_mod(num_rows, prim);

    int i = 0;
    int j = 0;
    int dim_null_space;
    int rank_count = 0;
    
    while(i < num_rows){
    	bool T = false;
    	while((j < num_cols) && !T){
    	    if(result[i][j].getValue() != 0){
    		T = true;
    	    }
	    
    	    else{
    		int max_row = i;
    		int max_value = 0;
    		for(int k = i+1; k < num_rows; k++){
		    
    		    if(result[k][j].getValue() > max_value){
    			max_row = k;
    			max_value = result[k][j].getValue();
    		    }
    		}
		
    		if(max_row != i){
		    X.swap(max_row,i);
		    result.swap(max_row,i);    
    		    T = true;
    		}

    		else{
    		    j++;
    		}
    	    }
    	}
	
    	if(T){
//	    std::cout<<"pivot in row = "<<i<<std::endl;
	    rank_count++;
	    int inv_no_modulus = invMult(result[i][j].getValue(), prim);
    	    Int d = Int(inv_no_modulus, prim);

	    std::vector<Int> auxXi = X[i];
	    X[i] = mult(d, auxXi);
	    std::vector<Int> auxRESULT = result[i];
	    result[i] = mult(d, auxRESULT);

	    for(int s = 0; s < num_rows; s++){
		
		if(i != s){
		    std::vector<Int> mult_aux_1 = mult(result[s][j], X[i]);
		    X[s] = sub(X[s], mult_aux_1);
		    std::vector<Int> mult_aux_2 = mult(result[s][j], result[i]);
		    result[s] = sub(result[s], mult_aux_2);
		}
	    }
    	}
    	i++;
    	j++;
    }

    dim_null_space = num_rows - rank_count;
//    std::cout<<"dim null space: "<<dim_null_space<<std::endl;
    
    for(int k = 0; k < dim_null_space; k++){
	basis.push_back(X[num_rows-1-k]);
//	std::cout<<"added to bassis"<<std::endl;
    }

    // for(auto&& v : basis){
    // 	for(auto&& elem: v){
    // 	    std::cout<<elem<<", ";
    // 	}
    // 	std::cout<<std::endl;
    // }
//    std::cout<<X<<std::endl;
    print(X, result);
    return basis;
}

Matrix& MatrixIn(int n){
    auto result = new Matrix(n, n);
    for(int i = 0; i < n; i++){
	for(int j = 0; j < n; j++){
	    if(i == j){
		(*result)[i][j] = Int(1);
	    }
	    else
		(*result)[i][j] = Int(0);
	}
    }
    return *result;
}

Matrix& MatrixIn_mod(int n, int mod){
    auto result = new Matrix(n, n);
    for(int i = 0; i < n; i++){
	for(int j = 0; j < n; j++){
	    if(i == j){
		(*result)[i][j] = Int(1, mod);
	    }
	    else
		(*result)[i][j] = Int(0, mod);
	}
    }
    return *result;
}



std::ostream & operator<<(std::ostream & os, Matrix & M){
    for(int i = 0; i < M.num_rows; i++){
	os<<"(";
	for(int j = 0; j < M.num_cols; j++){
	    if(j == M.num_cols -1){
		os<<std::setw(1)<<std::setfill(' ')<<M[i][j];
	    }
	    else
		os<<std::setw(1)<<std::setfill(' ')<<M[i][j]<<"  ";
	}
	os<<")"<<std::endl<<std::endl;
    }
    return os;
}


