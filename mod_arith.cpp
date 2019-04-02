#include <stdio.h>
#include <iostream>
#include <memory>
#include <vector>
#include <cmath>

#include "mod_arith.h"


//using namespace std;




template<class Type>
Type mcd(Type a, Type b){
    if (b == 0){
	return a;
    }
    else
	return mcd(b, a % b);
}

// Algoritmo Euclides extendido:


 






// Inverso multiplicativo
//binary better??


//Garner's Algorithm: Chinesse Remainder Theorem
//template<class Type>
	    



