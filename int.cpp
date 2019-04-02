#include <stdio.h>
#include <iostream>
#include <vector>
#include <stack>
#include <random>
#include <utility>
#include <assert.h>

#include "int.h"


std::default_random_engine & global_urng(){
    static std::default_random_engine u{};
    return u;
}

int random_number(int low, int high){
    static std::uniform_int_distribution<> d{};
    using parm_t = decltype(d)::param_type;
    return d(global_urng(), parm_t{low, high});
}

bool Int::check(){
    return (mod >= 0);
}

Int::Int(): val(0), mod(0){}

Int::Int(long long signed int m): val(m), mod(0){
    if(!check()) throw Invalid();
}

Int::Int(long long signed int m, long long signed int p): val(m), mod(p){
    if(!check()){
	throw Invalid();
    }
    
    if(p != 0){
	val = ((m % p)<0) ? ((m%p)+p) : (m%p);
    }
    
    else{
	val = m;
    }
    mod = p;
}

// Return opposite of Int
void Int::setNeg(){
    // if(*this.coef.val>0)
    *this = Int(-val,mod);
}

// copy constructor
Int::Int(const Int & a){
    val = a.val;
    mod = a.mod;
}

// assignment operator = default
Int & Int::operator=(const Int & a)=default;

// set functions
void Int::set_value(const long long signed int a){
    val = a;
    while(val >= mod){
	val -= mod;
    }
    while(val < 0){
	val += mod;
    }
}

void Int::set_mod(const long long signed int p){
    mod = p;
    val = val % mod;
}

bool Int::operator==(const Int & a) const{
    return ((val == a.val)&&(mod == a.mod));
}

bool Int::operator!=(const Int & a) const{
    return (val != a.val);
}

// simplify the code by ussing code like if x == 0
bool Int::operator==(long long signed int a) const{
    return (val == a);
}

bool Int::operator!=(long long signed int a) const{
    return (val != a);
}

Int Int::operator+(const Int & a) const{
    Int result;
    long long signed int aux = this->val + a.val;
    if(a.mod != 0){
	result.val = ((aux % a.mod)<0) ? ((aux % a.mod) + a.mod) : (aux % a.mod);
	result.mod = a.mod;
    }
    else{
	result.val = aux;
	result.mod = 0;
    }
    return result;      
}

Int Int::operator-(const Int & a) const{
    Int result;
    if(mod != 0){
	result.val = (((val-a.val)%mod)<0) ? (((val-a.val)%mod)+mod) : ((val-a.val)%mod) ;
	result.mod = mod;
    }
    else{
	result.val = val - a.val;
	result.mod = 0;
    }
    return result;
}

Int Int::operator*(const Int & a)const {
    Int result;
    long long signed int aux = val * a.val;
    if(mod != 0){
	result.val = ((aux % mod)<0) ? ((aux % mod) + mod) : (aux % mod);
	result.mod = mod;
    }
    else{
	result.val = aux;
	result.mod = 0;
    }
    return result;
}

Int Int::operator*(long long signed int a)const {
    Int result = Int(a*(this->val), this->mod);
    return result;
}

// Procedural functions
void sum(const Int & a, const Int & b, Int & result){
    assert(a.getMod() == b.getMod());
    result.set_mod(a.getMod());
    result.set_value(a.getValue() + b.getValue());
}

void sub(const Int & a, const Int & b, Int & result){
    assert(a.getMod() == b.getMod());
    result.set_mod(a.getMod());
    result.set_value(a.getValue() - b.getValue());    
}

void mult(const Int & a, const Int & b, Int & result){
    assert(a.getMod() == b.getMod());
    result.set_mod(a.getMod());
    result.set_value(a.getValue() * b.getValue());
}

// void power(const Int & a, const long long signed int n, Int & result){
//     std::vector<int> bin_n = intToBin(n);
//     int k = bin_n.size();
//     std::vector<Int> b;
//     b.push_back(*this);
//     Int aux;
    
//     for(int i = 1; i < k; i++){
// 	aux = b[i-1]*b[i-1];
// 	if(bin_n[i] == 1){	   
// 	    b.push_back(aux*(*this));
// 	}
	
// 	else{
// 	    b.push_back(aux);
// 	}
//     }
//     Int b1 = Int(b[k-1]);    
// }

void inverse_mod(const Int & a, Int & result){
    result.set_mod(a.getMod());
    
    result.set_value(invMult(a.getValue(), a.getMod()));
}

void div(const Int & a, const Int & b, Int & result){
    assert(a.getMod() == b.getMod());
    result.set_mod(a.getMod());
    inverse_mod(b, result);
    mult(a, result, result);
}

// Extended Euclides Algorithm for integers
void extEuclides(long long signed int a, long long signed int b, long long signed int auxA, long long signed int auxB, long long signed int mcd){
    long long signed int q, r;
    long long signed int x1, x2, y1, y2;
    
    if(b == 0){
        auxA = 1;
        auxB = 0;
        mcd = a;
    }
    
    else{
        x1 = 0;
        x2 = 1;
        y1 = 1;
        y2 = 0;

        while(b>0){
            q = a/b;
            r = a%b;
            a = b;
            b = r;
            auxA = x2 - (q*x1);
            auxB = y2 - (q*y1);
            x2 = x1;
            x1 = auxA;
            y2 = y1;
            y1 = auxB;
        }
        mcd = a;
        auxA = x2;
        auxB = y2;
    }
}

// modular inverse using Extended Euclides Algorithm
long long signed int invMult(const long long signed int a, const long long signed int mod){
    long long signed int p = mod;
    long long signed int aa = a;
    long long signed int result;
    long long signed int q, r;
    long long signed int x1, x2, y1, y2;
    long long signed int auxA, auxB;
    if(mod != 0){
        x1 = 0;
        x2 = 1;
        y1 = 1;
        y2 = 0;

        while(p > 0){
            q = aa / p;
            r = aa % p;
            aa = p;
            p = r;
            auxA = x2 - (q*x1);
            auxB = y2 - (q*y1);
            x2 = x1;
            x1 = auxA;
            y2 = y1;
            y1 = auxB;
        }
	
        auxA = x2;
        auxB = y2;
	
	if(aa == 1){
	    do{
		auxA += mod;
	    }
	    while(auxA <= 0);
	    
	    result = auxA;
	}
	else{	    
	    std::cout<<"not inverse? mod = "<<mod<<std::endl;
	}
    }
    return result;
}

// return inverse multiplication
Int Int::operator/(const Int & a)const {
    assert (a.getValue() != 0);
    Int result;
    if(mod!=0){
        long long signed int a_val;
	long long signed int a_mod;
        a_val = a.val;
        a_mod = a.mod;
	long long signed int inv_a = invMult(a_val, a_mod);
        Int new_inv_a = Int(inv_a, a_mod);
        result = (*this) * new_inv_a;
    }
    else{
        result.mod = 0;
        result.val = val/a.val;
    }
    return result;
}

// return the binary representation of an integer
std::vector<int> intToBin(const long long signed int a){
    std::vector<int> result;
    std::stack<int> stac;
    long long signed int b = a;
    while(b > 0){
	stac.push(b%2);
	b = b/2;
    }
    while(!stac.empty()){
	result.push_back(stac.top());
	stac.pop();
    }
    return result;
}

// repeated squaring algorithm for power
Int Int::operator^(const long long signed int n){
    std::vector<int> bin_n = intToBin(n);
    int k = bin_n.size();
    std::vector<Int> b;
    b.push_back(*this);
    Int aux;
    
    for(int i = 1; i < k; i++){
	aux = b[i-1]*b[i-1];
	if(bin_n[i] == 1){	   
	    b.push_back(aux*(*this));
	}
	
	else{
	    b.push_back(aux);
	}
    }
    Int b1 = Int(b[k-1]);
    return b1;
}

// Chinesse Remainder Theorem
long long signed int garnerRemainder(long long signed int u[], long long signed int m[], long long signed int l){    
    long long signed int c[l];
    long long signed int result;
    long long signed int v;
    for(int i = 1; i<=l; i++){
	c[i] = 1;
	for(int j = 0; j<i; j++){
	    v = invMult(m[j], m[i]);
	    c[i] = (v*c[i]) % m[i];
	    c[i] = (c[i]>0) ? c[i] : c[i]+m[i];
	}
    }

    v = u[0];
    result = v;
   
    for(int i = 1; i<=l; i++){
	long long signed int prod = 1;
	v = ((u[i] - result)*c[i]) % m[i];
	v = (v>0) ? v : v + m[i];
	for(int j = 0; j<=(i-1); j++){
	    prod = prod * m[j];
	}
	result = result + (v * prod);
    }
    return result;
}

std::vector<int> prime_divisors(const long long signed int n){
    long long signed int m = n;
    std::vector<int> primes;
    for(int i = 2; i <= m; i++){
	bool T = true;
	while(m % i == 0){
	    if(T == true){
		primes.push_back(i);
	    }
	    else{
		if(primes.back() != i){
		    primes.push_back(i);
		}
	    }
	    T = false;
	    m = m/i;
	}
    }
    return primes;
}

// for printing in standard output
std::ostream & operator <<(std::ostream &os, const Int & a){
    os<<a.getValue();
    return os;
}
